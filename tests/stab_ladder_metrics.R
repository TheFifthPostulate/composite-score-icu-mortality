# tests/stab_ladder_metrics.R -------------------------------------------------------
# THE STABILITY LADDER: CONSUMER (author, 2026-10-30). Reads the store written by
# tests/stab_ladder_gen.R and fits nothing.
#
# QUESTION. Retrain each arm on bootstrap bags of the training set: does the same
# explanation come back for the same patient? Four arms, one step apart:
#   score -> jgam_gated      joint fitting
#   jgam_gated -> jgam_ungated   reading missingness as evidence
#   jgam_ungated -> shap     the model class (boosted trees read through SHAP)
# Every arm's explanation is an n x 31 matrix: 19 channels, then 12 paired
# interventions, each intervention whole and counted once.
#
# MASKS. An entry an arm's own rule makes structurally zero (masks/, written by
# the generator) is set to 0 before anything is ranked or marked: a channel is
# available where measured; an intervention, for the score and the gated GAM,
# where any of its paired channels is measured, and for the ungated GAM and SHAP
# always.
#
# THE BOLD RULE OVER 31 COLUMNS. A column is bold for a patient when its
# contribution is larger than that of at least 85% of the stays where that
# column is available, in the same replicate (each replicate its own reference,
# as a retrained model recomputes the card's). Counting STRICTLY smaller values
# departs from the card's rule (smaller or equal) only at exact ties, and the
# channels are continuous. The intervention columns are not: every unexposed
# stay carries the same intervention contribution in a GAM, a block of often
# more than 85% of stays. Under the card's rule that whole block would sit at
# the block's upper percentile and be bold; under this one it is bold only if
# 85% of stays lie below it. The size of the largest tied block per column is
# reported (tables/tied_blocks).
#
# PER PAIR OF REPLICATES, OVER PATIENTS:
#   jaccard    mean Jaccard index of the two bold sets, over patients with a
#              bold column in either replicate
#   identical  share of patients whose two bold sets are identical, both empty
#              included
#   retention  share of one replicate's bold marks also bold in the other,
#              averaged over the two directions
#   top1       masked leading-column agreement, tie-tolerant at delta = 0 (the
#              definition of attr_topk_agree, written out for 31 columns, which
#              the 30-column bitmask helper cannot hold)
#   *_ch       the same for the 19 channel columns alone
# PRIMARY (author, 2026-10-30): the 19-channel figures (jaccard_ch,
# identical_ch) carry the cross-arm comparison. In the GAM arms an intervention
# column is tied at the unexposed value for 82-99% of stays, so only exposed
# stays can mark it; SHAP's intervention values are never tied, so it marks
# about one more column per patient and its 31-column agreement mixes two
# marking rules. The 31-column figures are reported beside, with that reason.
# CONTRASTS:
#   bag        every pair of bags, seed 0: the primary contrast. The joint GAMs
#              and the score are deterministic given the bag, so this is the
#              same estimand for all four arms, and the pairs are the same bag
#              pairs, compared pair by pair across arms.
#   bag_seed   SHAP only: bag i at seed 0 against bag j at seed 1, the combined
#              bag-and-seed reading of the attribution store.
#   vs_full    each bag against bag 0, the model fitted on the whole training set.
# BY MISSINGNESS: jaccard_ch, identical_ch (primary), jaccard, identical and top1
# within patients grouped by how many channels they have unmeasured (0, 1-2, 3 or
# more). The first version computed only the 31-column figures by band (fixed
# 2026-10-30, before any band figure entered the paper).
#
# Aggregates only (hard rule 1): no row of the store reaches the console.
#
#   Rscript tests/stab_ladder_metrics.R out/runs/stabgen_<id> [--allow-partial]
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args  <- commandArgs(trailingOnly = TRUE)
STORE <- args[!startsWith(args, "--")][1]
if (is.na(STORE) || !dir.exists(STORE)) stop("give the store: out/runs/stabgen_<id>", call. = FALSE)
ALLOW_PARTIAL <- "--allow-partial" %in% args
ARMS  <- c("score", "jgam_gated", "jgam_ungated", "shap")
SITES <- c("test", "eicu")
STEPS <- list(c("jgam_gated", "score"), c("jgam_ungated", "jgam_gated"), c("shap", "jgam_ungated"), c("shap", "score"))
PCT_BOLD   <- 85
MISS_BANDS <- c(0, 1, 3, Inf)
MISS_LAB   <- c("0", "1-2", ">=3")

# --- the store and its guards ---------------------------------------------------
ID <- qs2::qs_read(file.path(STORE, "identity.qs2"))
COLS <- ID$columns; NC <- length(COLS)
cfg  <- tar_read(cfg); sigs <- as.character(unlist(cfg$signals))
if (!identical(COLS[seq_along(sigs)], sigs) || NC != 31L) stop("the store's columns are not the 19 channels and 12 interventions", call. = FALSE)
st <- read_manifest(STORE)$status
if (!identical(st, "complete") && !ALLOW_PARTIAL) stop("the store is not closed complete (", st, "); run the generator's --close, or pass --allow-partial", call. = FALSE)
cpath <- function(arm, site, b, s) file.path(STORE, "contrib", arm, site, sprintf("b%02d_s%d.qs2", b, s))
excluded <- function(arm, b) file.exists(file.path(STORE, "excluded", sprintf("%s_b%02d.txt", arm, b)))
bags_of <- function(arm, s = 0L) {
  f <- list.files(file.path(STORE, "contrib", arm, "eicu"), sprintf("_s%d[.]qs2$", s))
  b <- as.integer(sub("^b(\\d+)_.*$", "\\1", f))
  b[vapply(b, function(x) all(file.exists(vapply(SITES, function(si) cpath(arm, si, x, s), ""))) && !excluded(arm, x), TRUE)]
}
BAGS <- Reduce(intersect, lapply(ARMS, bags_of))
BAGS <- sort(BAGS[BAGS > 0L])
if (!all(vapply(ARMS, function(a) 0L %in% bags_of(a), TRUE))) stop("bag 0 is missing for an arm", call. = FALSE)
SHAP_S1 <- sort(intersect(bags_of("shap", 1L), BAGS))
cat(sprintf("store %s (%s): %d bags common to all four arms; SHAP seed-1 bags %d\n", basename(STORE), st, length(BAGS), length(SHAP_S1)))
if (length(BAGS) < 2L) stop("fewer than two common bags", call. = FALSE)

# --- encodings ------------------------------------------------------------------
POW <- 2^(seq_len(NC) - 1L)                     # 31 bits fit an R integer
CH_BITS <- sum(POW[seq_along(sigs)])
POP8 <- vapply(0:255, function(v) sum(bitwAnd(v, 2L^(0:7)) > 0), integer(1))
popcount <- function(v) POP8[v %% 256L + 1L] + POP8[(v %/% 256L) %% 256L + 1L] +
  POP8[(v %/% 65536L) %% 256L + 1L] + POP8[v %/% 16777216L + 1L]
stopifnot(popcount(c(0L, 1L, 255L, 256L, .Machine$integer.max, 2^19 - 1)) == c(0, 1, 8, 1, 31, 19))
pack <- function(B) as.integer(B %*% POW)       # logical n x 31 -> bitmask

#' One replicate: masked bold set, leading column and its tie set, and the size
#' of each column's largest block of tied values among available stays.
encode <- function(C, M) {
  if (!identical(dim(C), dim(M)) || !identical(rownames(C), rownames(M)) || !identical(colnames(C), COLS)) stop("a replicate is not in the mask's shape and order", call. = FALSE)
  C[!M] <- 0
  B <- matrix(FALSE, nrow(C), NC); tied <- numeric(NC)
  for (j in seq_len(NC)) {
    w <- which(M[, j]); if (!length(w)) next
    x <- C[w, j]; ref <- sort(x)
    below <- findInterval(x, ref, left.open = TRUE)          # strictly smaller values
    B[w, j] <- round(100 * below / length(ref)) >= PCT_BOLD
    tied[j] <- max(tabulate(match(x, unique(x)))) / length(x)
  }
  A <- abs(C); top <- max.col(A, ties.method = "first"); mx <- A[cbind(seq_len(nrow(A)), top)]
  list(bold = pack(B), top = top, tie = pack(A == mx), tied = tied)
}
#' Agreement of two encodings, per patient.
agree <- function(a, b, bits = .Machine$integer.max) {
  A <- bitwAnd(a$bold, bits); Bb <- bitwAnd(b$bold, bits)
  inter <- popcount(bitwAnd(A, Bb)); uni <- popcount(bitwOr(A, Bb))
  list(inter = inter, uni = uni, na = popcount(A), nb = popcount(Bb), same = A == Bb)
}
top1_ok <- function(a, b) {
  # tie-tolerant at delta = 0: b's leading column is tied with a's maximum, and vice versa
  bitwAnd(a$tie, as.integer(POW[b$top])) != 0L & bitwAnd(b$tie, as.integer(POW[a$top])) != 0L
}
summ_pair <- function(g, ok, band) {
  j <- ifelse(g$uni > 0, g$inter / pmax(g$uni, 1L), NA_real_)
  out <- c(jaccard = mean(j, na.rm = TRUE), identical = mean(g$same),
           retention = mean(c(sum(g$inter) / max(sum(g$na), 1), sum(g$inter) / max(sum(g$nb), 1))),
           bold_per_patient = mean(c(g$na, g$nb)), top1 = mean(ok))
  byb <- do.call(rbind, lapply(MISS_LAB, function(l) { w <- band == l
    data.frame(band = l, share_patients = mean(w), jaccard = if (any(w & g$uni > 0)) mean(j[w], na.rm = TRUE) else NA_real_,
               identical = mean(g$same[w]), top1 = mean(ok[w])) }))
  list(out = out, byb = byb)
}

RAW <- list(); BAND <- list(); TIED <- list()
t_all <- start_timer()
for (site in SITES) {
  masks <- lapply(stats::setNames(ARMS, ARMS), function(a) qs2::qs_read(file.path(STORE, "masks", sprintf("%s_%s.qs2", a, site))))
  ids <- rownames(masks$score)
  if (!identical(attr_key_hash(ids), ID[[paste0(site, "_ids")]])) stop(site, ": the masks' stays are not the store's", call. = FALSE)
  meas <- masks$score[, sigs]
  for (a in ARMS) if (!identical(masks[[a]][, sigs], meas)) stop(site, ": ", a, "'s channel mask differs from the measured mask", call. = FALSE)
  band <- as.character(cut(rowSums(!meas), MISS_BANDS, MISS_LAB, right = FALSE))
  if (anyNA(band)) stop("a patient falls outside the missingness bands", call. = FALSE)
  for (a in ARMS) {
    ENC <- list()
    keys <- c(sprintf("b%02d_s0", c(0L, BAGS)), if (a == "shap") sprintf("b%02d_s1", SHAP_S1))
    for (k in keys) {
      b <- as.integer(substr(k, 2, 3)); s <- as.integer(substr(k, 6, 6))
      ENC[[k]] <- encode(qs2::qs_read(cpath(a, site, b, s)), masks[[a]])
    }
    TIED[[length(TIED) + 1L]] <- data.frame(site = site, arm = a, column = COLS,
                                            kind = ifelse(COLS %in% sigs, "channel", "intervention"),
                                            largest_tied_share = ENC[["b00_s0"]]$tied, stringsAsFactors = FALSE)
    pairs <- rbind(
      data.frame(contrast = "bag", t(utils::combn(BAGS, 2L)), s_i = 0L, s_j = 0L),
      data.frame(contrast = "vs_full", X1 = 0L, X2 = BAGS, s_i = 0L, s_j = 0L),
      if (a == "shap") { g <- expand.grid(X1 = BAGS, X2 = SHAP_S1); g <- g[g$X1 != g$X2, ]
        data.frame(contrast = "bag_seed", X1 = g$X1, X2 = g$X2, s_i = 0L, s_j = 1L) })
    for (q in seq_len(nrow(pairs))) {
      ea <- ENC[[sprintf("b%02d_s%d", pairs$X1[q], pairs$s_i[q])]]; eb <- ENC[[sprintf("b%02d_s%d", pairs$X2[q], pairs$s_j[q])]]
      ok <- top1_ok(ea, eb)
      all31 <- summ_pair(agree(ea, eb), ok, band)
      chs <- summ_pair(agree(ea, eb, CH_BITS), ok, band); ch <- chs$out
      byb <- all31$byb; byb$jaccard_ch <- chs$byb$jaccard; byb$identical_ch <- chs$byb$identical
      RAW[[length(RAW) + 1L]] <- data.frame(site = site, arm = a, contrast = pairs$contrast[q], bag_i = pairs$X1[q], bag_j = pairs$X2[q],
        jaccard = all31$out[["jaccard"]], identical = all31$out[["identical"]], retention = all31$out[["retention"]],
        bold_per_patient = all31$out[["bold_per_patient"]], top1 = all31$out[["top1"]],
        jaccard_ch = ch[["jaccard"]], identical_ch = ch[["identical"]], bold_per_patient_ch = ch[["bold_per_patient"]], stringsAsFactors = FALSE)
      BAND[[length(BAND) + 1L]] <- cbind(data.frame(site = site, arm = a, contrast = pairs$contrast[q], bag_i = pairs$X1[q], bag_j = pairs$X2[q]), byb)
    }
    cat(sprintf("  %s %s: %d replicates, %d pairs (%.1f min so far)\n", site, a, length(ENC), nrow(pairs), t_all()$elapsed_sec / 60))
    rm(ENC); invisible(gc(verbose = FALSE))
  }
}
raw <- do.call(rbind, RAW); bnd <- do.call(rbind, BAND); tied <- do.call(rbind, TIED)

# --- summaries --------------------------------------------------------------------
MET <- c("jaccard", "identical", "retention", "bold_per_patient", "top1", "jaccard_ch", "identical_ch", "bold_per_patient_ch")
q3 <- function(v) c(median = stats::median(v), p025 = unname(stats::quantile(v, 0.025)), p975 = unname(stats::quantile(v, 0.975)))
summ <- do.call(rbind, lapply(split(raw, paste(raw$site, raw$arm, raw$contrast)), function(d) do.call(rbind, lapply(MET, function(m)
  data.frame(site = d$site[1], arm = d$arm[1], contrast = d$contrast[1], metric = m, n_pairs = nrow(d), t(q3(d[[m]])))))))
# THE LADDER, PAIRED: the same bag pairs in every arm, difference pair by pair.
pk <- function(d) paste(d$site, d$bag_i, d$bag_j)
bagp <- raw[raw$contrast == "bag", ]
paired <- do.call(rbind, lapply(STEPS, function(stp) {
  A <- bagp[bagp$arm == stp[1], ]; B <- bagp[bagp$arm == stp[2], ]
  m <- merge(cbind(A, key = pk(A)), cbind(B, key = pk(B)), by = "key", suffixes = c("_a", "_b"))
  if (nrow(m) != nrow(A) || nrow(m) != nrow(B)) stop("the arms' bag pairs differ: ", paste(stp, collapse = " vs "), call. = FALSE)
  do.call(rbind, lapply(split(m, m$site_a), function(d) do.call(rbind, lapply(c("jaccard", "identical", "top1", "jaccard_ch", "identical_ch"), function(mt) {
    v <- d[[paste0(mt, "_a")]] - d[[paste0(mt, "_b")]]
    data.frame(site = d$site_a[1], step = paste(stp[1], "minus", stp[2]), metric = mt, n_bag_pairs = length(v), t(q3(v)), share_a_higher = mean(v > 0))
  }))))
}))
band_s <- stats::aggregate(cbind(share_patients, jaccard_ch, identical_ch, jaccard, identical, top1) ~ site + arm + contrast + band, data = bnd,
                           FUN = function(v) stats::median(v, na.rm = TRUE), na.action = stats::na.pass)
tied_s <- stats::aggregate(largest_tied_share ~ site + arm + kind, data = tied, FUN = function(v) c(median = stats::median(v), max = max(v)))

run <- new_run("stabmetrics", cfg, note = sprintf("stability ladder agreement over 31 columns, store %s; fits nothing", basename(STORE)))
save_table(run, summ, "agreement")
save_table(run, raw, "agreement_pairs")
save_table(run, paired, "ladder_paired")
save_table(run, band_s, "agreement_by_missingness")
save_table(run, tied, "tied_blocks")
finalize_run(run, extra = list(store = STORE, store_status = st, n_bags = length(BAGS), n_shap_seed1 = length(SHAP_S1), pct_bold = PCT_BOLD,
                                bold_rule = "strictly smaller share >= 85% of available stays, own reference",
                                primary = "19 channel columns (jaccard_ch, identical_ch); 31 columns secondary"))

options(width = 220)
r3 <- function(d) { d[] <- lapply(d, function(v) if (is.numeric(v)) round(v, 3) else v); d }
for (site in SITES) {
  cat(sprintf("\n== %s: agreement between bags (median [2.5th, 97.5th] over pairs)\n", site))
  w <- summ[summ$site == site & summ$metric %in% c("jaccard", "identical", "top1", "bold_per_patient", "jaccard_ch", "identical_ch"), ]
  w$value <- sprintf("%.3f [%.3f, %.3f]", w$median, w$p025, w$p975)
  wide <- stats::reshape(w[, c("arm", "contrast", "n_pairs", "metric", "value")], idvar = c("arm", "contrast", "n_pairs"), timevar = "metric", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  wide$arm <- factor(wide$arm, ARMS)
  print(wide[order(wide$contrast, wide$arm), ], row.names = FALSE)
  cat(sprintf("\n== %s: the ladder, paired on the same bag pairs\n", site))
  print(r3(paired[paired$site == site, c("step", "metric", "n_bag_pairs", "median", "p025", "p975", "share_a_higher")]), row.names = FALSE)
  cat(sprintf("\n== %s: by number of unmeasured channels (bag pairs, median over pairs)\n", site))
  bs <- band_s[band_s$site == site & band_s$contrast == "bag", ]; bs$arm <- factor(bs$arm, ARMS)
  print(r3(bs[order(bs$band, bs$arm), c("band", "arm", "share_patients", "jaccard_ch", "identical_ch", "jaccard", "identical", "top1")]), row.names = FALSE)
}
cat("\n== largest block of tied values per column, among available stays (bag 0): median and max over columns\n")
print(r3(do.call(data.frame, tied_s)), row.names = FALSE)
cat("\nrun:", run$path, "\n")
