# paper/make_prior_dominance.R -------------------------------------------------
# How often the Dirichlet prior outweighs a stay's own observed hours.
#
# Added 2026-10-25 (code review, D open question 1). Results and the
# likelihood-assumptions appendix quote "the prior outweighed a stay's own
# observed hours for 2.7% of measured channels in the MIMIC-IV training stays,
# mostly lactate (36%) and glucose (24%)", and no committed script computed it.
#
# THE QUANTITY. A channel-stay is flagged when the prior's weight in the
# posterior mean, alpha0 / (alpha0 + n), exceeds 0.5, i.e. n < alpha0: the same
# rule as the card's double dagger (run/patient_card.R, HEAVY_SHRINK_PRIOR_WEIGHT).
# alpha0 is the FINAL (full-train) concentration of the channel; n is the stay's
# observed hours. The population is every measured channel-stay of the MIMIC-IV
# training split (n_obs > 0), the rows the priors were fitted on.
#
#   Rscript paper/make_prior_dominance.R
#
# Reads the targets store (row-level, locally); writes one aggregate row per
# channel plus an overall row to paper/figs/coefs/prior_dominance.csv.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(targets))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

OUT <- "paper/figs/coefs/prior_dominance.csv"
THRESHOLD <- 0.5   # the card's rule; see run/patient_card.R

priors <- tar_read(priors); tabs <- tar_read(tabs); folds <- tar_read(folds)
sp <- priors$signal[priors$signal$role == "final", , drop = FALSE]
sp$alpha0 <- sp$alpha_low + sp$alpha_mid + sp$alpha_high
tr <- folds$stay_id[folds$split == "train"]
sf <- tabs$signal_features
sf <- sf[sf$stay_id %in% tr & sf$n_obs > 0, c("signal", "n_obs"), drop = FALSE]
a0 <- sp$alpha0[match(sf$signal, sp$signal)]
if (anyNA(a0)) stop("a measured channel has no final prior", call. = FALSE)
w <- a0 / (a0 + sf$n_obs)
flag <- w > THRESHOLD

per <- do.call(rbind, lapply(split(seq_len(nrow(sf)), sf$signal), function(i) data.frame(
  signal = sf$signal[i[1]], alpha0 = round(a0[i[1]], 4), n_measured = length(i),
  n_flagged = sum(flag[i]), share_flagged = round(mean(flag[i]), 4),
  stringsAsFactors = FALSE)))
per <- per[order(-per$share_flagged), ]
per$share_of_all_flags <- round(per$n_flagged / sum(flag), 4)
tot <- data.frame(signal = "ALL", alpha0 = NA_real_, n_measured = nrow(sf), n_flagged = sum(flag),
                  share_flagged = round(mean(flag), 4), share_of_all_flags = 1)
res <- rbind(per, tot)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(res, OUT, row.names = FALSE)
print(res, row.names = FALSE)
cat("wrote", OUT, "\n")
