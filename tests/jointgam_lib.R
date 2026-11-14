# tests/jointgam_lib.R -------------------------------------------------------------
# The joint GAM's machinery, moved out of tests/jointgam.R (2026-10-30) so the
# stability-ladder generator (tests/stab_ladder_gen.R) fits bags through the same
# code; fit_joint() then rebuilt on the discretisable "gz" smooth class (below).
# Header, construction and checks: tests/jointgam.R.
#
# The caller defines, before sourcing: cfg, bs_ (bam_settings(cfg)), feat (the
# design's column names), UNGATED, VARIANT, and, before calling fit_joint(),
# y_all (the outcome of every design row). Sourced into its own environment
# (sys.source) when one process needs both variants.
# ------------------------------------------------------------------------------

# --- the block map: every design column, its channel or intervention, its term ----
paired_of <- function(iv) cfg$signals[vapply(cfg$signals, function(s) iv %in% interventions_of(s, cfg), TRUE)]
ivs <- unique(unlist(lapply(cfg$signals, interventions_of, cfg = cfg)))

term_of <- function(sg, var) {
  tt <- attr(stats::terms(build_formula(sg, "full", cfg)), "term.labels")
  hit <- tt[vapply(tt, function(t) {
    if (startsWith(t, "s(")) identical(.call_args(t)$vars, var) else identical(t, var)
  }, TRUE)]
  if (length(hit) != 1L) stop(sprintf("%s: %d terms for variable %s", sg, length(hit), var), call. = FALSE)
  hit
}

blocks <- lapply(feat, function(col) {
  sg <- cfg$signals[startsWith(col, paste0(cfg$signals, "__"))]
  if (length(sg) == 1L) {
    var <- sub(paste0("^", sg, "__"), "", col)
    return(list(col = col, owner = sg, kind = "measurement", gate_of = sg, term = term_of(sg, var), var = var))
  }
  iv <- ivs[startsWith(col, paste0(ivs, "__"))]
  if (length(iv) != 1L) stop("design column ", col, " belongs to no channel or intervention", call. = FALSE)
  ps <- paired_of(iv)
  tms <- unique(vapply(ps, term_of, "", var = col))
  if (length(tms) != 1L) stop(col, ": its paired channels declare different terms", call. = FALSE)
  list(col = col, owner = iv, kind = "intervention", gate_of = ps, term = tms, var = col)
})
names(blocks) <- feat
# Every term of every `full` formula must be a design column, and nothing else.
n_terms <- sum(vapply(cfg$signals, function(sg) length(attr(stats::terms(build_formula(sg, "full", cfg)), "term.labels")), 1L))
n_cover <- sum(vapply(blocks, function(b) length(b$gate_of), 1L))
if (n_cover != n_terms) stop(sprintf("design covers %d channel terms, the formulas hold %d", n_cover, n_terms), call. = FALSE)
if (any(grepl("^ti\\(", unlist(lapply(blocks, `[[`, "term"))))) stop("an interaction term reached the additive joint GAM", call. = FALSE)

# Gates from the missingness pattern, asserted all-or-none per channel.
gates_of <- function(X, name) {
  g <- list()
  for (sg in cfg$signals) {
    cs <- feat[vapply(blocks, function(b) b$kind == "measurement" && b$owner == sg, TRUE)]
    na <- is.na(as.matrix(X[cs]))
    if (any(rowSums(na) %% length(cs) != 0)) stop(name, ": ", sg, " has partly missing measurement columns", call. = FALSE)
    g[[sg]] <- as.numeric(!na[, 1])
  }
  iv_cols <- feat[vapply(blocks, function(b) b$kind == "intervention", TRUE)]
  if (anyNA(X[iv_cols])) stop(name, ": an intervention column has missing values", call. = FALSE)
  for (b in blocks) if (b$kind == "intervention") {
    g[[b$col]] <- if (UNGATED) rep(1, nrow(X)) else as.numeric(Reduce(`|`, lapply(b$gate_of, function(s) g[[s]] == 1)))
  }
  g
}
# The ungated arm's missingness constants, one per channel (see header).
MISS <- if (UNGATED) paste0(cfg$signals, "__unmeasured") else character(0)
miss_signal <- function(nm) sub("__unmeasured$", "", nm)
gate_vec <- function(g, b) if (b$kind == "measurement") g[[b$owner]] else g[[b$col]]

# --- building blocks ---------------------------------------------------------------
#' One block's basis, built on its gated fitting rows. Returns the smooth object
#' (or the binary centre) with no row-level matrix kept.
make_block <- function(b, x, gate) {
  xg <- x[gate == 1]
  if (startsWith(b$term, "s(")) {
    spec <- eval(str2lang(sub(paste0("s(", b$var), "s(v", b$term, fixed = TRUE)))
    # The centring constraint is absorbed HERE, explicitly, as the column means
    # over the gated fitting rows. mgcv's `absorb.cons = TRUE` imposes the same
    # constraint, but only to about 1e-8 (lactate pi_minus: 1.45e-8 over the
    # rows, measured 2026-10-30), which failed the gate check's tolerance. The
    # centre decides what "not measured" means here, so it is made exact.
    sm <- mgcv::smoothCon(spec, data = data.frame(v = xg), absorb.cons = FALSE)[[1]]
    if (length(sm$S) != 1L) stop(b$col, ": expected one penalty", call. = FALSE)
    cm <- colMeans(mgcv::PredictMat(sm, data.frame(v = xg)))
    Z  <- qr.Q(qr(matrix(cm, ncol = 1)), complete = TRUE)[, -1, drop = FALSE]
    S  <- crossprod(Z, sm$S[[1]] %*% Z)
    sm$X <- NULL
    list(type = "smooth", sm = sm, Z = Z, S = (S + t(S)) / 2, fill = stats::median(xg))
  } else {
    # Only 0/1 indicators may enter as linear terms: the design stores factors
    # as integer codes (.feat_numeric), which a linear term would misread.
    if (!grepl("__present_at_admission$", b$var) || !all(x %in% c(0, 1))) {
      stop(b$col, ": a non-smooth term that is not a 0/1 admission indicator", call. = FALSE)
    }
    list(type = "binary", centre = mean(xg))
  }
}
#' A block's columns for any rows: the basis times the gate.
block_matrix <- function(bb, x, gate) {
  if (bb$type == "smooth") {
    x[gate == 0] <- bb$fill
    M <- mgcv::PredictMat(bb$sm, data.frame(v = x)) %*% bb$Z
  } else {
    x[gate == 0] <- bb$centre
    M <- matrix(x - bb$centre, ncol = 1)
  }
  M * gate
}

# --- the gated smooth class ----------------------------------------------------------
# DISCRETISED LIKE THE CHANNEL MODELS (2026-10-30). Each smooth block enters bam as
# s(V, by = G, bs = "gz", xt = list(sm, Z)): a smooth class whose basis is the
# block's own (its base smooth times its explicit centring Z, both built on the
# gated fitting rows by make_block()) and whose penalty is t(Z) S Z. The numeric
# `by` multiplies each row by its gate, so the term is exactly zero off the gate,
# and the class declares no further constraint (C has no rows), so no free level
# can absorb missingness. Because the terms are s() terms, bam discretises them
# as it does the channel models' (bam.discrete in config/config.yml); the earlier
# paraPen construction had no s() terms and bam ignored `discrete`.
# Prediction never uses the discretised representation: predict_joint() and
# contrib_joint() evaluate the undiscretised basis, whose centring is exact.
smooth.construct.gz.smooth.spec <- function(object, data, knots) {
  sm <- object$xt$sm; Z <- object$xt$Z
  X <- mgcv::PredictMat(sm, stats::setNames(data.frame(data[[object$term]]), sm$term)) %*% Z
  S <- crossprod(Z, sm$S[[1]] %*% Z); S <- (S + t(S)) / 2
  object$X <- X; object$S <- list(S)
  object$rank <- qr(S)$rank; object$null.space.dim <- ncol(X) - object$rank
  object$C <- matrix(0, 0, ncol(X)); object$bs.dim <- ncol(X); object$df <- ncol(X)
  object$base <- sm; object$Z <- Z
  class(object) <- "gz.smooth"
  object
}
Predict.matrix.gz.smooth <- function(object, data) {
  mgcv::PredictMat(object$base, stats::setNames(data.frame(data[[object$term]]), object$base$term)) %*% object$Z
}
registerS3method("smooth.construct", "gz.smooth.spec", smooth.construct.gz.smooth.spec, envir = asNamespace("mgcv"))
registerS3method("Predict.matrix", "gz.smooth", Predict.matrix.gz.smooth, envir = asNamespace("mgcv"))

#' Fit one joint GAM on rows `fit` of design X. Returns the stripped fit.
#' Coefficients come in bam's order: the intercept, the parametric blocks (the
#' admission indicators and, ungated, the unmeasured constants) in `use` order,
#' then the smooth blocks in `use` order; `blocks` of the fit is stored in that
#' order, so predict_joint() reads them by position.
fit_joint <- function(X, g, fit, role, fold_id, use = c(names(blocks), MISS)) {
  t0 <- proc.time()[["elapsed"]]
  y_fit <- y_all[fit]
  dat <- list(y = y_fit); par <- list(); smo <- list()
  for (nm in use) {
    if (nm %in% MISS) {
      # Ungated arm: 1[channel unmeasured], never centred, so its coefficient
      # reads directly as the shift for "not measured". The 19 columns enter
      # together as one matrix term with ONE shared ridge penalty, its weight
      # estimated with the smoothing parameters (author, 2026-10-30): unpenalised,
      # the 9 training stays missing SpO2 alone, all survivors, separated and
      # sent the SpO2 and heart-rate constants to about -26 and +26.
      u <- 1 - g[[miss_signal(nm)]][fit]
      if (sum(u) == 0) stop(nm, ": the channel is measured in every fitting stay", call. = FALSE)
      par[[nm]] <- list(type = "missing", name = nm, signal = miss_signal(nm), ncol = 1L, n_unmeasured = sum(u))
      dat$MISSM <- cbind(dat$MISSM, u)
      next
    }
    if (length(dat$MISSM)) stop("the missingness constants must come after every other block in `use`", call. = FALSE)
    b <- blocks[[nm]]; gt <- gate_vec(g, b)[fit]
    if (sum(gt) < 50) stop(nm, ": fewer than 50 gated fitting stays", call. = FALSE)
    x  <- X[[b$col]][fit]
    bb <- make_block(b, x, gt)
    M  <- block_matrix(bb, x, gt)
    # gate check: exactly 0 off the gate, column means 0 over the gated rows
    if (any(M[gt == 0, ] != 0)) stop(nm, ": block is not zero on ungated rows", call. = FALSE)
    if (max(abs(colMeans(M[gt == 1, , drop = FALSE]))) > 1e-8) stop(nm, ": block is not centred over its gated rows", call. = FALSE)
    if (bb$type == "binary") {
      par[[nm]] <- c(bb, list(name = nm, ncol = 1L))
      dat[[paste0("P", length(par))]] <- as.numeric(M)
    } else {
      smo[[nm]] <- c(bb, list(name = nm, ncol = ncol(M)))
      i <- length(smo); x[gt == 0] <- bb$fill
      dat[[paste0("V", i)]] <- x; dat[[paste0("G", i)]] <- gt
    }
  }
  n_bin <- sum(vapply(par, function(z) z$type == "binary", TRUE))
  if (n_bin && length(dat$MISSM)) {
    if (any(vapply(par[seq_len(n_bin)], function(z) z$type != "binary", TRUE))) stop("binary blocks must precede the missingness constants", call. = FALSE)
  }
  pp <- if (length(dat$MISSM)) list(MISSM = list(diag(ncol(dat$MISSM)))) else NULL
  rhs <- c(if (n_bin) paste0("P", seq_len(n_bin)), if (length(dat$MISSM)) "MISSM",
           if (length(smo)) sprintf("s(V%d, by = G%d, bs = \"gz\", xt = list(sm = smo[[%d]]$sm, Z = smo[[%d]]$Z))",
                                    seq_along(smo), seq_along(smo), seq_along(smo), seq_along(smo)))
  f <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")), env = environment())
  b <- mgcv::bam(f, data = dat, family = bs_$family, method = bs_$method, discrete = bs_$discrete,
                 nthreads = bs_$nthreads, select = bs_$select, gamma = bs_$gamma, gc.level = bs_$gc_level,
                 paraPen = pp, na.action = stats::na.fail)
  diag <- list(signal = "joint", model = VARIANT, role = role, fold = fold_id, converged = .bam_converged(b))
  .accept_outcome_fit(b, diag, cfg)
  built <- c(par, smo)
  cf <- stats::coef(b)
  if (length(cf) != 1L + sum(vapply(built, `[[`, 1L, "ncol"))) stop("coefficient count does not match the blocks", call. = FALSE)
  if (length(b$smooth) != length(smo) ||
      any(vapply(b$smooth, function(s) as.numeric(s$last.para - s$first.para + 1), 1) != vapply(smo, `[[`, 1L, "ncol"))) {
    stop("bam's smooth terms do not match the smooth blocks", call. = FALSE)
  }
  # The missingness constants with their standard errors (model parameters only).
  mk <- which(vapply(par, function(z) z$type == "missing", TRUE))
  miss <- if (length(mk)) data.frame(signal = vapply(par[mk], `[[`, "", "signal"),
                                     n_unmeasured = vapply(par[mk], `[[`, 1, "n_unmeasured"),
                                     coef = unname(cf[1L + mk]), se = sqrt(diag(b$Vp))[1L + mk],
                                     stringsAsFactors = FALSE, row.names = NULL) else NULL
  fit_obj <- list(blocks = built, coef = cf, edf = sum(b$edf), sp = b$sp, converged = diag$converged,
                  discrete_effective = isTRUE(bs_$discrete), iter = b$iter, criterion = b$gcv.ubre, miss = miss,
                  miss_sp = if (!is.null(pp)) unname(b$sp[1]) else NA_real_,   # paraPen weights come first in b$sp
                  Vp = unname(b$Vp),   # Bayesian covariance, coefficient order; bands for the atlas (2026-10-30)
                  n = length(fit), deaths = sum(y_fit), seconds = proc.time()[["elapsed"]] - t0)
  # The stripped fit must reproduce bam's own undiscretised prediction on the
  # fitting rows: this checks the coefficient order and the gz basis.
  eta_b  <- as.numeric(stats::predict(b, newdata = dat, type = "link", discrete = FALSE))
  eta_in <- predict_joint(fit_obj, X, g, fit)
  if (max(abs(eta_in - eta_b)) > 1e-8) stop("manual linear predictor differs from bam's (", format(max(abs(eta_in - eta_b))), ")", call. = FALSE)
  fit_obj
}

#' Linear predictor of a stripped fit for rows `idx` of design X with gates g.
predict_joint <- function(fo, X, g, idx) {
  eta <- rep(unname(fo$coef[1]), length(idx)); at <- 1L
  for (key in names(fo$blocks)) {
    bb <- fo$blocks[[key]]
    if (bb$type == "missing") {
      M <- matrix(1 - g[[bb$signal]][idx], ncol = 1)
    } else {
      b <- blocks[[bb$name]]
      M <- block_matrix(bb, X[[b$col]][idx], gate_vec(g, b)[idx])
    }
    eta <- eta + as.numeric(M %*% fo$coef[at + seq_len(bb$ncol)]); at <- at + bb$ncol
  }
  if (anyNA(eta)) stop("a joint GAM prediction is NA", call. = FALSE)
  eta
}

#' Per-column contributions of a stripped fit for rows `idx`: a length(idx) x
#' length(cols) matrix. Every block adds to the column of its owner: a channel's
#' measurement terms (and, ungated, its unmeasured constant) to the channel, an
#' intervention's terms to the intervention, once. The intercept belongs to no
#' column and is returned as attr "intercept"; the row sums plus the intercept
#' are asserted equal to predict_joint().
contrib_joint <- function(fo, X, g, idx, cols) {
  C <- matrix(0, length(idx), length(cols), dimnames = list(NULL, cols)); at <- 1L
  for (key in names(fo$blocks)) {
    bb <- fo$blocks[[key]]
    if (bb$type == "missing") {
      M <- matrix(1 - g[[bb$signal]][idx], ncol = 1); own <- bb$signal
    } else {
      b <- blocks[[bb$name]]
      M <- block_matrix(bb, X[[b$col]][idx], gate_vec(g, b)[idx]); own <- b$owner
    }
    if (!own %in% cols) stop("contrib_joint: block owner ", own, " is not a column", call. = FALSE)
    C[, own] <- C[, own] + as.numeric(M %*% fo$coef[at + seq_len(bb$ncol)]); at <- at + bb$ncol
  }
  ic <- unname(fo$coef[1])
  d <- max(abs(rowSums(C) + ic - predict_joint(fo, X, g, idx)))
  if (d > 1e-8) stop("contrib_joint: contributions do not sum to the linear predictor (", format(d), ")", call. = FALSE)
  attr(C, "intercept") <- ic
  C
}
