# tests/lambda_binom_check.R ---------------------------------------------------
# Residual check of the binomial intensity models (Equation eq:lambda_binom) for
# the vasopressor and inotrope agent counts, the one random-variable model whose
# assumption had no check in the diagnostics (results_appendix.tex, assumption 6;
# author, 2026-11-01).
#
# Under the model, the agent count A of an exposed stay is Binomial(M, p(D)) with
# logit p(D) = c0 + c1 log D, and the covariate is the Pearson residual
# (A - M p) / sqrt(M p (1 - p)) (lambda_value() in R/04b_conditional.R). Checks,
# per intervention and site, on the exposed stays:
#   - mean observed and mean expected count, and the mean Pearson residual
#     (zero at MIMIC-IV by the score equations of the logistic fit; at eICU the
#     size of a conditional-mean shift);
#   - dispersion, the mean squared Pearson residual (one under the model);
#   - the observed and expected shares of stays at each count 1..M.
# Parameters: the FINAL priors of the frozen bundle named in config/external.yml,
# applied to the MIMIC-IV training stays and to eICU. Application, not fitting
# (hard rule 9). Before any number is written, the MIMIC-IV mean raw residual
# must be zero to 1e-6, which confirms the parameters and the rows match the fit.
# Aggregates only (hard rule 1).
#
#   Rscript tests/lambda_binom_check.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(arrow); library(yaml); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

rc        <- yaml::read_yaml("config/internal.yml")
cfg_local <- load_config(rc$config)
xcfg      <- yaml::read_yaml("config/external.yml")
bpath     <- cfg_req(xcfg, "bundle")
bundle    <- load_bundle(bpath, cfg = cfg_local, strict = TRUE)
cfg_m     <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
tabs_m    <- load_tables(cfg_m$paths, cfg_m, site = "mimic", verbose = FALSE)
folds     <- assign_folds(tabs_m$cohort, cfg_m)
train_ids <- folds$stay_id[folds$split == "train"]
cfg_e     <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
tabs_e    <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)

ip <- bundle$priors$intervention
ip <- ip[ip$role == "final" & ip$family == "binomial", , drop = FALSE]
if (!nrow(ip)) stop("no binomial intensity models in the bundle's final priors", call. = FALSE)

rows <- list(); lev <- list(); dq <- list()
for (i in seq_len(nrow(ip))) {
  par <- as.list(ip[i, , drop = FALSE]); iv <- par$intervention
  sp  <- lambda_spec_of(iv, cfg_m)
  if (is.null(sp) || !identical(sp$family, "binomial") || !identical(as.integer(sp$M), as.integer(par$M)))
    stop(iv, ": bundle and design disagree on the binomial model", call. = FALSE)
  for (site in c("mimic_train", "eicu")) {
    z <- if (site == "mimic_train") tabs_m$intervention_features[tabs_m$intervention_features$intervention == iv &
                                      as.character(tabs_m$intervention_features$stay_id) %in% as.character(train_ids), ]
         else tabs_e$intervention_features[tabs_e$intervention_features$intervention == iv, ]
    D <- z[[sp$exposure]] * sp$exposure_scale; A <- z[[sp$accumulation]]
    ex <- is.finite(D) & D > 0 & is.finite(A)
    d <- D[ex]; a <- pmin(A[ex], par$M)
    pe <- stats::plogis(par$c0 + par$c1 * log(d)); mu <- par$M * pe
    r  <- lambda_value(A[ex], d, par)
    if (site == "mimic_train" && abs(mean(a - mu)) > 1e-6)
      stop(sprintf("%s: MIMIC-IV mean raw residual %.2e, expected 0; parameters and rows do not match the fit", iv, mean(a - mu)), call. = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(intervention = iv, site = site, M = par$M, n_exposed = sum(ex),
      n_above_pool = sum(A[ex] > par$M), mean_observed = mean(a), mean_expected = mean(mu),
      mean_resid = mean(r), dispersion = mean(r^2), sd_resid = stats::sd(r), stringsAsFactors = FALSE)
    # Conditional mean across exposure: mean observed and expected count within quartiles of D (cut at MIMIC-IV
    # training quartiles of the same intervention, so the two sites share bins).
    if (site == "mimic_train") qd <- unique(stats::quantile(d, c(0, 0.25, 0.5, 0.75, 1)))
    bin <- cut(d, breaks = unique(c(-Inf, qd[-c(1, length(qd))], Inf)), labels = FALSE)
    for (b in sort(unique(bin))) {
      w <- bin == b
      dq[[length(dq) + 1L]] <- data.frame(intervention = iv, site = site, exposure_bin = b,
        hours_lo = min(d[w]), hours_hi = max(d[w]), n = sum(w), mean_observed = mean(a[w]), mean_expected = mean(mu[w]),
        mean_resid = mean(r[w]), stringsAsFactors = FALSE)
    }
    exp_share <- vapply(0:par$M, function(k) mean(stats::dbinom(k, par$M, pe)), 0)
    lev[[length(lev) + 1L]] <- data.frame(intervention = iv, site = site, count = 0:par$M,
      observed_share = vapply(0:par$M, function(k) mean(a == k), 0), expected_share = exp_share,
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows); lv <- do.call(rbind, lev); dqt <- do.call(rbind, dq)

run <- new_run("lambdabinom", cfg_local, note = "residual check of the binomial agent-count models, frozen final priors")
save_table(run, tab, "lambda_binom_check")
save_table(run, lv, "lambda_binom_levels")
save_table(run, dqt, "lambda_binom_by_exposure")
finalize_run(run, extra = list(bundle = bpath))
options(width = 200)
print(transform(tab, mean_observed = round(mean_observed, 3), mean_expected = round(mean_expected, 3),
                mean_resid = round(mean_resid, 3), dispersion = round(dispersion, 3), sd_resid = round(sd_resid, 3)), row.names = FALSE)
print(transform(lv, observed_share = round(observed_share, 4), expected_share = round(expected_share, 4)), row.names = FALSE)
print(transform(dqt, mean_observed = round(mean_observed, 3), mean_expected = round(mean_expected, 3), mean_resid = round(mean_resid, 3)), row.names = FALSE)
cat("run:", run$path, "\n")
