# tests/ebm_export.R --------------------------------------------------------------
# THE EBM COMPARATOR, STEP 1: write the design tables the EBM is fitted on.
# (Author, 2026-10-28; decisions in docs/ebm_comparator_20261028.md.)
#
# The EBM gets EXACTLY the input of the `xgb_feat` booster: the 99 layer-1
# covariates built by `xgb_design_feat()`, with NA where a channel is unmeasured
# (the EBM handles missing values natively, as the booster does). The design is
# fold-dependent (`delta`, `lambda` are fitted), so, as `xgb_oof_perfold()` does,
# one table is built per held-out fold with that fold's out-of-fold priors:
#
#   design_fold<f>   all training stays, oof priors of fold f; the EBM for fold f
#                    trains on the rows outside f and predicts the rows in f
#   design_final     all training stays, the bundle's frozen final priors
#   design_eicu      the eICU cohort, frozen final priors (apply only)
#   design_test      the MIMIC-IV held-out test stays, frozen final priors (the
#                    EBM's single first look, author 2026-10-20; scored once by
#                    tests/ebm_score.R --test)
#   rows             stay order, outcome, fold, patient group of the training rows
#
# Column sets are asserted identical across every table. The final priors used
# here are asserted to be the bundle's. The tables are row-level derived data:
# they stay under out/cache/ebm/ and nothing row-level is printed (hard rule 1).
#
#   Rscript tests/ebm_export.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(targets); library(arrow); library(yaml); library(mgcv) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

OUT <- "out/cache/ebm"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cfg    <- tar_read(cfg)
tabs   <- tar_read(tabs)
folds  <- tar_read(folds)
priors <- tar_read(priors)
tr     <- tar_read(train_ids)
rc     <- yaml::read_yaml("config/internal.yml")
bundle <- load_bundle(rc$test_look$bundle, cfg = cfg, strict = TRUE, verbose = FALSE)

# The final priors of the targets store must be the bundle's: the final EBM is
# fitted on one and applied through the other.
fin <- function(p) lapply(cfg$signals, function(sg) priors_for(p, sg, "final", NA_integer_))
if (!isTRUE(all.equal(fin(priors), fin(bundle$priors), tolerance = 0))) {
  stop("the targets store's final priors differ from the bundle's", call. = FALSE)
}

ids  <- as.character(tr)
y    <- tabs$cohort$mortality[match(ids, as.character(tabs$cohort$stay_id))]
fold <- folds$fold[match(ids, as.character(folds$stay_id))]
grp  <- patient_group_of(tabs$cohort, cfg, tr)
if (anyNA(y) || anyNA(fold) || anyNA(grp)) stop("a training stay lacks an outcome, fold or patient group", call. = FALSE)

write_design <- function(X, name) {
  d <- as.data.frame(X, check.names = FALSE); rownames(d) <- NULL
  arrow::write_parquet(d, file.path(OUT, paste0(name, ".parquet")))
  c(rows = nrow(X), cols = ncol(X), na_share = round(mean(is.na(X)), 4))
}

ref <- NULL; info <- list()
for (f in sort(unique(fold))) {
  X <- xgb_design_feat(tabs, cfg, priors, tr, role = "oof", fold = f)
  if (!identical(rownames(X), ids)) stop("fold ", f, " design is not in training-stay order", call. = FALSE)
  if (is.null(ref)) ref <- colnames(X) else if (!identical(colnames(X), ref)) stop("fold ", f, " has a different column set", call. = FALSE)
  info[[paste0("design_fold", f)]] <- write_design(X, paste0("design_fold", f))
}
Xf <- xgb_design_feat(tabs, cfg, bundle$priors, tr, role = "final", feature_names = ref)
if (!identical(colnames(Xf), ref) || !identical(rownames(Xf), ids)) stop("final design is misaligned", call. = FALSE)
info$design_final <- write_design(Xf, "design_final")
arrow::write_parquet(data.frame(row = seq_along(ids), y = as.integer(y), fold = as.integer(fold),
                                group = as.character(grp)), file.path(OUT, "rows.parquet"))

# The MIMIC-IV held-out test stays: built here, scored once by ebm_score.R --test.
test_ids <- folds$stay_id[folds$split == "test"]
Xt <- xgb_design_feat(tabs, cfg, bundle$priors, test_ids, role = "final", feature_names = ref)
if (!identical(colnames(Xt), ref)) stop("test design is misaligned", call. = FALSE)
info$design_test <- write_design(Xt, "design_test")
saveRDS(as.character(test_ids), file.path(OUT, "test_ids.rds"))

# eICU, with the bundle's design snapshot and frozen final priors (apply only).
ec     <- yaml::read_yaml("config/external.yml")
cfg_e  <- bundle_cfg(bundle, ec$paths)
tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
e_ids  <- tabs_e$cohort$stay_id
Xe <- xgb_design_feat(tabs_e, cfg_e, bundle$priors, e_ids, role = "final", feature_names = ref)
if (!identical(colnames(Xe), ref)) stop("eICU design is misaligned", call. = FALSE)
info$design_eicu <- write_design(Xe, "design_eicu")
saveRDS(as.character(e_ids), file.path(OUT, "eicu_ids.rds"))
saveRDS(ids, file.path(OUT, "train_ids.rds"))

writeLines(ref, file.path(OUT, "feature_names.txt"))
summ <- do.call(rbind, lapply(names(info), function(k) data.frame(table = k, t(info[[k]]))))
write.csv(summ, file.path(OUT, "export_summary.csv"), row.names = FALSE)
cat(sprintf("training stays %d, deaths %d, folds %d; columns %d\n", length(ids), sum(y), length(unique(fold)), length(ref)))
print(summ, row.names = FALSE)
