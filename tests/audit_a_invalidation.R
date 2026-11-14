# tests/audit_a_invalidation.R -----------------------------------------------
# AUDIT A of docs/v2_audit_plan_20260903.md section 2, and invariance item 2:
# THE INPUT-INVALIDATION MATRIX. The highest-value item in the plan, because it
# closes the class that has actually bitten.
#
# THE FAILURE IT EXISTS FOR. `mimic_files` was a `format = "file"` target that
# nothing depended on. Re-extracting a parquet therefore rebuilt `mimic_files`
# and NOTHING ELSE: `tabs` was served from a cache that predated the new
# columns, and the pipeline died a hundred lines away validating a stale table.
# One line of this matrix -- "touching a parquet rebuilds 1 target" -- would
# have said so immediately. It also catches the REVERSE failure, a config edit
# that fails to invalidate the fits it should, which nothing currently detects
# at all.
#
# TWO ARMS, AND THEY ANSWER DIFFERENT QUESTIONS.
#
#   STATIC     For each tracked input, the set of targets REACHABLE from its
#              file target in `tar_network()`. This is what the graph says
#              would rebuild. It needs no store and runs in seconds.
#   EMPIRICAL  For each tracked input, change the file's bytes, ask
#              `tar_outdated()` what it now considers stale, and restore. This
#              is what `targets` ACTUALLY does. It requires an up-to-date store
#              and is refused without one, because an already-stale target
#              would be attributed to the perturbation.
#
# The two agreeing is the point. Static alone would not notice a `targets`
# behaviour that departs from the graph; empirical alone gives a number with no
# statement of what it should have been.
#
# HOW THE PERTURBATION IS MADE SAFE. Each file is COPIED to `<file>.audit_bak`
# first, then one byte is appended, then the copy is restored and the md5 is
# verified against the pre-perturbation md5. Restoration is registered with
# `on.exit()` so an error or an interrupt still restores. If a previous run was
# killed hard, `.audit_bak` files are left behind and this script refuses to
# start until it has restored from them -- which it does automatically. The
# five parquets total about 31 MB, so the copy is a fraction of a second and
# the file is never left in a modified state for longer than one
# `tar_outdated()` call.
#
# APPENDING A BYTE MAKES A PARQUET TEMPORARILY UNREADABLE, because the format's
# footer is at the end of the file. That is safe here for one reason and it
# should be stated rather than assumed: `tar_outdated()` HASHES files, it does
# not parse them, and it runs no target command. Nothing reads the perturbed
# bytes.
#
# NO DATA IS READ, AND NOTHING LEAVES THE MACHINE (hard rule 1). `file.copy()`
# writes the backup into the SAME directory on the same disk; the perturbation
# is one byte appended in place; `tools::md5sum` hashes locally. No row, column
# or value is ever loaded, nothing is transmitted, and the only thing this
# script prints is file basenames and counts of outdated targets.
#
# THE EMPIRICAL ARM IS THEREFORE OPT-IN, behind `--touch-data`. It is safe and
# it is verified, but writing to a file under a credentialed DUA is a thing to
# be asked for rather than defaulted to.
#
#   Rscript tests/audit_a_invalidation.R              static only, ~15 s
#   Rscript tests/audit_a_invalidation.R --touch-data static + empirical, ~4 min
#   Rscript tests/audit_a_invalidation.R --bless      write the expectation file
#
# The expectation lives in `tests/audit_a_expected.csv`, one row per
# (input, target). It is written by `--bless` and reviewed by a human once;
# after that, any drift in either direction is a failure of this script.
# ----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  root <- "."
  source(file.path(root, "tests", "audit_common.R"))
  library(targets)
}))

args    <- commandArgs(trailingOnly = TRUE)
# THE EMPIRICAL ARM IS OPT-IN. It modifies files under `data/` in place --
# safely, reversibly, and verified by md5 -- but they are under a PhysioNet
# credentialed DUA and a script that writes to them must be asked for
# explicitly rather than run by default. The static arm answers the question
# this audit exists for (which targets a change would invalidate) and touches
# nothing.
do_emp  <- ("--touch-data" %in% args) && !("--static" %in% args)
do_bless <- "--bless" %in% args
EXPECT  <- file.path(root, "tests", "audit_a_expected.csv")

# --- the tracked inputs ------------------------------------------------------
# Derived from the graph rather than listed, so an input added to `_targets.R`
# enters this matrix without anybody remembering to add it here.

man <- targets::tar_manifest(fields = c("name", "command", "format"),
                             callr_function = NULL)
net <- targets::tar_network(targets_only = TRUE, callr_function = NULL)
edge <- as.data.frame(net$edges, stringsAsFactors = FALSE)
file_targets <- man$name[man$format == "file"]

cfg <- yaml::read_yaml(file.path(root, "config", "config.yml"))
# Which files each file target covers. `mimic_files` expands to the five
# parquets, and each is probed separately: a parquet missing from
# `paths$mimiciv` would not be tracked at all, and only a per-file probe says so.
input_files <- list(
  config_file  = file.path(root, "config", "config.yml"),
  pairing_file = file.path(root, "config", "pairing.csv"),
  domains_file = file.path(root, "config", "domains.csv"),
  mimic_files  = file.path(root, unlist(cfg$paths$mimiciv, use.names = FALSE)),
  # THE CODE, tracked since 2026-09-03 (audit finding F3). Probed on ONE file
  # rather than all twenty: unlike the parquets, which are five separate
  # extractions any of which could be re-run alone, these are one directory
  # listing and one file is enough to show the target invalidates. The static
  # arm covers the reachable set for all of them.
  source_files = file.path(root, "R", "00_utils.R")
)
missing_tracked <- setdiff(file_targets, names(input_files))
if (length(missing_tracked)) {
  stop("audit A: file target(s) with no declared file list: ",
       paste(missing_tracked, collapse = ", "),
       ". Add them to `input_files` above.", call. = FALSE)
}

# --- static arm --------------------------------------------------------------

descendants <- function(seed) {
  seen <- character(0); st <- seed
  while (length(st)) {
    n <- st[1]; st <- st[-1]
    if (n %in% seen) next
    seen <- c(seen, n)
    st <- c(st, setdiff(edge$to[edge$from == n], seen))
  }
  sort(setdiff(seen, seed))
}
static <- lapply(stats::setNames(nm = names(input_files)), descendants)

# --- declared invariants -----------------------------------------------------
# These hold regardless of what the blessed expectation file says, and they are
# the assertions that carry the meaning. The expectation file catches drift;
# these catch the specific failures the project has already had or would care
# most about.

FIT_TARGETS <- c("tabs", "priors", "layer1_oof_branch", "layer1_final_raw",
                 "final_models", "l_oof", "bundle")

inv <- list(
  list(id = "A1",
       claim = "every tracked input has at least one dependent",
       ok = all(lengths(static) > 0L),
       detail = paste(names(static), lengths(static), sep = "=", collapse = ", ")),
  list(id = "A2",
       claim = "a re-extraction (mimic_files) rebuilds every fit and the bundle",
       ok = all(FIT_TARGETS %in% static$mimic_files),
       detail = paste("missing:", paste(setdiff(FIT_TARGETS, static$mimic_files),
                                        collapse = ", "))),
  list(id = "A3",
       claim = "a config edit (config_file) rebuilds every fit and the bundle",
       ok = all(FIT_TARGETS %in% static$config_file),
       detail = paste("missing:", paste(setdiff(FIT_TARGETS, static$config_file),
                                        collapse = ", "))),
  list(id = "A4",
       claim = "a pairing.csv edit rebuilds every fit and the bundle",
       ok = all(FIT_TARGETS %in% static$pairing_file),
       detail = paste("missing:", paste(setdiff(FIT_TARGETS, static$pairing_file),
                                        collapse = ", "))),
  list(id = "A5",
       claim = "domains.csv reaches the bundle",
       ok = "bundle" %in% static$domains_file,
       detail = paste(static$domains_file, collapse = ", ")),
  # Added with the fix for audit finding F3. The bundle carries `source_hashes`
  # as provenance; if an edit to `R/` did not reach the bundle, that provenance
  # could describe code other than the code that built it.
  list(id = "A5b",
       claim = "an edit to R/ reaches the bundle's provenance",
       ok = all(c("code_hashes", "bundle") %in% static$source_files),
       detail = paste(static$source_files, collapse = ", ")),
  list(id = "A6",
       claim = "every target is downstream of at least one tracked input",
       ok = length(setdiff(man$name, c(names(input_files),
                                       unlist(static, use.names = FALSE)))) == 0L,
       detail = paste("orphans:", paste(setdiff(man$name,
                        c(names(input_files), unlist(static, use.names = FALSE))),
                        collapse = ", ")))
)
inv_df <- do.call(rbind, lapply(inv, function(z) data.frame(
  id = z$id, claim = z$claim, verdict = if (isTRUE(z$ok)) "PASS" else "FAIL",
  detail = if (isTRUE(z$ok)) "" else z$detail, stringsAsFactors = FALSE)))

# --- empirical arm -----------------------------------------------------------

.md5 <- function(p) unname(tools::md5sum(p))

#' Restore any `.audit_bak` left by a hard-killed previous run, before anything
#' else happens. A leftover backup means a file is in a perturbed state, and
#' starting a new probe on top of that would compare against a corrupted input.
recover_backups <- function(paths) {
  n <- 0L
  for (p in paths) {
    bak <- paste0(p, ".audit_bak")
    if (file.exists(bak)) {
      file.copy(bak, p, overwrite = TRUE)
      unlink(bak)
      n <- n + 1L
      message("audit A: RECOVERED ", p, " from a leftover backup")
    }
  }
  n
}

#' Perturb one file, ask `tar_outdated()` what went stale, restore, verify.
probe <- function(path) {
  bak <- paste0(path, ".audit_bak")
  before <- .md5(path)
  size0  <- file.info(path)$size
  if (!file.copy(path, bak, overwrite = TRUE)) {
    stop("audit A: could not back up ", path, call. = FALSE)
  }
  on.exit({
    if (file.exists(bak)) {
      file.copy(bak, path, overwrite = TRUE)
      unlink(bak)
    }
  }, add = TRUE)

  con <- file(path, open = "ab")
  writeBin(as.raw(10L), con)          # one newline byte; hash changes, nothing reads it
  close(con)

  out <- targets::tar_outdated(callr_function = NULL)

  file.copy(bak, path, overwrite = TRUE)
  unlink(bak)
  after <- .md5(path)
  if (!identical(before, after) || !identical(size0, file.info(path)$size)) {
    stop("audit A: FAILED TO RESTORE ", path,
         " -- restore it from ", bak, " by hand.", call. = FALSE)
  }
  sort(out)
}

empirical <- NULL
emp_note  <- ""
if (do_emp) {
  all_paths <- unlist(input_files, use.names = FALSE)
  miss <- all_paths[!file.exists(all_paths)]
  if (length(miss)) {
    emp_note <- paste0("SKIPPED: input file(s) absent: ", paste(miss, collapse = ", "))
  } else {
    recover_backups(all_paths)
    base_out <- targets::tar_outdated(callr_function = NULL)
    if (length(base_out)) {
      emp_note <- paste0(
        "SKIPPED: the store is not up to date (", length(base_out),
        " target(s) already outdated: ", paste(utils::head(base_out, 8), collapse = ", "),
        "). Run `targets::tar_make()` first -- an already-stale target would be ",
        "attributed to the perturbation and the matrix would be meaningless.")
    } else {
      empirical <- list()
      for (tgt in names(input_files)) {
        for (p in input_files[[tgt]]) {
          got <- probe(p)
          empirical[[length(empirical) + 1L]] <- list(
            input = tgt, file = basename(p), n = length(got), targets = got)
          cat(sprintf("  probed %-46s -> %2d target(s) outdated\n", basename(p), length(got)))
        }
      }
    }
  }
}

# --- static vs empirical -----------------------------------------------------

cmp <- NULL
if (!is.null(empirical)) {
  cmp <- do.call(rbind, lapply(empirical, function(z) {
    s <- static[[z$input]]
    # The file target itself rebuilds too, and is not in its own descendant set.
    e <- setdiff(z$targets, z$input)
    data.frame(input = z$input, file = z$file,
               n_static = length(s), n_empirical = length(e),
               agree = setequal(s, e),
               only_empirical = paste(setdiff(e, s), collapse = ", "),
               only_static = paste(setdiff(s, e), collapse = ", "),
               stringsAsFactors = FALSE)
  }))
}

# --- the blessed expectation -------------------------------------------------

exp_tbl <- do.call(rbind, lapply(names(static), function(k)
  if (length(static[[k]])) data.frame(input = k, target = static[[k]],
                                      stringsAsFactors = FALSE) else NULL))

drift <- NULL
if (do_bless || !file.exists(EXPECT)) {
  utils::write.csv(exp_tbl, EXPECT, row.names = FALSE)
  bless_note <- sprintf("expectation %s to %s (%d rows). REVIEW IT.",
                        if (do_bless) "re-written" else "written for the first time",
                        EXPECT, nrow(exp_tbl))
} else {
  old <- utils::read.csv(EXPECT, colClasses = "character")
  key_new <- paste(exp_tbl$input, exp_tbl$target)
  key_old <- paste(old$input, old$target)
  add <- exp_tbl[!key_new %in% key_old, , drop = FALSE]
  rem <- old[!key_old %in% key_new, , drop = FALSE]
  drift <- rbind(
    if (nrow(add)) data.frame(change = "now rebuilds", add, stringsAsFactors = FALSE),
    if (nrow(rem)) data.frame(change = "no longer rebuilds", rem, stringsAsFactors = FALSE))
  bless_note <- if (is.null(drift) || !nrow(drift))
    sprintf("matches the blessed expectation (%d rows).", nrow(exp_tbl))
  else sprintf("DRIFT against %s: %d change(s).", EXPECT, nrow(drift))
}

# --- report ------------------------------------------------------------------

md <- c("# Audit A — the input-invalidation matrix", "",
        "Generated by `tests/audit_a_invalidation.R`. For each tracked input,",
        "the set of targets that would rebuild if it changed — computed from the",
        "graph (static) and, where the store is up to date, measured by perturbing",
        "the file and asking `tar_outdated()` (empirical).", "",
        "No data is read: files under `data/` are copied and hashed as bytes",
        "(hard rule 1).", "")

md <- c(md, "## Declared invariants", "", md_table(inv_df), "")

md <- c(md, "## Static — targets reachable from each input", "",
        md_table(data.frame(input = names(static),
                            n_targets = lengths(static),
                            targets = vapply(static, paste, character(1), collapse = " "),
                            stringsAsFactors = FALSE)), "")

md <- c(md, "## Empirical — `tar_outdated()` after perturbing each file", "")
if (!is.null(cmp)) {
  md <- c(md, md_table(cmp[, c("input", "file", "n_static", "n_empirical",
                               "agree", "only_empirical", "only_static")]), "")
} else {
  md <- c(md, paste0("_", if (nzchar(emp_note)) emp_note else "not run (--static)", "_"), "")
}

md <- c(md, "## Expectation file", "", paste0("`", EXPECT, "` — ", bless_note), "")
if (!is.null(drift) && nrow(drift)) md <- c(md, md_table(drift), "")

audit_write(md, "audit_a_invalidation.md", root)

cat("\n=== AUDIT A: input-invalidation matrix ===\n\n")
print(inv_df[, c("id", "verdict", "claim")], row.names = FALSE)
cat("\nstatic reachable sets:\n")
print(data.frame(input = names(static), n_targets = lengths(static)), row.names = FALSE)
cat("\nempirical:\n")
if (!is.null(cmp)) print(cmp[, c("input", "file", "n_static", "n_empirical", "agree")],
                         row.names = FALSE) else cat("  ", emp_note, "\n")
cat("\nexpectation: ", bless_note, "\n", sep = "")
if (!is.null(drift) && nrow(drift)) print(drift, row.names = FALSE)
if (any(inv_df$verdict == "FAIL")) {
  cat("\nFAILING INVARIANTS:\n")
  print(inv_df[inv_df$verdict == "FAIL", ], row.names = FALSE)
}
