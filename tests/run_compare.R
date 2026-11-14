# tests/run_compare.R ----------------------------------------------------------
# Compare two run directories table by table. FITS NOTHING, READS AGGREGATES.
#
# Added 2026-10-25 for the code-review reruns. A rerun meant to change a named
# set of tables (the severity curves, after the frozen cut points were wired)
# must leave every other table identical, and "identical" is checked here
# rather than asserted: every CSV under tables/ and diagnostics/ present in
# either run is compared cell by cell after reading, and the verdict per file
# is IDENTICAL, DIFFERENT (with the columns that differ), ONLY_OLD or ONLY_NEW.
#
#   Rscript tests/run_compare.R out/runs/test_OLD out/runs/test_NEW [expected_regex]
#
# `expected_regex` names the files that are allowed to differ; the script exits
# non-zero when any other file differs or disappears. Aggregate CSVs only: the
# .rds objects of a run can hold row-level values and are never opened here.
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: run_compare.R OLD NEW [expected_regex]", call. = FALSE)
OLD <- args[1]; NEW <- args[2]; EXPECT <- if (length(args) >= 3L) args[3] else "^$"

list_csv <- function(d) {
  f <- unlist(lapply(c("tables", "diagnostics"), function(s)
    if (dir.exists(file.path(d, s))) file.path(s, list.files(file.path(d, s), pattern = "\\.csv$")) else character(0)))
  sort(f)
}
fo <- list_csv(OLD); fn <- list_csv(NEW)
rows <- list()
for (f in union(fo, fn)) {
  st <- if (!f %in% fn) "ONLY_OLD" else if (!f %in% fo) "ONLY_NEW" else {
    a <- read.csv(file.path(OLD, f), stringsAsFactors = FALSE, check.names = FALSE)
    b <- read.csv(file.path(NEW, f), stringsAsFactors = FALSE, check.names = FALSE)
    if (identical(a, b)) "IDENTICAL" else {
      cols <- union(names(a), names(b))
      dif <- cols[vapply(cols, function(cn) !identical(a[[cn]], b[[cn]]), logical(1))]
      paste0("DIFFERENT [", paste(utils::head(dif, 8), collapse = ", "),
             if (length(dif) > 8) ", ..." else "", "]",
             if (nrow(a) != nrow(b)) sprintf(" rows %d -> %d", nrow(a), nrow(b)) else "")
    }
  }
  rows[[length(rows) + 1L]] <- data.frame(file = f, status = st, expected = grepl(EXPECT, f),
                                          stringsAsFactors = FALSE)
}
R <- do.call(rbind, rows)
cat(sprintf("%s\n  vs %s\n\n", OLD, NEW))
cat(sprintf("  identical: %d   different: %d   only old: %d   only new: %d\n\n",
            sum(R$status == "IDENTICAL"), sum(startsWith(R$status, "DIFFERENT")),
            sum(R$status == "ONLY_OLD"), sum(R$status == "ONLY_NEW")))
z <- R[R$status != "IDENTICAL", , drop = FALSE]
if (nrow(z)) for (i in seq_len(nrow(z))) cat(sprintf("  %-9s %-60s %s\n",
  if (z$expected[i]) "expected" else "UNEXPECTED", z$file[i], z$status[i]))
bad <- z[!z$expected & z$status != "ONLY_NEW", , drop = FALSE]
if (nrow(bad)) {
  cat(sprintf("\n%d unexpected difference(s).\n", nrow(bad))); quit(status = 1L)
}
cat("\nno unexpected difference.\n")
