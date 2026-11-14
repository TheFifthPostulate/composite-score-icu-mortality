# tests/audit_all.R ----------------------------------------------------------
# Run every mechanical audit of docs/v2_audit_plan_20260903.md section 2, in
# the plan's order, and report which produced findings.
#
# Each audit runs in its own R process, so one failing does not take the rest
# down with it -- the point of the suite is a findings LIST, and an audit that
# errors is itself a finding rather than a reason to stop.
#
#   Rscript tests/audit_all.R           B, C, D, E, then A static      ~1 min
#   Rscript tests/audit_all.R --full    the same, plus A's empirical arm, which
#                                       WRITES TO AND RESTORES FILES UNDER data/  ~5 min
#
# Reports land in out/audit_*.md. Invariance items 1, 3 and 5 are separate and
# cost minutes to hours: see tests/invariance_suite.R.
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args

audits <- list(
  list(id = "B", script = "tests/audit_b_config.R",      args = character(0),
       what = "config reachability"),
  list(id = "C", script = "tests/audit_c_targets.R",     args = character(0),
       what = "target dependency hygiene, hard rules 7 and 9"),
  list(id = "D", script = "tests/audit_d_columns.R",     args = character(0),
       what = "hard-coded column lists against CONTRACT"),
  list(id = "E", script = "tests/audit_e_fitsettings.R", args = character(0),
       what = "fit-setting provenance"),
  list(id = "A", script = "tests/audit_a_invalidation.R",
       args = if (full) "--touch-data" else character(0),
       what = "input-invalidation matrix")
)

rows <- list()
for (a in audits) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("AUDIT ", a$id, " -- ", a$what, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")
  st <- system2("Rscript", c(a$script, a$args), stdout = "", stderr = "")
  rows[[length(rows) + 1L]] <- data.frame(
    audit = a$id, what = a$what, status = if (identical(st, 0L)) "ran" else "ERRORED",
    stringsAsFactors = FALSE)
}

cat("\n", strrep("=", 78), "\n", sep = "")
print(do.call(rbind, rows), row.names = FALSE)
cat("\nreports: out/audit_a_invalidation.md, out/audit_b_config.md,\n")
cat("         out/audit_c_targets.md, out/audit_d_columns.md,\n")
cat("         out/audit_e_fitsettings.md\n")
cat("\nAudits F and invariance item 6 are section 6 step 5 and are not built.\n")
