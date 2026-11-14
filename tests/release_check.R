# tests/release_check.R ---------------------------------------------------------
# What the repository would publish, and whether any of it must not go out (2026-11-13).
#
# Lists the files git would publish: those already staged or committed, plus those `git add .` would stage
# (not ignored by .gitignore). Fails, with exit status 1, if any of them
#   - lies under a folder that holds credentialed data or objects built from it (data/, out/, _targets/, ...),
#   - has a data-file extension (.rds, .qs2, .parquet, .zip, a rendered notebook, ...),
#   - is a CSV not on the reviewed list below (a new CSV needs a small-count review before it is added),
#   - is larger than 5 MB, or
#   - contains an absolute path into a user folder.
# It prints file paths and counts only, never file contents. Run it from the project root before every
# commit. Before `git init` it uses a temporary git directory and leaves the project untouched.
#
#   Rscript tests/release_check.R            # summary and verdict
#   Rscript tests/release_check.R --list     # also every path it would publish
# ------------------------------------------------------------------------------
if (!file.exists(".gitignore") || !dir.exists("R")) stop("run from the project root, next to .gitignore", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)

FORBIDDEN_DIRS <- c("data/", "out/", "_targets/", "notebooks/", "demo/", ".venv-ebm/", ".claude/",
                    "thinking_archive/", "paper/_local_credentialed_do_not_publish/")
FORBIDDEN_TEXT <- c("/_backup", "credentialed")                  # anywhere in a path
DATA_EXT <- "\\.(rds|rda|rdata|qs|qs2|parquet|feather|arrow|fst|zip|gz|7z|tar|pkl|pickle|db|sqlite|duckdb|xls|xlsx|sav|dta|nb\\.html)$"
SAFE_CSV <- c("config/pairing.csv", "config/domains.csv",          # design declarations
              "tests/audit_a_expected.csv",                        # config and target names
              "tests/invariance_determinism_snapshot.csv")         # 16-character hashes
MAX_MB   <- 5
TEXT_EXT <- "\\.(R|py|sql|yml|yaml|md|txt|csv|tex)$"
USERPATH <- "[A-Za-z]:[/\\\\]+Users[/\\\\]|/Users/[A-Za-z]|/home/[A-Za-z]"

git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status"))) stop("git ", paste(c(...), collapse = " "), " failed:\n", paste(out, collapse = "\n"), call. = FALSE)
  out
}
if (dir.exists(".git")) {
  files <- git("-c", "core.quotepath=off", "ls-files", "--cached", "--others", "--exclude-standard")
  how <- "staged or committed, plus what `git add .` would stage"
} else {
  tmp <- tempfile("release_check_git")
  git("init", "-q", "--bare", tmp)
  files <- git("-c", "core.quotepath=off", "--git-dir", tmp, "--work-tree", ".", "ls-files", "--others", "--exclude-standard")
  unlink(tmp, recursive = TRUE)
  how <- "what `git add .` would stage after `git init` (no repository yet)"
}
files <- sort(unique(files[nzchar(files)]))
size  <- file.size(files)

problems <- character(0)
flag <- function(f, why) problems <<- c(problems, sprintf("%-60s %s", f, why))
for (i in seq_along(files)) {
  f <- files[i]; lf <- tolower(f)
  if (any(startsWith(f, FORBIDDEN_DIRS))) flag(f, "inside a folder that must stay private")
  else if (any(vapply(FORBIDDEN_TEXT, grepl, TRUE, x = lf, fixed = TRUE))) flag(f, "a backup or credentialed path")
  if (grepl(DATA_EXT, lf)) flag(f, "a data-file extension")
  if (grepl("\\.csv$", lf) && !f %in% SAFE_CSV) flag(f, "a CSV not on the reviewed list (SAFE_CSV)")
  if (!is.na(size[i]) && size[i] > MAX_MB * 1e6) flag(f, sprintf("%.1f MB, over %d MB", size[i] / 1e6, MAX_MB))
  if (!is.na(size[i]) && grepl(TEXT_EXT, f) && any(grepl(USERPATH, readLines(f, warn = FALSE))))
    flag(f, "contains an absolute path into a user folder")
}

top <- ifelse(grepl("/", files), sub("/.*", "/", files), "(top level)")
cat(sprintf("release check: %d files, %.2f MB -- %s\n", length(files), sum(size, na.rm = TRUE) / 1e6, how))
tab <- table(top); cat(sprintf("  %-14s %4d\n", names(tab), as.integer(tab)), sep = "")
if ("--list" %in% args) cat(paste0("    ", files), sep = "\n")
if (length(problems)) {
  cat("\nNOT READY. These would be published and must not be:\n", paste0("  ", problems, "\n"), sep = "")
  cat("Fix .gitignore (or `git rm --cached <path>` for a staged file), then run this again.\n")
  quit(save = "no", status = 1L)
}
cat("\nready: nothing on the forbidden list would be published\n")
