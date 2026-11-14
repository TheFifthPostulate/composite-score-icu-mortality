# tests/audit_common.R -------------------------------------------------------
# Shared machinery for the mechanical audits of docs/v2_audit_plan_20260903.md
# section 2. Sourced by tests/audit_b_config.R, audit_c_targets.R,
# audit_d_columns.R and audit_e_fitsettings.R. Runs nothing on its own.
#
# WHY THIS FILE EXISTS AT ALL. Every audit in the plan is a question about a
# RELATIONSHIP between a declaration and the code that consumes it, so every
# audit needs the same two primitives: an inventory of the declarations, and a
# search over the code that does not lie about where a name occurs. The second
# is the harder one and is the reason this is not four independent greps.
#
# THE COMMENT PROBLEM, WHICH IS THE WHOLE POINT. A config key named in a
# comment is NOT read by anything. `signals_not_transportable` is documented at
# length in config.yml and discussed in R/ comments, and it was still dead. A
# naive grep would have called it reachable and the audit would have reported
# clean. So `scan_sources()` strips comments with a character walk that tracks
# string state, rather than with a regex, because `"#"` inside a string literal
# is not a comment and `# "` inside a comment does not open one.
#
# NOTHING HERE READS data/ OR _targets/objects/ (hard rule 1). These audits are
# static analysis over source and config text only. Audit A is the one that
# needs the store, and it lives in its own script.
# ----------------------------------------------------------------------------

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# --- source inventory -------------------------------------------------------

#' Every source file the audits consider, tagged by the LAYER it sits in.
#'
#' The layer is what separates a "read" verdict from a "vestigial" one: a key
#' read only by `tests/` is read only by a path that no runner traverses, which
#' is exactly the state `smooth_k`'s constructs-off entries are in.
audit_sources <- function(root = ".") {
  f <- c(
    list.files(file.path(root, "R"),     pattern = "[.][Rr]$", full.names = TRUE),
    list.files(file.path(root, "run"),   pattern = "[.][Rr]$", full.names = TRUE),
    list.files(file.path(root, "tests"), pattern = "[.][Rr]$", full.names = TRUE),
    file.path(root, "_targets.R")
  )
  f <- f[file.exists(f)]
  # The audits themselves are not part of the pipeline and must not make a key
  # look reachable just by naming it in a findings table.
  f <- f[!grepl("audit_[^/]*[.][Rr]$", gsub("\\", "/", f, fixed = TRUE))]
  p <- gsub("\\", "/", f, fixed = TRUE)
  layer <- ifelse(grepl("/tests/", p, fixed = TRUE), "tests",
           ifelse(grepl("/run/",   p, fixed = TRUE), "run",
           ifelse(grepl("_targets.R", p, fixed = TRUE), "targets", "lib")))
  data.frame(path = p, layer = layer, stringsAsFactors = FALSE)
}

#' Strip R comments from one line, respecting string literals.
#'
#' Returns the code part of the line. A `#` inside "..." or '...' is data; the
#' first `#` outside one starts a comment. Backslash escapes are honoured.
#' Backticked names are treated as a third quote type so that a `#` inside one
#' does not open a comment either.
.strip_comment_line <- function(line) {
  if (!grepl("#", line, fixed = TRUE)) return(line)
  ch <- strsplit(line, "", fixed = TRUE)[[1]]
  if (!length(ch)) return("")
  quote <- ""
  i <- 1L
  n <- length(ch)
  bs <- "\\"
  while (i <= n) {
    c1 <- ch[i]
    if (nzchar(quote)) {
      if (identical(c1, bs)) { i <- i + 2L; next }
      if (identical(c1, quote)) quote <- ""
    } else {
      if (c1 %in% c("\"", "'", "`")) quote <- c1
      else if (identical(c1, "#")) {
        if (i == 1L) return("")
        return(paste(ch[seq_len(i - 1L)], collapse = ""))
      }
    }
    i <- i + 1L
  }
  line
}

#' Read every source file, comments removed, as a long data frame of code lines.
#'
#' One row per line so that a hit can be reported as file:line, which is what
#' makes an adjudication cheap.
scan_sources <- function(root = ".") {
  src <- audit_sources(root)
  out <- lapply(seq_len(nrow(src)), function(i) {
    raw <- readLines(src$path[i], warn = FALSE)
    code <- vapply(raw, .strip_comment_line, character(1), USE.NAMES = FALSE)
    data.frame(path = src$path[i], layer = src$layer[i],
               line = seq_along(raw), code = code, raw = raw,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' A PCRE fragment matching `x` literally. `\Q...\E` rather than hand-escaping,
#' because hand-escaping a bracket class is exactly the sort of thing that goes
#' wrong silently and makes an audit report clean when it is not.
.lit <- function(x) paste0("\\Q", x, "\\E")

#' Every line on which `name` is USED as a name, not merely as a substring.
#'
#' Matches the four ways this codebase reaches a config entry: `$name`,
#' `[["name"]]`, a bare quoted `"name"` (which covers `setdiff(..., "name")`
#' and `switch` arms), and a bare symbol `name` (which covers a local variable
#' bound from the key). The bare-symbol form is deliberately included and is
#' the reason a hit is ADJUDICATED rather than trusted: it over-matches, and an
#' over-matching reachability audit fails safe while an under-matching one
#' declares a dead key live.
grep_name <- function(code, name) {
  q <- .lit(name)
  pat <- paste0("([$]\\s*", q, "\\b)",
                "|(\\[\\[\\s*[\"']", q, "[\"']\\s*\\]\\])",
                "|([\"']", q, "[\"'])",
                "|(\\b", q, "\\b)")
  grepl(pat, code$code, perl = TRUE)
}

#' As `grep_name()` but WITHOUT the bare-symbol arm: only an explicit `$name`,
#' `[["name"]]` or quoted `"name"` counts. Used to separate a genuine config
#' read from a coincidental local variable of the same name.
grep_name_strict <- function(code, name) {
  q <- .lit(name)
  pat <- paste0("([$]\\s*", q, "\\b)",
                "|(\\[\\[\\s*[\"']", q, "[\"']\\s*\\]\\])",
                "|([\"']", q, "[\"'])")
  grepl(pat, code$code, perl = TRUE)
}

#' Which SYNTACTIC ARM matched, per line, for one name.
#'
#' The arm matters because it separates a real config read from a coincidence.
#' A nested config key is normally reached as `$key`; a hit on the quoted arm
#' alone means the name appeared as a bare string literal somewhere, which for
#' a short generic key like `out` or `seed` is far more likely to be an
#' unrelated `"out"` than a config read. `paths.out` is the worked example: it
#' scored "reachable" on the quoted arm alone, from `out_dir <- "out"` in a
#' test, and it is in fact read by nothing.
grep_arms <- function(code, name) {
  q <- .lit(name)
  list(
    dollar  = grepl(paste0("[$]\\s*", q, "\\b"), code$code, perl = TRUE),
    bracket = grepl(paste0("\\[\\[\\s*[\"']", q, "[\"']\\s*\\]\\]"), code$code, perl = TRUE),
    quoted  = grepl(paste0("[\"']", q, "[\"']"), code$code, perl = TRUE),
    symbol  = grepl(paste0("\\b", q, "\\b"), code$code, perl = TRUE)
  )
}

#' Every `cfg_req()` / `cfg_flag()` call, as the dotted config path it reads.
#'
#' WHY THIS EXISTS. Those accessors were introduced on 2026-09-03 to remove the
#' `cfg$key %||% <literal>` double declaration (audit findings F9, F10, F11).
#' They reach a key by QUOTED ARGUMENT rather than by `$` or `[[`, so the
#' accessor matcher cannot see them, and audit B immediately reported three
#' live keys as dead or vestigial. That is the audit working -- a matcher that
#' does not understand an idiom must report a miss rather than a pass -- and
#' the fix is to teach it the idiom rather than to loosen the matcher.
#'
#' AST-based, not line-based, so a call split across lines is still found.
#' Returns one row per (path, file, line, layer).
cfg_req_paths <- function(root = ".") {
  src <- audit_sources(root)
  out <- list()
  for (i in seq_len(nrow(src))) {
    ex <- tryCatch(parse(src$path[i], keep.source = TRUE), error = function(e) NULL)
    if (is.null(ex)) next
    refs <- utils::getSrcref(ex)
    for (j in seq_along(ex)) {
      line <- if (!is.null(refs) && length(refs) >= j) as.integer(refs[[j]])[1] else NA_integer_
      walk <- function(x) {
        if (is.call(x)) {
          fn <- tryCatch(as.character(x[[1]])[1], error = function(e) "")
          if (fn %in% c("cfg_req", "cfg_flag")) {
            args <- as.list(x)[-1]
            nms  <- names(args)
            keep <- vapply(seq_along(args), function(k) {
              if (!is.null(nms) && nzchar(nms[k] %||% "")) return(FALSE)
              a <- tryCatch(args[[k]], error = function(e) NULL)
              is.character(a) && length(a) == 1L
            }, logical(1))
            path <- unlist(args[keep], use.names = FALSE)
            if (length(path)) {
              out[[length(out) + 1L]] <<- data.frame(
                path = paste(path, collapse = "."), key = path[length(path)],
                file = src$path[i], line = line, layer = src$layer[i],
                stringsAsFactors = FALSE)
            }
          }
        }
        if (is.call(x) || is.pairlist(x) || is.expression(x) || is.list(x)) {
          l <- as.list(x)
          for (k in seq_along(l)) {
            if (identical(unname(l[k]), list(quote(expr = )))) next
            walk(l[[k]])
          }
        }
        invisible()
      }
      walk(ex[[j]])
    }
  }
  if (!length(out)) {
    return(data.frame(path = character(0), key = character(0), file = character(0),
                      line = integer(0), layer = character(0)))
  }
  unique(do.call(rbind, out))
}

#' Every `cfg$<key>` reference in the given layers, as a character vector.
#'
#' This is the inventory `BUNDLE_DESIGN_KEYS` claims to be derived from, and
#' having it as a function is what turns that claim into a checkable assertion
#' rather than a comment.
cfg_refs <- function(code, layers = "lib") {
  z <- code[code$layer %in% layers, ]
  m <- gregexpr("(?<=\\bcfg[$])[A-Za-z_][A-Za-z0-9_.]*", z$code, perl = TRUE)
  sort(unique(unlist(regmatches(z$code, m))))
}

#' Verdict for one declared name, given the layers it is reached from.
#'
#' read       reached from R/, run/ or _targets.R
#' vestigial  reached only from tests/
#' dead       reached from nowhere
reach_verdict <- function(layers) {
  if (!length(layers)) return("dead")
  if (any(layers %in% c("lib", "run", "targets"))) "read" else "vestigial"
}

#' Print a data frame as a GitHub markdown table.
md_table <- function(df, digits = 4) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) df[[j]] <- format(round(df[[j]], digits), trim = TRUE)
    df[[j]] <- as.character(df[[j]])
    df[[j]][is.na(df[[j]])] <- ""
    df[[j]] <- gsub("|", "\\|", df[[j]], fixed = TRUE)
  }
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- vapply(seq_len(nrow(df)), function(i)
    paste0("| ", paste(unlist(df[i, ], use.names = FALSE), collapse = " | "), " |"),
    character(1))
  c(hdr, sep, rows)
}

#' Write a report to out/ and echo its headline counts. Audits are scripts, not
#' library code, so building a path here does not violate hard rule 9.
audit_write <- function(lines, file, root = ".") {
  dir.create(file.path(root, "out"), showWarnings = FALSE, recursive = TRUE)
  p <- file.path(root, "out", file)
  writeLines(lines, p)
  cat(sprintf("\nwrote %s  (%d lines)\n", p, length(lines)))
  invisible(p)
}
