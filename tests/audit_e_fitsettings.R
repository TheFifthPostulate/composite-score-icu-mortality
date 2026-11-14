# tests/audit_e_fitsettings.R ------------------------------------------------
# AUDIT E of docs/v2_audit_plan_20260903.md section 2: fit-setting provenance.
#
# THE QUESTION. Every model-fitting knob should come from config and should
# appear in exactly one place. The failure class is a hard-coded value sitting
# beside a config key that declares the same thing: the two agree today, they
# are edited apart tomorrow, and which answer you get depends on which entry
# point you came in through.
#
# THE CODEBASE'S IDIOM IS `cfg$group$key %||% <literal>`, and that idiom is most
# of the audit. It is a good idiom -- the code keeps running against an older
# config -- but it means EVERY knob has two declarations: the config entry and
# the fallback literal. They must agree, and nothing currently checks that they
# do. When they disagree, the divergence is invisible until someone deletes or
# misspells the config key, at which point the pipeline silently fits a
# different model rather than failing.
#
#   E1   Every `<cfg path> %||% <literal>` in `R/`, resolved through local
#        aliases (`bs <- cfg$bam`, `x <- cfg$xgboost`), with the literal
#        compared against what `config/config.yml` actually declares.
#          agree     the fallback restates the declared value.
#          DISAGREE  they differ: a knob with two answers, and which one is in
#                    force depends on whether the config key is present.
#          no key    the fallback is the ONLY declaration, so the knob is in
#                    force and is not reviewable in config.
#
#   E2   Every argument of every `bam()` / `gam()` / `scam()` / `xgb.train()` /
#        `xgboost()` call, classified config-derived or literal. A literal here
#        is a setting no config edit can reach.
#
#   E3a  Named knob arguments taking a literal, in ANY call, with the CALLEE
#        reported -- which is what separates a model setting from an unrelated
#        `stats::cor(method = "spearman")`.
#
#   E3b  Smooth specifications built as STRINGS. `R/05_formula.R` writes
#        `sprintf("s(%s, bs = \"%s\", k = %d)", ...)`, so `bs` and `k` are never
#        arguments and no AST pass reaches them. This project's basis dimension
#        lives there, so a text pass is not a fallback here: it is the only
#        instrument that reaches the setting at all.
#
#   E3c  Function FORMAL DEFAULTS whose name matches a config key, compared
#        against the declared value. Same two-declarations problem as E1 wearing
#        a different syntax: `score_metrics(n_boot = 200L)` and
#        `metrics.n_boot: 200` are two places one number is written down.
#
# READS NO DATA (hard rule 1): config text and source text only.
#
#   Rscript tests/audit_e_fitsettings.R  ~5 s, writes out/audit_e_fitsettings.md
# ----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  root <- "."
  source(file.path(root, "tests", "audit_common.R"))
}))

cfg <- yaml::read_yaml(file.path(root, "config", "config.yml"))

# --- AST helpers -------------------------------------------------------------

walk_ast <- function(x, f) {
  f(x)
  if (is.call(x) || is.pairlist(x) || is.expression(x) || is.list(x)) {
    l <- as.list(x)
    for (i in seq_along(l)) {
      if (identical(unname(l[i]), list(quote(expr = )))) next
      walk_ast(l[[i]], f)
    }
  }
  invisible()
}

.fname <- function(x) tryCatch(as.character(x[[1]])[1], error = function(e) "")

#' The callee name, seeing through `pkg::fn`.
callee_of <- function(n) {
  fn <- .fname(n)
  if (is.call(n[[1]]) && .fname(n[[1]]) %in% c("::", ":::")) {
    fn <- tryCatch(as.character(n[[1]][[3]]), error = function(e) "")
  }
  fn
}

#' A `$` chain as a character path: `cfg$bam$k` -> c("cfg","bam","k").
#' NULL for anything that is not a pure chain of symbols.
dollar_path <- function(x) {
  if (is.symbol(x)) return(as.character(x))
  if (!is.call(x) || !identical(.fname(x), "$") || length(x) < 3L) return(NULL)
  lhs <- dollar_path(x[[2]])
  if (is.null(lhs)) return(NULL)
  rhs <- x[[3]]
  if (!is.symbol(rhs) && !is.character(rhs)) return(NULL)
  c(lhs, as.character(rhs))
}

is_literal <- function(x) {
  (is.numeric(x) || is.character(x) || is.logical(x)) && length(x) == 1L
}

fmt <- function(v) {
  if (is.null(v)) return("<null>")
  v <- unlist(v, use.names = FALSE)
  if (!length(v)) return("<empty>")
  paste(format(v, trim = TRUE), collapse = ",")
}

cfg_at <- function(path) {
  v <- cfg
  for (p in path) {
    if (!is.list(v) || is.null(v[[p]])) return(NULL)
    v <- v[[p]]
  }
  v
}

files <- sort(list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE))

#' Iterate top-level expressions of every file under R/, with a start line and
#' the name of the thing being defined.
for_each_top <- function(f_each) {
  for (f in files) {
    ex <- parse(f, keep.source = TRUE)
    refs <- utils::getSrcref(ex)
    for (j in seq_along(ex)) {
      top <- ex[[j]]
      line <- if (!is.null(refs) && length(refs) >= j) as.integer(refs[[j]])[1] else NA_integer_
      owner <- if (is.call(top) && .fname(top) %in% c("<-", "=", "<<-") &&
                   length(top) >= 3L && (is.symbol(top[[2]]) || is.character(top[[2]])))
        as.character(top[[2]]) else "<top level>"
      f_each(basename(f), line, owner, top)
    }
  }
}

# --- E1: `%||%` fallbacks against config ------------------------------------

e1 <- list()
for_each_top(function(file, line, owner, top) {
  # Local aliases: `bs <- cfg$bam`, `x <- cfg$xgboost %||% list()`.
  alias <- list()
  walk_ast(top, function(n) {
    if (!is.call(n) || !.fname(n) %in% c("<-", "=") || length(n) < 3L) return(invisible())
    if (!is.symbol(n[[2]])) return(invisible())
    rhs <- n[[3]]
    if (is.call(rhs) && identical(.fname(rhs), "%||%") && length(rhs) >= 3L) rhs <- rhs[[2]]
    p <- dollar_path(rhs)
    if (!is.null(p) && length(p) > 1L && identical(p[1], "cfg")) {
      alias[[as.character(n[[2]])]] <<- p[-1]
    }
    invisible()
  })

  walk_ast(top, function(n) {
    if (!is.call(n) || !identical(.fname(n), "%||%") || length(n) < 3L) return(invisible())
    lhs <- n[[2]]; rhs <- n[[3]]
    if (!is_literal(rhs)) return(invisible())
    p <- dollar_path(lhs)
    if (is.null(p) || length(p) < 2L) return(invisible())
    key <- if (identical(p[1], "cfg")) p[-1]
           else if (!is.null(alias[[p[1]]])) c(alias[[p[1]]], p[-1])
           else return(invisible())
    declared <- cfg_at(key)
    verdict <- if (is.null(declared)) "no key"
               else if (identical(fmt(declared), fmt(rhs))) "agree" else "DISAGREE"
    e1[[length(e1) + 1L]] <<- data.frame(
      file = file, line = line, owner = owner,
      key = paste(key, collapse = "."),
      fallback = fmt(rhs), declared = fmt(declared),
      verdict = verdict, stringsAsFactors = FALSE)
    invisible()
  })
})
e1 <- if (length(e1)) unique(do.call(rbind, e1)) else
  data.frame(file = character(0), line = integer(0), owner = character(0),
             key = character(0), fallback = character(0),
             declared = character(0), verdict = character(0))

# --- E2: fitter arguments ----------------------------------------------------

FITTERS <- c("bam", "gam", "scam", "xgb.train", "xgboost", "xgb.cv")

e2 <- list()
for_each_top(function(file, line, owner, top) {
  walk_ast(top, function(n) {
    if (!is.call(n)) return(invisible())
    fn <- callee_of(n)
    if (!fn %in% FITTERS) return(invisible())
    args <- as.list(n)[-1]
    anm  <- names(args) %||% rep("", length(args))
    for (i in seq_along(args)) {
      if (identical(unname(args[i]), list(quote(expr = )))) next
      a <- args[[i]]
      vars <- all.vars(a)
      src <- if (!length(vars)) "LITERAL"
             else if (length(intersect(vars, c("cfg", "bs", "x", "p", "ct")))) "config-derived"
             else "argument/derived"
      e2[[length(e2) + 1L]] <<- data.frame(
        file = file, line = line, call = fn,
        arg = if (nzchar(anm[i])) anm[i] else paste0("[", i, "]"),
        value = paste(deparse(a), collapse = " "),
        source = src, stringsAsFactors = FALSE)
    }
    invisible()
  })
})
e2 <- if (length(e2)) do.call(rbind, e2) else
  data.frame(file = character(0), line = integer(0), call = character(0),
             arg = character(0), value = character(0), source = character(0))

# --- E3a: knob arguments taking a literal -----------------------------------

KNOBS <- c("k", "bs", "method", "discrete", "select", "family", "nrounds",
           "nrounds_max", "eta", "max_depth", "subsample", "colsample_bytree",
           "min_child_weight", "lambda", "nthreads", "nthread", "gc.level",
           "gc_level", "early_stopping_rounds", "valid_frac", "n_bins", "n_boot",
           "sp", "gamma", "knots")

MODEL_CALLEES <- c(FITTERS, "s", "te", "ti", "optim", "nlm", "glm", "lm",
                   "predict", "predict.gam", "predict.bam", "polr")

e3 <- list()
for_each_top(function(file, line, owner, top) {
  walk_ast(top, function(n) {
    if (!is.call(n)) return(invisible())
    fn <- callee_of(n)
    args <- as.list(n)[-1]
    anm <- names(args)
    if (is.null(anm)) return(invisible())
    for (i in seq_along(args)) {
      if (!nzchar(anm[i]) || !(anm[i] %in% KNOBS)) next
      if (identical(unname(args[i]), list(quote(expr = )))) next
      if (!is_literal(args[[i]])) next
      e3[[length(e3) + 1L]] <<- data.frame(
        knob = anm[i], callee = fn, file = file, line = line,
        value = fmt(args[[i]]), model_setting = fn %in% MODEL_CALLEES,
        stringsAsFactors = FALSE)
    }
    invisible()
  })
})
e3 <- if (length(e3)) do.call(rbind, e3) else
  data.frame(knob = character(0), callee = character(0), file = character(0),
             line = integer(0), value = character(0), model_setting = logical(0))

# --- E3b: smooth specifications built as strings -----------------------------

src <- scan_sources(root)
lib <- src[src$layer == "lib", ]
hit <- grepl("s\\s*\\(.*(bs\\s*=|k\\s*=)", lib$code, perl = TRUE) &
       grepl("[\"']", lib$code, perl = TRUE)
e3b <- data.frame(file = basename(lib$path[hit]), line = lib$line[hit],
                  code = trimws(lib$code[hit]), stringsAsFactors = FALSE)
# A literal digit after `k =` inside a smooth spec is a basis dimension written
# into the code; a `%d` is the config-driven form that `smooth_k_of()` fills.
e3b$literal_k <- grepl("k\\s*=\\s*[0-9]", e3b$code, perl = TRUE)

# --- E3c: formal defaults that restate a config value ------------------------

flat_cfg <- local({
  out <- list()
  add <- function(prefix, v) {
    for (nm in names(v)) {
      val <- v[[nm]]
      key <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
      if (is.list(val) && !is.null(names(val))) add(key, val) else out[[key]] <<- val
    }
  }
  add("", cfg)
  out
})
leaf_of <- function(k) sub("^.*[.]", "", k)

e3c <- list()
for_each_top(function(file, line, owner, top) {
  walk_ast(top, function(n) {
    if (!is.call(n) || !identical(.fname(n), "function") || length(n) < 2L) return(invisible())
    fm <- n[[2]]
    if (!is.pairlist(fm) && !is.list(fm)) return(invisible())
    l <- as.list(fm)
    for (i in seq_along(l)) {
      nm <- names(l)[i]
      if (is.null(nm) || !nzchar(nm)) next
      if (identical(unname(l[i]), list(quote(expr = )))) next
      d <- l[[i]]
      if (!is_literal(d)) next
      keys <- names(flat_cfg)[vapply(names(flat_cfg), function(k)
        identical(leaf_of(k), nm), logical(1))]
      if (!length(keys)) next
      for (k in keys) {
        e3c[[length(e3c) + 1L]] <<- data.frame(
          file = file, line = line, owner = owner, formal = nm,
          default = fmt(d), key = k, declared = fmt(flat_cfg[[k]]),
          verdict = if (identical(fmt(d), fmt(flat_cfg[[k]]))) "agree" else "DISAGREE",
          stringsAsFactors = FALSE)
      }
    }
    invisible()
  })
})
e3c <- if (length(e3c)) unique(do.call(rbind, e3c)) else
  data.frame(file = character(0), line = integer(0), owner = character(0),
             formal = character(0), default = character(0), key = character(0),
             declared = character(0), verdict = character(0))

# --- report ------------------------------------------------------------------

bad1 <- e1[e1$verdict != "agree", , drop = FALSE]
lit2 <- e2[e2$source == "LITERAL", , drop = FALSE]
bad3c <- e3c[e3c$verdict != "agree", , drop = FALSE]

md <- c("# Audit E — fit-setting provenance", "",
        "Generated by `tests/audit_e_fitsettings.R`. Source and config text only.", "",
        "`line` is the START LINE OF THE ENCLOSING TOP-LEVEL DEFINITION, not the",
        "exact line of the hit: R keeps source references on top-level expressions",
        "and not on the calls inside them. Use it with `owner` to find the hit.", "")

md <- c(md, "## E1 — `%||%` fallbacks against the declared config value", "",
        "The codebase reaches every knob as `cfg$group$key %||% <literal>`, so",
        "every knob has two declarations. `DISAGREE` means they differ and which",
        "is in force depends on whether the config key is present. `no key` means",
        "the fallback is the only declaration.", "",
        sprintf("%d fallback(s): %d agree, %d disagree, %d with no config key.",
                nrow(e1), sum(e1$verdict == "agree"),
                sum(e1$verdict == "DISAGREE"), sum(e1$verdict == "no key")), "",
        if (nrow(bad1)) md_table(bad1[order(bad1$verdict, bad1$key), ]) else
          "**PASS** — every fallback restates a declared config value.", "")

md <- c(md, "### fallbacks that agree (reference)", "",
        if (any(e1$verdict == "agree"))
          md_table(e1[e1$verdict == "agree", c("file", "line", "key", "fallback")])
        else "_none_", "")

md <- c(md, "", "## E2 — arguments of every fitter call", "",
        sprintf("%d argument(s) across %d call site(s); %d literal.",
                nrow(e2), length(unique(paste(e2$file, e2$line))), nrow(lit2)), "",
        "### literal arguments", "",
        if (nrow(lit2)) md_table(lit2[, c("file", "line", "call", "arg", "value")])
        else "**PASS** — no fitter argument is a bare literal.", "",
        "### all fitter arguments", "",
        if (nrow(e2)) md_table(e2[, c("file", "line", "call", "arg", "value", "source")])
        else "_none_", "")

md <- c(md, "", "## E3a — knob arguments taking a literal, by callee", "",
        "`model_setting` marks a callee that fits or predicts. A",
        "`method = \"spearman\"` on `stats::cor` matches the knob name and is not",
        "a fit setting.", "",
        if (nrow(e3)) md_table(e3[order(!e3$model_setting, e3$knob, e3$file), ])
        else "_none_", "")

md <- c(md, "## E3b — smooth specifications built as strings", "",
        "`R/05_formula.R` builds every model smooth with `sprintf()`, so `bs` and",
        "`k` are never arguments and no AST pass reaches them. `literal_k = TRUE`",
        "is a basis dimension written into the code rather than taken from",
        "`config/smooth_k`.", "",
        if (nrow(e3b)) md_table(e3b[order(!e3b$literal_k, e3b$file, e3b$line), ])
        else "_none_", "")

md <- c(md, "## E3c — formal defaults that restate a config value", "",
        "MATCHED BY LEAF NAME, so every row needs adjudicating: a formal called",
        "`n_bins` is not necessarily the config's `metrics.n_bins`, and a",
        "`seed = 1L` default is a placeholder that every call site overrides with",
        "`cfg$seed`. What the check is for is the case where the default IS the",
        "value in force at some entry point, which is invisible from either",
        "declaration alone.", "",
        sprintf("%d default(s) share a name with a config key: %d agree, %d disagree.",
                nrow(e3c), sum(e3c$verdict == "agree"), nrow(bad3c)), "",
        if (nrow(bad3c)) md_table(bad3c) else
          "**PASS** — every such default restates the declared value.", "",
        "### defaults that agree (reference)", "",
        if (any(e3c$verdict == "agree"))
          md_table(e3c[e3c$verdict == "agree", c("file", "line", "owner", "formal",
                                                 "default", "key")])
        else "_none_", "")

audit_write(md, "audit_e_fitsettings.md", root)

cat("\n=== AUDIT E: fit-setting provenance ===\n")
cat(sprintf("\nE1  %d fallbacks: %d agree, %d DISAGREE, %d no config key\n",
            nrow(e1), sum(e1$verdict == "agree"), sum(e1$verdict == "DISAGREE"),
            sum(e1$verdict == "no key")))
if (nrow(bad1)) print(bad1[order(bad1$verdict, bad1$key), ], row.names = FALSE)
cat("\nE2  fitter-call arguments that are bare literals:\n")
if (nrow(lit2)) print(lit2[, c("file", "line", "call", "arg", "value")], row.names = FALSE) else cat("  none\n")
cat("\nE3a knob literals on a fitting/predicting callee:\n")
if (any(e3$model_setting)) print(e3[e3$model_setting, c("knob", "callee", "file", "line", "value")],
                                 row.names = FALSE) else cat("  none\n")
cat("\nE3a knob literals on other callees (expected, listed for completeness): ",
    sum(!e3$model_setting), " site(s)\n", sep = "")
cat("\nE3b smooth specifications built as strings:\n")
if (nrow(e3b)) print(e3b[order(!e3b$literal_k, e3b$file, e3b$line), ], row.names = FALSE) else cat("  none\n")
cat("\nE3c formal defaults that restate a config value: ", nrow(e3c),
    " (", nrow(bad3c), " disagree)\n", sep = "")
if (nrow(bad3c)) print(bad3c, row.names = FALSE)
