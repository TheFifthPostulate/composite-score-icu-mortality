# tests/audit_d_columns.R ----------------------------------------------------
# AUDIT D of docs/v2_audit_plan_20260903.md section 2: hard-coded column lists
# against `CONTRACT`.
#
# THE QUESTION IT EXISTS FOR. `.measured_train_rows()` held a literal `need`
# vector and then wrote `sf[keep, need]`. Both halves were correct in
# isolation: the vector named real columns and the subset took exactly those.
# The defect existed only relative to a CALLER that expected other columns to
# survive, and when the `se_` replicate statistics were added to the extraction
# they reached the parquet, passed validation, and never reached the estimator.
# Nobody reading the function could see it. The plan's question is how many
# siblings that vector has.
#
# WHAT IS COLLECTED. Every literal `c("a", "b", ...)` of length three or more,
# all of whose elements are character, anywhere under `R/`. Each is scored
# against every table in `CONTRACT` and attributed to the best match.
#
# THE CLASSIFICATION THAT MATTERS. A literal column list is one of two things,
# and they have opposite risk:
#
#   REQUIREMENT   used only as `setdiff(need, names(df))` -- an assertion that
#                 the columns are present. Adding a column to the extraction
#                 changes nothing and nothing is lost. Safe.
#   SUBSET        used inside `[` -- `df[, need]`. Adding a column to the
#                 extraction silently drops it here. This is the failure class.
#
# The audit decides between them mechanically: for each variable bound to a
# literal character vector, it asks whether that variable ever appears in an
# index position of a `[` call within the same function body. A literal used
# inline inside `[` is caught the same way.
#
# The omitted-columns column is the payload: for a SUBSET hit against a
# CONTRACT table, it names exactly which contract columns do NOT survive.
#
# READS NO DATA. `CONTRACT` is a declaration in `R/02_validate.R`, and this
# audit reads source text (hard rule 1).
#
#   Rscript tests/audit_d_columns.R        ~5 s, writes out/audit_d_columns.md
# ----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  root <- "."
  source(file.path(root, "tests", "audit_common.R"))
  source(file.path(root, "R", "02_validate.R"))
}))

contract_cols <- lapply(CONTRACT, as.character)
all_contract  <- unique(unlist(contract_cols, use.names = FALSE))

MIN_LEN <- 3L   # the plan's threshold: "length three or more"

# --- AST helpers -------------------------------------------------------------

.is_char_vec_call <- function(x) {
  if (!is.call(x)) return(FALSE)
  fn <- tryCatch(as.character(x[[1]])[1], error = function(e) "")
  if (!identical(fn, "c")) return(FALSE)
  args <- as.list(x)[-1]
  if (length(args) < MIN_LEN) return(FALSE)
  all(vapply(args, function(a) is.character(a) && length(a) == 1L, logical(1)))
}

.char_vec_value <- function(x) as.character(unlist(as.list(x)[-1]))

#' Walk an expression, calling `f(node)` on every node. `l[i]` guards the
#' empty-argument slot in `z[, drop = FALSE]`, which throws when forced.
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

#' Symbols that appear in an INDEX position of a `[` or `[[` call, plus any
#' literal character vector used there. This is what separates a subset from a
#' requirement check.
index_symbols <- function(expr) {
  out <- character(0)
  walk_ast(expr, function(n) {
    if (!is.call(n)) return(invisible())
    fn <- tryCatch(as.character(n[[1]])[1], error = function(e) "")
    if (!fn %in% c("[", "[[")) return(invisible())
    args <- as.list(n)[-c(1, 2)]           # drop `[` and the object
    for (i in seq_along(args)) {
      if (identical(unname(args[i]), list(quote(expr = )))) next
      out <<- c(out, all.vars(args[[i]]))
    }
    invisible()
  })
  unique(out)
}

#' Every `[` index expression in `expr`, deparsed, paired with the variables it
#' mentions. This is what makes a SUBSET hit self-adjudicating: `sf[keep, need]`
#' and `sf[keep, c(need, opt)]` both bind `need` to an index position, and only
#' the first one drops what `need` omits. Printing the expression settles it
#' without anybody opening the file.
index_exprs <- function(expr) {
  out <- list()
  walk_ast(expr, function(n) {
    if (!is.call(n)) return(invisible())
    fn <- tryCatch(as.character(n[[1]])[1], error = function(e) "")
    if (!fn %in% c("[", "[[")) return(invisible())
    args <- as.list(n)[-c(1, 2)]
    for (i in seq_along(args)) {
      if (identical(unname(args[i]), list(quote(expr = )))) next
      out[[length(out) + 1L]] <<- list(
        txt  = paste(deparse(args[[i]]), collapse = " "),
        vars = all.vars(args[[i]]))
    }
    invisible()
  })
  out
}

#' Literal character vectors used directly inside a `[` index position.
index_literals <- function(expr) {
  out <- list()
  walk_ast(expr, function(n) {
    if (!is.call(n)) return(invisible())
    fn <- tryCatch(as.character(n[[1]])[1], error = function(e) "")
    if (!fn %in% c("[", "[[")) return(invisible())
    args <- as.list(n)[-c(1, 2)]
    for (i in seq_along(args)) {
      if (identical(unname(args[i]), list(quote(expr = )))) next
      walk_ast(args[[i]], function(m) {
        if (.is_char_vec_call(m)) out[[length(out) + 1L]] <<- .char_vec_value(m)
        invisible()
      })
    }
    invisible()
  })
  out
}

# --- collect -----------------------------------------------------------------

rows <- list()
files <- sort(list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE))
for (f in files) {
  ex <- parse(f, keep.source = TRUE)
  refs <- utils::getSrcref(ex)
  for (j in seq_along(ex)) {
    top  <- ex[[j]]
    line <- if (!is.null(refs) && length(refs) >= j) as.integer(refs[[j]])[1] else NA_integer_
    owner <- if (is.call(top) && as.character(top[[1]])[1] %in% c("<-", "=", "<<-") &&
                 (is.symbol(top[[2]]) || is.character(top[[2]])))
      as.character(top[[2]]) else "<top level>"

    idx_syms <- index_symbols(top)
    idx_lits <- index_literals(top)
    idx_exprs <- index_exprs(top)

    walk_ast(top, function(n) {
      if (!.is_char_vec_call(n)) return(invisible())
      v <- .char_vec_value(n)
      hit <- sum(v %in% all_contract)
      if (hit < 2L) return(invisible())

      # Which contract table does it look like?
      ov <- vapply(contract_cols, function(cc) sum(v %in% cc), integer(1))
      tbl <- names(ov)[which.max(ov)]
      cc  <- contract_cols[[tbl]]

      used_as_index <- any(vapply(idx_lits, function(z) identical(z, v), logical(1)))
      rows[[length(rows) + 1L]] <<- list(
        file = basename(f), line = line, owner = owner,
        n = length(v), n_contract = hit, table = tbl,
        pct_of_table = round(sum(v %in% cc) / length(cc), 3),
        literal_index = used_as_index,
        omitted = paste(setdiff(cc, v), collapse = " "),
        value = paste(v, collapse = " "))
      invisible()
    })

    # Variables bound to a literal character vector inside this top-level
    # expression, and whether that variable is used as an index anywhere in it.
    walk_ast(top, function(n) {
      if (!is.call(n)) return(invisible())
      fn <- tryCatch(as.character(n[[1]])[1], error = function(e) "")
      if (!fn %in% c("<-", "=", "<<-")) return(invisible())
      if (!is.symbol(n[[2]])) return(invisible())
      rhs <- n[[3]]
      if (!.is_char_vec_call(rhs)) return(invisible())
      v <- .char_vec_value(rhs)
      if (sum(v %in% all_contract) < 2L) return(invisible())
      nm <- as.character(n[[2]])
      for (k in seq_along(rows)) {
        if (identical(rows[[k]]$value, paste(v, collapse = " ")) &&
            identical(rows[[k]]$owner, owner) && identical(rows[[k]]$file, basename(f))) {
          rows[[k]]$bound_to <<- nm
          rows[[k]]$var_index <<- nm %in% idx_syms
          rows[[k]]$index_expr <<- paste(unique(vapply(
            Filter(function(z) nm %in% z$vars, idx_exprs),
            function(z) z$txt, character(1))), collapse = " | ")
        }
      }
      invisible()
    })
  }
}

d <- do.call(rbind, lapply(rows, function(r) data.frame(
  file = r$file, line = r$line, owner = r$owner,
  n = r$n, n_contract = r$n_contract, table = r$table,
  pct_of_table = r$pct_of_table,
  bound_to = r$bound_to %||% "",
  index_expr = r$index_expr %||% "",
  role = if (isTRUE(r$literal_index) || isTRUE(r$var_index)) "SUBSET" else "requirement",
  omitted = r$omitted, value = r$value,
  stringsAsFactors = FALSE)))
if (is.null(d)) d <- data.frame()

subsets <- if (nrow(d)) d[d$role == "SUBSET", , drop = FALSE] else d
reqs    <- if (nrow(d)) d[d$role != "SUBSET", , drop = FALSE] else d

# --- report ------------------------------------------------------------------

md <- c("# Audit D — hard-coded column lists against `CONTRACT`", "",
        "Generated by `tests/audit_d_columns.R`. Every literal `c(\"a\", \"b\", ...)`",
        sprintf("of length >= %d under `R/` with at least two elements naming a", MIN_LEN),
        "`CONTRACT` column, attributed to the contract table it best matches.", "",
        "`line` is the START LINE OF THE ENCLOSING TOP-LEVEL DEFINITION, not the",
        "exact line of the hit: R keeps source references on top-level expressions",
        "and not on the calls inside them. Use it with `owner` to find the hit.", "",
        "**SUBSET** means the list is used inside `[` — adding a column to the",
        "extraction silently drops it here, which is the `.measured_train_rows()`",
        "failure. **requirement** means it is used only to assert presence, which",
        "is safe. `omitted` names the contract columns of that table which the list",
        "does NOT carry.", "",
        sprintf("%d list(s) found: %d SUBSET, %d requirement.",
                nrow(d), nrow(subsets), nrow(reqs)), "")

md <- c(md, "## SUBSET — a list that silently drops what it does not name", "",
        if (nrow(subsets)) md_table(subsets[, c("file", "line", "owner", "bound_to",
                                                "index_expr", "table", "n",
                                                "n_contract", "omitted")])
        else "_none_", "")

md <- c(md, "## requirement — presence assertions, safe", "",
        if (nrow(reqs)) md_table(reqs[, c("file", "line", "owner", "table",
                                          "n", "n_contract", "value")])
        else "_none_", "")

audit_write(md, "audit_d_columns.md", root)

cat("\n=== AUDIT D: hard-coded column lists ===\n")
cat(sprintf("%d list(s): %d SUBSET, %d requirement\n\n",
            nrow(d), nrow(subsets), nrow(reqs)))
cat("SUBSET (silently drops what it does not name):\n")
if (nrow(subsets)) {
  for (i in seq_len(nrow(subsets))) {
    cat("\n  ", subsets$file[i], ":", subsets$line[i], "  ", subsets$owner[i], "\n", sep = "")
    cat("    list       ", subsets$value[i], "\n", sep = "")
    cat("    index expr ",
        if (nzchar(subsets$index_expr[i])) subsets$index_expr[i] else "<inline literal>",
        "\n", sep = "")
    cat("    omits from ", subsets$table[i], ": ",
        if (nzchar(subsets$omitted[i])) subsets$omitted[i] else "(nothing)", "\n", sep = "")
  }
} else cat("  none\n")
cat("\nrequirement lists:\n")
if (nrow(reqs)) print(reqs[, c("file", "line", "owner", "table", "n", "n_contract")],
                      row.names = FALSE) else cat("  none\n")
