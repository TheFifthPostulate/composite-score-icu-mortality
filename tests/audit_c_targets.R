# tests/audit_c_targets.R ----------------------------------------------------
# AUDIT C of docs/v2_audit_plan_20260903.md section 2: target dependency
# hygiene, plus the two hard rules that are otherwise enforced by convention
# alone.
#
# THE QUESTION IT EXISTS FOR. `mimic_files` was tracked as a `format = "file"`
# target and nothing depended on it, so re-extracting a parquet rebuilt
# `mimic_files`, changed nothing downstream, and `tabs` served a cached table
# that predated the new columns. The graph was not wrong in any line; it was
# wrong in an EDGE THAT WAS ABSENT. This audit reads the edge set as data.
#
#   C1   Which targets have no dependents? A `format = "file"` target with no
#        dependents is the `mimic_files` signature and is a hard failure: it
#        exists solely to invalidate something and it invalidates nothing.
#   C1b  A non-file target with no dependents is a LEAF, which is legitimate
#        only if a runner exports it. A leaf that is neither consumed nor
#        exported is computed for nothing. Checked against `config/internal.yml`
#        in both directions, so an export naming a target that does not exist
#        also fails here rather than at export time.
#   C2   Does any target's command reference a name the graph does not carry as
#        a dependency? Normally empty; non-empty exactly when a name is reached
#        in a way static analysis cannot see, which is when a stale value gets
#        served.
#   C3   HARD RULE 7, made mechanical. No target's VALUE may contain a
#        timestamp. Checked over the TRANSITIVE CLOSURE of functions the
#        command can reach, because a clock two calls deep is still a clock.
#        Unseeded RNG is checked in the same pass and for the same reason: a
#        draw outside `with_seed()` makes a value change between runs exactly
#        the way a clock does. A draw INSIDE `with_seed()` is not reported --
#        the walker refuses to descend into one.
#   C4   HARD RULE 9, made mechanical. Only `R/11_run.R` may build a path or
#        read the clock, and no target may reach a function defined there.
#   C5   Is any function defined in more than one file under `R/`? `tar_source`
#        sources in filename order and the last definition wins, so a duplicate
#        makes behaviour depend on source order and makes every file-based
#        attribution ambiguous.
#
# READS THE GRAPH, NOT THE STORE. `tar_network()` and `tar_manifest()` are
# metadata operations; no target value is loaded (hard rule 1).
#
#   Rscript tests/audit_c_targets.R        ~20 s, writes out/audit_c_targets.md
# ----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  root <- "."
  source(file.path(root, "tests", "audit_common.R"))
  library(targets)
}))

net  <- targets::tar_network(targets_only = FALSE, callr_function = NULL)
man  <- targets::tar_manifest(fields = c("name", "command", "pattern", "format"),
                              callr_function = NULL)
edge <- as.data.frame(net$edges, stringsAsFactors = FALSE)
vert <- as.data.frame(net$vertices, stringsAsFactors = FALSE)

tgt_names <- man$name

# --- expression views --------------------------------------------------------
# `.vars()` (all.vars) is the dependency view: `cfg$paths` yields `cfg`, not
# `paths`, because a list field name is not a dependency and counting it as one
# manufactures false findings. `.names()` (all.names) is the call view, which is
# what "does this reach Sys.time" needs.

.parse1 <- function(txt) tryCatch(parse(text = txt), error = function(z) NULL)
.vars <- function(txt) {
  e <- .parse1(txt); if (is.null(e)) character(0)
  else unique(unlist(lapply(as.list(e), all.vars)))
}
.names <- function(txt) {
  e <- .parse1(txt); if (is.null(e)) character(0)
  else unique(unlist(lapply(as.list(e), all.names)))
}

#' Names from `needles` that occur OUTSIDE any call to a `guard` function.
#'
#' The walker refuses to descend into a guarded call, which is how a seeded
#' `sample.int()` is distinguished from an unseeded one without the audit
#' having to know anything about seeding.
#'
#' `l[i]` rather than `l[[i]]` for the emptiness test: a call with an empty
#' argument slot (`z[, drop = FALSE]`) holds R's missing-argument object there,
#' and passing that to a function throws.
collect_unguarded <- function(x, needles, guard = character(0), inside = FALSE) {
  out <- character(0)
  if (is.symbol(x)) {
    s <- as.character(x)
    if (nzchar(s) && !inside && s %in% needles) out <- s
    return(out)
  }
  if (is.call(x)) {
    fn <- tryCatch(as.character(x[[1]]), error = function(e) character(0))
    if (length(fn) && fn[1] %in% guard) inside <- TRUE
  }
  if (is.call(x) || is.pairlist(x) || is.expression(x) || is.list(x)) {
    l <- as.list(x)
    for (i in seq_along(l)) {
      if (identical(unname(l[i]), list(quote(expr = )))) next
      out <- c(out, collect_unguarded(l[[i]], needles, guard, inside))
    }
  }
  unique(out)
}

unguarded_in <- function(txt, needles, guard = character(0)) {
  e <- .parse1(txt)
  if (is.null(e)) return(character(0))
  unique(unlist(lapply(as.list(e), collect_unguarded,
                       needles = needles, guard = guard)))
}

# --- C1 / C1b: dependents and exports ---------------------------------------

n_dependents <- vapply(tgt_names, function(n) sum(edge$from == n), integer(1))
n_deps       <- vapply(tgt_names, function(n) sum(edge$to   == n), integer(1))

int <- yaml::read_yaml(file.path(root, "config", "internal.yml"))
exported <- unique(c(unlist(int$export$tables), unlist(int$export$diagnostics),
                     unlist(int$export$objects), unlist(int$export$bundle)))
exported <- as.character(exported)

c1 <- data.frame(target = tgt_names, format = man$format,
                 n_deps = as.integer(n_deps),
                 n_dependents = as.integer(n_dependents),
                 exported = tgt_names %in% exported,
                 stringsAsFactors = FALSE)
c1$leaf <- c1$n_dependents == 0L
c1$verdict <- ifelse(c1$leaf & c1$format == "file", "FAIL: file target, no dependents",
              ifelse(c1$leaf & !c1$exported, "FAIL: leaf, not exported",
              ifelse(c1$leaf, "leaf, exported", "ok")))

c1_fail <- c1[grepl("^FAIL", c1$verdict), , drop = FALSE]
c1_leaf <- c1[c1$verdict == "leaf, exported", , drop = FALSE]
export_not_target <- setdiff(exported, tgt_names)

# --- C2 ---------------------------------------------------------------------

env <- new.env(parent = globalenv())
for (f in list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE)) {
  sys.source(f, envir = env)
}
globals <- ls(env, all.names = TRUE)

c2 <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
  nm  <- man$name[i]
  sym <- .vars(man$command[i])
  # A pattern's `map(x)` argument is a dependency expressed in `pattern`, not
  # in `command`, so it is added here rather than reported as missing.
  sym <- c(sym, if (!is.na(man$pattern[i])) .vars(man$pattern[i]) else character(0))
  want <- setdiff(sym[sym %in% tgt_names | sym %in% globals], nm)
  miss <- setdiff(want, edge$from[edge$to == nm])
  if (!length(miss)) return(NULL)
  data.frame(target = nm, missing = paste(sort(miss), collapse = ", "),
             stringsAsFactors = FALSE)
}))
if (is.null(c2)) c2 <- data.frame(target = character(0), missing = character(0))

# --- function -> file, the call graph, and C5 -------------------------------

fn_files <- list()   # name -> every file that defines it, in source order
fn_body  <- list()   # name -> the LAST definition, which is the one in force
for (f in sort(list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE))) {
  for (e in parse(f)) {
    if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=", "<<-") &&
        (is.symbol(e[[2]]) || is.character(e[[2]])) &&
        is.call(e[[3]]) && identical(as.character(e[[3]][[1]])[1], "function")) {
      nm <- as.character(e[[2]])
      fn_files[[nm]] <- c(fn_files[[nm]], basename(f))
      fn_body[[nm]]  <- paste(deparse(e[[3]]), collapse = "\n")
    }
  }
}
fn_file <- vapply(fn_files, function(z) z[length(z)], character(1))

c5 <- do.call(rbind, lapply(names(fn_files), function(nm) {
  if (length(fn_files[[nm]]) < 2L) return(NULL)
  data.frame(name = nm, files = paste(fn_files[[nm]], collapse = ", "),
             in_force = fn_files[[nm]][length(fn_files[[nm]])],
             stringsAsFactors = FALSE)
}))
if (is.null(c5)) c5 <- data.frame(name = character(0), files = character(0),
                                  in_force = character(0))

# --- C6: does the store's view of `R/` match the disk? ----------------------
#
# THE GAP THIS CLOSES, FOUND 2026-09-06 BY ADDING A FILE. `_targets.R` tracks
# the library as a `format = "file"` target whose command is a `list.files("R")`,
# which audit finding F3 introduced so that the contents of `R/` would be a
# DECLARED dependency of `bundle` rather than a hidden read. It does that job
# for the files it already knows about: edit one and `code_hashes` and `bundle`
# rebuild.
#
# IT DOES NOT NOTICE A FILE THAT DID NOT EXIST WHEN THE TARGET LAST RAN.
# `targets` decides whether a `format = "file"` target is stale by hashing the
# paths RECORDED IN ITS METADATA; it does not re-run the command to see whether
# the enumeration has changed. So adding `R/14_attribution_eval.R` left the
# metadata holding 20 paths against 21 on disk, `code_hashes` missing an entry,
# and `tar_outdated()` reporting NOTHING outdated -- while `tar_source("R")`
# loads the new file into every target's evaluation environment regardless.
#
# A new file that only adds functions is harmless, which is what made this easy
# to miss. One that REDEFINES an existing name silently changes what every
# target computes while the graph reports itself up to date. C5 above would name
# the collision, but only if someone ran it, and nothing would have told them to.
#
# There is no pure-`targets` fix: a directory listing cannot be a declared
# dependency of the command that produces it. So it is a check.
#
# AND `tar_make()` ALONE DOES NOT FIX IT, which is the part that matters when
# this check fires. The target is not outdated, so `tar_make()` skips it -- it
# reported "75 skipped" with a file on disk the store had never seen. The remedy
# is to invalidate the enumeration explicitly:
#
#   Rscript -e 'targets::tar_invalidate(source_files); targets::tar_make()'
#
# which rebuilds `source_files`, `code_hashes`, `bundle`, `bundle_contents` and
# `bundle_checks` in about 16 seconds and refits nothing.
#
# STILL METADATA, NOT A VALUE. `tar_meta()$path` holds a file target's tracked
# paths, so this reads the graph like everything else in the file and loads no
# target value (hard rule 1).
c6_disk <- sort(basename(list.files(file.path(root, "R"), pattern = "[.]R$")))
c6_store <- tryCatch(
  sort(basename(unlist(targets::tar_meta(names = "source_files")$path))),
  error = function(e) character(0))
c6 <- data.frame(file = sort(union(setdiff(c6_disk, c6_store),
                                   setdiff(c6_store, c6_disk))),
                 stringsAsFactors = FALSE)
c6$state <- if (nrow(c6)) ifelse(c6$file %in% c6_disk,
                                 "on disk, not in the store",
                                 "in the store, not on disk") else character(0)

call_cache <- new.env(parent = emptyenv())
closure_of <- function(seeds) {
  seen <- character(0)
  stack <- intersect(seeds, names(fn_body))
  while (length(stack)) {
    n <- stack[1]; stack <- stack[-1]
    if (n %in% seen) next
    seen <- c(seen, n)
    if (!exists(n, envir = call_cache, inherits = FALSE)) {
      assign(n, intersect(.names(fn_body[[n]]), names(fn_body)), envir = call_cache)
    }
    stack <- c(stack, setdiff(get(n, envir = call_cache), seen))
  }
  seen
}

# --- C3: hard rule 7 --------------------------------------------------------

# WHAT MAKES A TARGET'S VALUE UNSTABLE, in three separate classes, because
# they have different remedies and conflating them produces a false finding.
#
#   CLOCK    genuinely non-deterministic: the value differs between two runs of
#            identical code on an identical tree. Never acceptable in a target.
#   RNG      the same, unless the draw sits inside `with_seed()`.
#   FS_READ  reads the filesystem. NOT a defect on its own: a target may
#            legitimately read a file, and `targets` handles that correctly
#            PROVIDED the file is tracked as a `format = "file"` dependency.
#            The defect is an UNDECLARED read, which is what `build_bundle()`
#            calling `.source_hashes()` was (audit finding F3) and what
#            `code_hashes` reading `source_files` is not. Reported by C3b.
CLOCK <- c("Sys.time", "Sys.Date", "Sys.getenv", "date", "proc.time",
           "system.time", "Sys.setenv", "tempfile", "file.info")
RNG   <- c("runif", "rnorm", "sample", "sample.int", "rbinom", "rpois",
           "rmultinom", "rgamma", "rbeta")
FS_READ <- c("md5sum", "readRDS", "read.csv", "read_parquet", "readLines",
             "list.files", "file.exists", "qs_read")

# A function can be seeded by its CALLER rather than in its own body --
# `assign_folds()` writes `with_seed(cfg$seed, .allocate(...))`, so `.allocate`
# draws under a fixed stream even though nothing in `.allocate` mentions a
# seed. A purely lexical check would report it, so the guard is PROPAGATED:
# a function counts as seed-guarded when every call site of it, anywhere in
# `R/` or in a target command, sits inside a `with_seed()`. That is the
# condition under which its draws are in fact reproducible.
.seed_guarded <- local({
  callers <- c(unlist(fn_body, use.names = FALSE), man$command)
  fnames  <- names(fn_body)
  free <- unique(unlist(lapply(callers, function(b)
    unguarded_in(b, fnames, guard = "with_seed"))))
  seen <- unique(unlist(lapply(callers, function(b) intersect(.names(b), fnames))))
  setdiff(seen, free)
})

c3 <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
  nm    <- man$name[i]
  clos  <- closure_of(intersect(.names(man$command[i]), names(fn_body)))
  rows <- c(
    { h <- c(unguarded_in(man$command[i], CLOCK),
             unguarded_in(man$command[i], RNG, guard = "with_seed"))
      if (length(h)) paste0("<command>: ", paste(h, collapse = ",")) else NULL },
    unlist(lapply(clos, function(f) {
      h <- c(unguarded_in(fn_body[[f]], CLOCK),
             if (f %in% .seed_guarded) character(0)
             else unguarded_in(fn_body[[f]], RNG, guard = "with_seed"))
      if (length(h)) paste0(f, "() [", fn_file[f], "]: ", paste(h, collapse = ",")) else NULL
    })))
  if (!length(rows)) return(NULL)
  data.frame(target = nm, via = paste(rows, collapse = "; "),
             stringsAsFactors = FALSE)
}))
if (is.null(c3)) c3 <- data.frame(target = character(0), via = character(0))

# --- C4: hard rule 9 --------------------------------------------------------

# --- C3b: filesystem reads, and whether they are DECLARED --------------------
# A target that reads a file is fine when the file is tracked. `file_dep` is
# TRUE when the target has a `format = "file"` target among its transitive
# dependencies, which is what lets `targets` see the read at all.

file_targets <- man$name[man$format == "file"]
ancestors <- function(nm) {
  seen <- character(0); st <- edge$from[edge$to == nm]
  while (length(st)) {
    n <- st[1]; st <- st[-1]
    if (n %in% seen) next
    seen <- c(seen, n)
    st <- c(st, setdiff(edge$from[edge$to == n], seen))
  }
  seen
}

c3b <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
  nm   <- man$name[i]
  clos <- closure_of(intersect(.names(man$command[i]), names(fn_body)))
  rows <- unlist(lapply(clos, function(f) {
    h <- unguarded_in(fn_body[[f]], FS_READ)
    if (length(h)) paste0(f, "() [", fn_file[f], "]: ", paste(h, collapse = ",")) else NULL
  }))
  h0 <- unguarded_in(man$command[i], FS_READ)
  if (length(h0)) rows <- c(paste0("<command>: ", paste(h0, collapse = ",")), rows)
  if (!length(rows)) return(NULL)
  data.frame(target = nm,
             file_dep = length(intersect(c(nm, ancestors(nm)), file_targets)) > 0L,
             via = paste(rows, collapse = "; "), stringsAsFactors = FALSE)
}))
if (is.null(c3b)) c3b <- data.frame(target = character(0), file_dep = logical(0),
                                    via = character(0))

# HARD RULE 9 AS WRITTEN: BUILDING a path, or reading the clock. Opening a file
# whose path was passed in is not a breach -- that is exactly what R/01_load.R
# is for, and the rule's own wording says the caller passes paths in from
# config. Those calls are collected separately as context rather than counted
# as violations, because reporting them as breaches made the check cry wolf on
# the loader and buried the one file that genuinely broke the rule.
#
# `proc.time` and `system.time` ARE DELIBERATELY ABSENT from IO, and this is
# not an oversight to be tidied up later. Rule 9's subject is the clock: a
# datetime read in `R/` is what would let two apply sites stop being the same
# code path. Those two functions read process counters and return a DURATION,
# never a datetime, which is why `benchmark_fit()` (`R/06_layer1.R`) and
# `attribution_floor()` (`R/13b`) may time themselves. They stay in CLOCK
# above, because a duration in a TARGET's value is still unstable and that is
# hard rule 7, a different rule with a different remedy (audit finding F12).
IO <- c("file.path", "dir.create", "normalizePath", "list.files", "setwd",
        "getwd", "Sys.time", "Sys.Date")
IO_CONTEXT <- c("readRDS", "saveRDS", "read.csv", "write.csv", "writeLines",
                "readLines", "file.exists", "file.remove", "unlink", "sink",
                "png", "pdf", "md5sum", "read_parquet", "qs_read")

src <- scan_sources(root)
lib <- src[src$layer == "lib", ]
c4a <- do.call(rbind, lapply(unique(lib$path), function(p) {
  z <- lib[lib$path == p, ]
  hits <- unlist(lapply(IO, function(fn) {
    if (any(grepl(paste0("\\b", .lit(fn), "\\s*\\("), z$code, perl = TRUE))) fn else NULL
  }))
  if (!length(hits)) return(NULL)
  data.frame(file = basename(p), calls = paste(sort(hits), collapse = ", "),
             stringsAsFactors = FALSE)
}))
if (is.null(c4a)) c4a <- data.frame(file = character(0), calls = character(0))
c4a$exempt <- c4a$file == "11_run.R"

c4ctx <- do.call(rbind, lapply(unique(lib$path), function(p) {
  z <- lib[lib$path == p, ]
  hits <- unlist(lapply(IO_CONTEXT, function(fn) {
    if (any(grepl(paste0("\\b", .lit(fn), "\\s*\\("), z$code, perl = TRUE))) fn else NULL
  }))
  if (!length(hits)) return(NULL)
  data.frame(file = basename(p), calls = paste(sort(hits), collapse = ", "),
             stringsAsFactors = FALSE)
}))
if (is.null(c4ctx)) c4ctx <- data.frame(file = character(0), calls = character(0))

# Reaching a function DEFINED in R/11_run.R. `%||%` is excluded by name because
# it is defined in R/00_utils.R as well and the duplicate is reported by C5;
# attributing it to the run layer would bury every real hit under 27 copies of
# the same non-finding.
run_fns <- setdiff(names(fn_file)[fn_file == "11_run.R"],
                   names(fn_files)[lengths(fn_files) > 1L])
c4b <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
  hit <- intersect(closure_of(intersect(.names(man$command[i]), names(fn_body))),
                   run_fns)
  if (!length(hit)) return(NULL)
  data.frame(target = man$name[i], reaches = paste(sort(hit), collapse = ", "),
             stringsAsFactors = FALSE)
}))
if (is.null(c4b)) c4b <- data.frame(target = character(0), reaches = character(0))

# --- report -----------------------------------------------------------------

md <- c("# Audit C — target dependency hygiene, and hard rules 7 and 9", "",
        "Generated by `tests/audit_c_targets.R` from `tar_network()` and",
        "`tar_manifest()`. Metadata only; no target value is loaded.", "",
        sprintf("%d targets, %d graph vertices, %d edges.",
                nrow(man), nrow(vert), nrow(edge)), "")

md <- c(md, "## C1 — targets with no dependents", "",
        "A `format = \"file\"` target exists to invalidate something; one with no",
        "dependents invalidates nothing. A non-file leaf is legitimate only if a",
        "runner exports it.", "",
        if (nrow(c1_fail)) md_table(c1_fail[, c("target", "format", "n_deps", "verdict")])
        else "**PASS** — every leaf is either a file target with dependents or an exported leaf.",
        "")

md <- c(md, "### exported names that are not targets", "",
        if (length(export_not_target))
          paste0("**FAIL** — `", paste(export_not_target, collapse = "`, `"), "`")
        else "**PASS** — every name in `config/internal.yml` `export:` is a target.", "")

md <- c(md, "### leaves that are exported (expected)", "",
        if (nrow(c1_leaf)) paste0("`", paste(c1_leaf$target, collapse = "`, `"), "`") else "_none_", "")

md <- c(md, "", "## C2 — command references not carried as dependencies", "",
        if (nrow(c2)) md_table(c2)
        else "**PASS** — every name a command references is an incoming edge.", "")

md <- c(md, "## C3 — hard rule 7: clock or unseeded RNG reachable from a target value", "",
        "Transitive over the call graph. A draw inside `with_seed()` is not",
        "reported: the walker refuses to descend into one. The guard is also",
        "PROPAGATED — a function every one of whose call sites sits inside a",
        "`with_seed()` is treated as seeded, which is what `.allocate()` is.", "",
        if (nrow(c3)) md_table(c3)
        else "**PASS** — no target reaches a clock or an unseeded draw.", "")

md <- c(md, "## C3b — targets that read the filesystem", "",
        "Not a defect on its own. `file_dep = TRUE` means the target has a",
        "`format = \"file\"` target among its transitive dependencies, so",
        "`targets` can see the read. `file_dep = FALSE` is the",
        "`build_bundle()` / `.source_hashes()` failure: a hidden input.", "",
        if (nrow(c3b)) md_table(c3b) else "_none_", "")

md <- c(md, "## C4a — hard rule 9: path building and clock reads in `R/`", "",
        "`R/11_run.R` is the only exempt file. This is the rule as written:",
        "BUILDING a path, or reading the clock. Opening a file whose path was",
        "passed in is not a breach — that is what `R/01_load.R` is for — and",
        "those calls are listed below as context.", "",
        if (nrow(c4a)) md_table(c4a[order(!c4a$exempt, c4a$file), ]) else
          "**PASS** — no file under `R/` builds a path or reads the clock.", "",
        "### context: I/O on paths that were passed in", "",
        if (nrow(c4ctx)) md_table(c4ctx) else "_none_", "")

md <- c(md, "## C4b — targets reaching a function defined in `R/11_run.R`", "",
        if (nrow(c4b)) md_table(c4b)
        else "**PASS** — the graph never reaches the run/export layer.", "")

md <- c(md, "## C5 — functions defined in more than one file under `R/`", "",
        "`tar_source()` sources in filename order, so the last definition wins.", "",
        if (nrow(c5)) md_table(c5) else "**PASS** — every function has one definition.", "")

md <- c(md, "## C6 — the store's view of `R/` against the disk", "",
        "A `format = \"file\"` target is checked by hashing the paths recorded in",
        "its metadata, not by re-running its command, so `source_files` does not",
        "notice a file ADDED to `R/` since it last ran. `tar_source(\"R\")` loads",
        "that file into every target's environment regardless, and",
        "`tar_outdated()` reports nothing. `tar_make()` ALONE DOES NOT FIX IT --",
        "the target is not outdated, so it is skipped. The remedy is",
        "`targets::tar_invalidate(source_files); targets::tar_make()`, which",
        "rebuilds `source_files`, `code_hashes`, `bundle`, `bundle_contents` and",
        "`bundle_checks` in about 16 seconds and refits nothing.", "",
        if (nrow(c6)) md_table(c6) else
          "**PASS** — the store's file list matches `R/` on disk.", "")

audit_write(md, "audit_c_targets.md", root)

cat("\n=== AUDIT C: target hygiene ===\n")
cat(sprintf("%d targets, %d edges\n\n", nrow(man), nrow(edge)))
cat("C1  leaves that fail:\n")
if (nrow(c1_fail)) print(c1_fail[, c("target", "format", "n_deps", "verdict")], row.names = FALSE) else cat("  none\n")
cat("\nC1  exported names that are not targets: ",
    if (length(export_not_target)) paste(export_not_target, collapse = ", ") else "none", "\n", sep = "")
cat("\nC2  command names not carried as dependencies:\n")
if (nrow(c2)) print(c2, row.names = FALSE) else cat("  none\n")
cat("\nC3  clock / unseeded RNG reachable from a target value:\n")
if (nrow(c3)) print(c3, row.names = FALSE) else cat("  none\n")
cat("\nC3b targets that read the filesystem (file_dep FALSE is the defect):\n")
if (nrow(c3b)) print(c3b, row.names = FALSE) else cat("  none\n")
cat("\nC4a R/ files building paths or reading the clock (11_run.R exempt):\n")
if (nrow(c4a)) print(c4a[order(!c4a$exempt, c4a$file), ], row.names = FALSE) else cat("  none\n")
cat("\nC4b targets reaching R/11_run.R:\n")
if (nrow(c4b)) print(c4b, row.names = FALSE) else cat("  none\n")
cat("\nC5  duplicate function definitions:\n")
if (nrow(c5)) print(c5, row.names = FALSE) else cat("  none\n")
cat("\nC6  R/ files the store has not seen (run tar_make() to resynchronise):\n")
if (nrow(c6)) print(c6, row.names = FALSE) else cat("  none\n")
