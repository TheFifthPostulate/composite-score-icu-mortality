# run/reproduce.R ---------------------------------------------------------------
# Reproduce the paper: run the steps of config/reproduce.yml in order, hand each step the run directories its
# upstream steps created, and stop at the first failure (2026-11-13).
#
# WHAT A STEP IS. One or more commands (an interpreter named in the config, a script, its options), an
# optional run-directory prefix, optional runner-config edits and optional `needs`. `{id}` in an option is
# replaced by the run directory step `id` created in this chain. After a step's commands, exactly one new
# directory with the step's prefix must exist under the run root of R/11_run.R, and its manifest must say
# `status: complete`; that directory is the step's result. Steps without a prefix write caches or paper
# figures and are done when they exit 0.
#
# DEPENDENCIES. Before anything runs, every `{id}` of every selected step must have a run: done or bound in
# this chain, or produced by an earlier selected step. Otherwise the driver stops and says which. Steps also
# read what earlier steps leave without naming it (the targets store, the bundle path in the runner configs,
# out/cache/ebm, out/cache/jointgam, paper/figs/coefs); `needs` declares those, and when a partial selection
# skips one the driver prints a note and runs anyway, since a cache from earlier work may be what is meant.
#
# RESUMABLE. The state (which steps are done, the run each created, how many of a step's commands finished)
# is kept in <state_dir>/state.yml and written after every command. Rerunning the driver skips finished steps
# and finished commands, so a chain interrupted after ten hours resumes where it stopped. --redo forgets a
# step; its earlier run directory is kept, as every run is immutable. --bind records an existing complete run
# as a step's result; it is remembered until --redo forgets it. There is no "use the cited runs" shortcut:
# consumers of one producer cite different generations of it, so only an explicit --bind is unambiguous.
#
# CONFIG EDITS. The runners read the bundle and a few run paths from config/internal.yml, external.yml,
# external_demo.yml and gam_qc.yml by design. A step's `configs` set those keys to its run, and the driver
# keeps them in line with the state: they are set when the step finishes or is bound, and again whenever the
# driver starts with that step done (after --restore-configs too). Only the value on that line changes;
# comments, quoting and line endings stay. Each file is backed up once under <state_dir>/config_backup/
# before its first edit, and --restore-configs puts the backups back. Every change is printed.
#
# CONSOLE. Step names, run directories, config changes, durations and, on failure, the last lines of that
# step's log. The scripts' own output goes to <state_dir>/logs/; it is aggregate only by the scripts' contract
# (hard rule 1). No clock time is read or written (hard rule 10): durations come from start_timer().
#
#   Rscript run/reproduce.R --list                         # steps, tiers, state
#   Rscript run/reproduce.R --dry-run                      # what it would run, resolved; writes nothing
#   Rscript run/reproduce.R                                # everything not yet done, in order
#   Rscript run/reproduce.R --tier checks,core             # some tiers
#   Rscript run/reproduce.R --from gam_qc --to card        # a range, in list order
#   Rscript run/reproduce.R --only strata --bind internal=out/runs/internal_<stamp>   # name an upstream run
#   Rscript run/reproduce.R --redo gam_qc,card             # forget these steps, then run what is selected
#   Rscript run/reproduce.R --restore-configs              # put the runner configs back
#   options: --config <file> (default config/reproduce.yml); --bind may be repeated
# ------------------------------------------------------------------------------
if (!file.exists("R/11_run.R")) stop("run/reproduce.R runs from the project root", call. = FALSE)
suppressPackageStartupMessages(library(yaml))
source("R/00_utils.R")    # %||%, start_timer()
source("R/11_run.R")      # read_manifest(); the run root is new_run()'s default

ROOT  <- eval(formals(new_run)$root)
TIERS <- c("checks", "core", "comparators", "heterogeneity", "reproducibility", "figures")
DONE  <- c("done", "bound")

# --- arguments ----------------------------------------------------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
VALUED <- c("--config", "--tier", "--only", "--from", "--to", "--redo", "--bind")
SWITCH <- c("--list", "--dry-run", "--restore-configs")
.i <- 1L
while (.i <= length(args)) {
  a <- args[.i]
  if (a %in% VALUED) { if (.i == length(args)) stop(a, " needs a value", call. = FALSE); .i <- .i + 2L }
  else if (a %in% SWITCH) .i <- .i + 1L
  else stop("unknown argument ", a, " (see the header of run/reproduce.R)", call. = FALSE)
}
.opt <- function(flag) { i <- match(flag, args); if (is.na(i)) NULL else args[i + 1L] }
.csv <- function(x) if (is.null(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]
CONFIG <- .opt("--config") %||% "config/reproduce.yml"
DRY    <- "--dry-run" %in% args
BINDS  <- args[which(args == "--bind") + 1L]

M         <- yaml::read_yaml(CONFIG)
STATE_DIR <- M$state_dir %||% stop(CONFIG, " has no state_dir", call. = FALSE)
STATE_F   <- file.path(STATE_DIR, "state.yml")
BACKUP    <- file.path(STATE_DIR, "config_backup")
STEPS     <- M$steps
IDS       <- vapply(STEPS, function(s) as.character(s$id %||% ""), "")
names(STEPS) <- IDS
REF       <- "[{]([a-z0-9_]+)[}]"
refs_in   <- function(x) unique(gsub("[{}]", "", unlist(regmatches(x, gregexpr(REF, x)))))
step_text <- function(s) c(unlist(s$run), vapply(s$configs %||% list(), function(cf) as.character(cf$value %||% ""), ""))

# --- restore the runner configs and stop ----------------------------------------------------------------------
if ("--restore-configs" %in% args) {
  fs <- if (dir.exists(BACKUP)) list.files(BACKUP, recursive = TRUE) else character(0)
  if (!length(fs)) { cat("no backed-up runner config under", BACKUP, "\n"); quit(save = "no") }
  for (f in fs) {
    if (!file.copy(file.path(BACKUP, f), f, overwrite = TRUE, copy.date = TRUE)) stop("could not restore ", f, call. = FALSE)
    unlink(file.path(BACKUP, f)); cat("restored", f, "\n")
  }
  cat("the next run of the driver sets the edited keys to this chain's runs again\n")
  quit(save = "no")
}

# --- validate the config before running anything -------------------------------------------------------------
problems <- character(0)
bad <- function(...) problems <<- c(problems, paste0(...))
for (id in IDS[!grepl("^[a-z][a-z0-9_]*$", IDS)])
  bad("step id '", id, "': ids are lowercase letters, digits and _ (YAML reads y, n, yes, no, on, off as true/false)")
if (anyDuplicated(IDS)) bad("duplicated step ids: ", paste(unique(IDS[duplicated(IDS)]), collapse = ", "))
for (k in seq_along(STEPS)) {
  s <- STEPS[[k]]; id <- IDS[k]
  if (!isTRUE(s$tier %in% TIERS)) bad(id, ": tier must be one of ", paste(TIERS, collapse = ", "))
  if (!is.null(s$prefix) && !grepl("^[a-z0-9]+$", s$prefix)) bad(id, ": prefix must be lowercase alphanumeric")
  if (!length(s$run)) bad(id, ": no commands")
  for (cmd in unlist(s$run)) {
    tok <- strsplit(cmd, " +")[[1]]
    if (!tok[1] %in% names(M$interpreters)) bad(id, ": unknown interpreter '", tok[1], "'")
    if (length(tok) < 2L || !file.exists(tok[2])) bad(id, ": script not found: ", tok[2])
  }
  for (r in refs_in(step_text(s))) {
    j <- match(r, IDS)
    if (is.na(j)) bad(id, ": {", r, "} names no step")
    else if (j > k) bad(id, ": {", r, "} names a later step")
    else if (is.null(STEPS[[j]]$prefix)) bad(id, ": {", r, "} names a step that creates no run")
  }
  for (r in unlist(s$needs)) if (!isTRUE(match(r, IDS) < k)) bad(id, ": needs '", r, "', which is not an earlier step")
  for (cf in s$configs %||% list()) {
    if (is.null(cf$file) || is.null(cf$key) || is.null(cf$value)) bad(id, ": a config edit needs file, key and value")
    else if (!file.exists(cf$file)) bad(id, ": config file not found: ", cf$file)
  }
}
if (length(problems)) stop(CONFIG, " is not valid:\n  ", paste(problems, collapse = "\n  "), call. = FALSE)

# --- state: load, forget (--redo), bind (--bind) ---------------------------------------------------------------
STATE <- if (file.exists(STATE_F)) (yaml::read_yaml(STATE_F)$steps %||% list()) else list()
save_state <- function() {
  if (DRY) return(invisible())
  dir.create(STATE_DIR, showWarnings = FALSE, recursive = TRUE)
  yaml::write_yaml(list(config = CONFIG, steps = STATE), STATE_F)
}
is_done <- function(id) isTRUE(STATE[[id]]$status %in% DONE)
for (r in .csv(.opt("--redo"))) {
  if (!r %in% IDS) stop("--redo: no step ", r, call. = FALSE)
  if (!is.null(STATE[[r]]$run)) cat(sprintf("redo %s: forgetting %s (the directory is kept)\n", r, STATE[[r]]$run))
  STATE[[r]] <- NULL
}
for (b in BINDS) {
  kv <- strsplit(b, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L || !kv[1] %in% IDS) stop("--bind takes <step id>=<run directory>, got ", b, call. = FALSE)
  pfx <- STEPS[[kv[1]]]$prefix
  dir <- sub("[/\\\\]+$", "", kv[2])
  if (is.null(pfx)) stop("--bind ", kv[1], ": that step creates no run", call. = FALSE)
  if (!dir.exists(dir) || !grepl(paste0("^", pfx, "_[0-9]{8}T[0-9]{6}(-[0-9]+)?$"), basename(dir)))
    stop("--bind ", b, ": not an existing ", pfx, "_ run directory", call. = FALSE)
  if (!identical(tryCatch(read_manifest(dir)$status, error = function(e) NULL), "complete"))
    stop("--bind ", b, ": that run is not marked complete", call. = FALSE)
  STATE[[kv[1]]] <- list(status = "bound", run = dir)
}
if (length(BINDS) || !is.null(.opt("--redo"))) save_state()

# --- which steps ----------------------------------------------------------------------------------------------
pick <- function(x, what) { j <- match(x, IDS); if (anyNA(j)) stop(what, ": no step ", paste(x[is.na(j)], collapse = ", "), call. = FALSE); j }
SEL <- seq_along(STEPS)
if (!is.null(.opt("--tier"))) {
  tt <- .csv(.opt("--tier")); if (!all(tt %in% TIERS)) stop("--tier takes ", paste(TIERS, collapse = ", "), call. = FALSE)
  SEL <- SEL[vapply(STEPS[SEL], function(s) s$tier %in% tt, TRUE)]
}
if (!is.null(.opt("--from"))) SEL <- SEL[SEL >= pick(.opt("--from"), "--from")]
if (!is.null(.opt("--to")))   SEL <- SEL[SEL <= pick(.opt("--to"), "--to")]
if (!is.null(.opt("--only"))) SEL <- sort(pick(.csv(.opt("--only")), "--only"))

if ("--list" %in% args) {
  cat(sprintf("%-3s %-24s %-16s %-8s %s\n", "", "step", "tier", "state", "run"))
  for (k in seq_along(STEPS)) {
    st <- STATE[[IDS[k]]]
    cat(sprintf("%-3s %-24s %-16s %-8s %s\n", if (k %in% SEL) "*" else "", IDS[k], STEPS[[k]]$tier,
                st$status %||% "-", if (is.null(st$run)) "" else basename(st$run)))
  }
  quit(save = "no")
}

# --- interpreters, commands, run directories --------------------------------------------------------------------
interp_exe <- function(name) {
  exe <- M$interpreters[[name]]$exe %||% ""
  if (!nzchar(exe)) return(file.path(R.home("bin"), "Rscript"))
  if (!file.exists(exe)) stop(sprintf("interpreter '%s' not found at %s (interpreters in %s)", name, exe, CONFIG), call. = FALSE)
  normalizePath(exe, winslash = "/")
}
run_command <- function(name, argv, log) {
  it  <- M$interpreters[[name]]
  env <- unlist(it$env %||% list())
  if (length(env)) {
    old <- Sys.getenv(names(env), unset = NA, names = TRUE)     # one name is returned unnamed otherwise
    do.call(Sys.setenv, as.list(env))
    on.exit(for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n)))
  }
  # system2() quotes the command itself on Windows (quoting it again fails); the arguments are quoted here.
  system2(interp_exe(name), shQuote(c(unlist(it$args), argv)), stdout = log, stderr = log)
}
run_dirs <- function(prefix) {
  d <- list.dirs(ROOT, recursive = FALSE, full.names = TRUE)
  d[grepl(paste0("^", prefix, "_[0-9]{8}T[0-9]{6}(-[0-9]+)?$"), basename(d))]    # new_run() suffixes a clash
}
# `{id}` -> that step's run. In a dry run an unresolved one prints as <id> (made earlier in this selection) or
# <id: no run>; a real run never reaches this with one unresolved, because the check below stops it first.
resolve <- function(x, self, self_run, planned = character(0)) {
  for (r in refs_in(x)) {
    v <- if (identical(r, self)) self_run else if (is_done(r)) STATE[[r]]$run else NULL
    if (is.null(v)) {
      if (!DRY) stop(sprintf("step %s: {%s} has no run in this chain", self, r), call. = FALSE)
      v <- if (identical(r, self) || r %in% planned) paste0("<", r, ">") else paste0("<", r, ": no run>")
    }
    x <- gsub(paste0("{", r, "}"), v, x, fixed = TRUE)
  }
  x
}

# --- runner-config edits ---------------------------------------------------------------------------------------
set_config_key <- function(file, key, value, dry = FALSE) {
  parts <- strsplit(key, ".", fixed = TRUE)[[1]]
  get_key <- function() { x <- yaml::read_yaml(file); for (p in parts) x <- x[[p]]; x }
  if (identical(get_key(), value)) return(invisible(FALSE))
  raw <- readBin(file, "raw", file.size(file))
  eol <- if (any(raw == as.raw(13L))) "\r\n" else "\n"
  final_eol <- length(raw) > 0L && raw[length(raw)] == as.raw(10L)
  ln <- readLines(file, warn = FALSE, encoding = "UTF-8")
  if (length(parts) == 1L) {
    hit <- grep(paste0("^", parts, ":"), ln)
  } else if (length(parts) == 2L) {
    top <- grep(paste0("^", parts[1], ":[[:space:]]*(#.*)?$"), ln)
    if (length(top) != 1L) stop(file, ": no single block '", parts[1], ":'", call. = FALSE)
    end <- top + 1L
    while (end <= length(ln) && grepl("^([[:space:]]|#|$)", ln[end])) end <- end + 1L
    blk <- seq.int(top + 1L, length.out = end - top - 1L)
    hit <- blk[grepl(paste0("^[[:space:]]+", parts[2], ":"), ln[blk])]
    if (length(hit)) { ind <- nchar(sub("^([[:space:]]*).*", "\\1", ln[hit])); hit <- hit[ind == min(ind)] }
  } else stop("config keys are one or two levels deep: ", key, call. = FALSE)
  if (length(hit) != 1L) stop(file, ": key '", key, "' not found exactly once", call. = FALSE)
  LINE <- "^([^:]*:[[:space:]]*)([^#]*?)([[:space:]]*#.*)?$"
  old  <- sub(LINE, "\\2", ln[hit], perl = TRUE)
  q    <- substr(old, 1L, 1L); q <- if (q %in% c("\"", "'")) q else ""
  if (dry) { cat(sprintf("    config %s  %s: %s -> %s (dry run, not written)\n", file, key, old, value)); return(invisible(TRUE)) }
  bk <- file.path(BACKUP, file)
  if (!file.exists(bk)) {
    dir.create(dirname(bk), showWarnings = FALSE, recursive = TRUE)
    if (!file.copy(file, bk, copy.date = TRUE)) stop("could not back up ", file, call. = FALSE)
  }
  ln[hit] <- sub(LINE, paste0("\\1", q, gsub("\\", "\\\\", value, fixed = TRUE), q, "\\3"), ln[hit], perl = TRUE)
  writeBin(charToRaw(enc2utf8(paste0(paste(ln, collapse = eol), if (final_eol) eol else ""))), file)
  if (!identical(get_key(), value)) {
    writeBin(raw, file)
    stop(file, ": ", key, " did not read back as ", value, "; the file was left as it was", call. = FALSE)
  }
  cat(sprintf("    config %s  %s: %s -> %s\n", file, key, old, value))
  invisible(TRUE)
}
apply_configs <- function(k, dry = FALSE)
  for (cf in STEPS[[k]]$configs %||% list())
    set_config_key(cf$file, cf$key, resolve(cf$value, IDS[k], STATE[[IDS[k]]]$run), dry = dry)

# --- check the selection before running anything ---------------------------------------------------------------
cat(sprintf("%s: %d of %d steps selected; state in %s%s\n", CONFIG, length(SEL), length(STEPS), STATE_F,
            if (DRY) " (dry run: nothing is run or written)" else ""))
planned <- character(0); missing <- character(0); notes <- character(0); interps <- character(0)
for (k in SEL) {
  id <- IDS[k]; s <- STEPS[[k]]
  if (is_done(id)) next
  for (r in setdiff(refs_in(step_text(s)), id))
    if (!is_done(r) && !r %in% planned) missing <- c(missing, sprintf("%s uses the run of %s", id, r))
  for (r in unlist(s$needs))
    if (!is_done(r) && !r %in% planned) notes <- c(notes, sprintf("%s reads what %s leaves, and %s has not run in this chain", id, r, r))
  interps <- union(interps, vapply(unlist(s$run), function(cmd) strsplit(cmd, " +")[[1]][1], "", USE.NAMES = FALSE))
  planned <- c(planned, id)
}
no_exe <- Filter(function(n) inherits(try(interp_exe(n), silent = TRUE), "try-error"), interps)
if (length(no_exe)) missing <- c(missing, sprintf("interpreter '%s' not found at %s", no_exe,
                                                  vapply(no_exe, function(n) M$interpreters[[n]]$exe, "")))
if (length(notes)) cat(paste0("note: ", notes, "\n"), sep = "")
if (length(missing)) {
  msg <- paste0("  ", missing, collapse = "\n")
  if (!DRY) stop("cannot start:\n", msg, "\nselect those steps too, or pass --bind <step>=<run directory>", call. = FALSE)
  cat("missing (a real run would stop here):\n", msg, "\n", sep = "")
}

# --- bring the runner configs in line with the state (done and bound steps) ------------------------------------
for (k in seq_along(STEPS)) if (is_done(IDS[k]) && length(STEPS[[k]]$configs)) apply_configs(k, dry = DRY)

# --- run ---------------------------------------------------------------------------------------------------------
fmt_secs <- function(s) if (s < 120) paste(round(s), "s") else if (s < 7200) sprintf("%.1f min", s / 60) else sprintf("%.1f h", s / 3600)
fail <- function(id, msg, state, log = NULL) {
  STATE[[id]] <<- c(list(status = "failed"), state); save_state()
  cat("\n  FAILED: ", id, ": ", msg, "\n", sep = "")
  if (!is.null(log) && file.exists(log)) {
    cat("  log:", log, "\n"); tl <- utils::tail(readLines(log, warn = FALSE), 8L)
    if (length(tl)) cat(paste0("    | ", tl), sep = "\n")
  }
  quit(save = "no", status = 1L)
}
TOTAL <- start_timer()
planned <- character(0)
for (k in SEL) {
  s <- STEPS[[k]]; id <- IDS[k]; st <- STATE[[id]] %||% list()
  if (is_done(id)) { cat(sprintf("  %-24s %s %s\n", id, st$status, if (is.null(st$run)) "" else basename(st$run))); next }
  cmds <- unlist(s$run); n_done <- as.integer(st$commands_done %||% 0L); cur <- st$run
  if (DRY) {
    cat(sprintf("  %-24s %s%s\n", id, s$tier, if (identical(st$status, "failed")) ", failed earlier: resumes" else ""))
    for (j in seq_along(cmds))
      cat(sprintf("      %s%s\n", if (j <= n_done) "(done) " else "", resolve(cmds[j], id, cur, planned)))
    for (cf in s$configs %||% list())
      cat(sprintf("      then config %s  %s: %s\n", cf$file, cf$key, resolve(cf$value, id, cur, planned)))
    planned <- c(planned, id)
    next
  }
  dir.create(file.path(STATE_DIR, "logs"), showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("  %-24s %s\n", id,
              if (n_done >= length(cmds)) "its commands have run; checking its run again"
              else sprintf("running (%d command%s%s)", length(cmds), if (length(cmds) == 1L) "" else "s",
                           if (n_done > 0L) sprintf(", resuming at %d", n_done + 1L) else "")))
  t_step <- start_timer()
  for (j in seq_along(cmds)) {
    if (j <= n_done) next
    tok    <- strsplit(cmds[j], " +")[[1]]
    argv   <- vapply(tok[-1], resolve, "", self = id, self_run = cur, USE.NAMES = FALSE)
    log    <- file.path(STATE_DIR, "logs", if (length(cmds) == 1L) paste0(id, ".log") else sprintf("%s_%d.log", id, j))
    before <- if (is.null(s$prefix)) character(0) else run_dirs(s$prefix)
    status <- run_command(tok[1], argv, log)
    new    <- if (is.null(s$prefix)) character(0) else setdiff(run_dirs(s$prefix), before)
    if (!identical(as.integer(status), 0L))
      fail(id, sprintf("command %d of %d exited with status %s%s; rerunning the driver resumes at this command",
                       j, length(cmds), status,
                       if (length(new)) paste0(". It left an unfinished run directory, kept: ", paste(new, collapse = ", ")) else ""),
           list(run = cur, commands_done = j - 1L), log)
    if (length(new) > 1L || (length(new) == 1L && !is.null(cur)))
      fail(id, sprintf("command %d created %s, but a step creates exactly one %s_ run. Inspect them, then --bind %s=<dir> or --redo %s",
                       j, paste(c(cur, new), collapse = " and "), s$prefix, id, id),
           list(run = cur, commands_done = j - 1L), log)
    if (length(new) == 1L) cur <- new
    STATE[[id]] <- list(status = "running", run = cur, commands_done = j); save_state()
  }
  if (!is.null(s$prefix)) {
    if (is.null(cur))
      fail(id, sprintf("its commands finished but created no %s_ run. --redo %s runs it again from the first command", s$prefix, id),
           list(run = NULL, commands_done = length(cmds)))
    if (!identical(tryCatch(read_manifest(cur)$status, error = function(e) NULL), "complete"))
      fail(id, sprintf("%s is not marked complete. --redo %s runs it again, or --bind %s=<dir> names another run", cur, id, id),
           list(run = cur, commands_done = length(cmds)))
  }
  secs <- t_step()$elapsed_sec
  STATE[[id]] <- list(status = "done", run = cur, commands_done = length(cmds), seconds = round(secs)); save_state()
  cat(sprintf("  %-24s done in %s%s\n", id, fmt_secs(secs), if (is.null(cur)) "" else paste0(", ", cur)))
  apply_configs(k)
}
if (!DRY) cat(sprintf("finished: %d step(s) selected, %s\n", length(SEL), fmt_secs(TOTAL()$elapsed_sec)))
