# tests/audit_b_config.R -----------------------------------------------------
# AUDIT B of docs/v2_audit_plan_20260903.md section 2: config reachability.
#
# THE QUESTION. For every declaration in `config/config.yml`, is there any live
# code path that reads it? A declared control that nothing reads is worse than
# no declaration at all, because it reads as being in force. That is not a
# hypothetical: `signals_not_transportable` has been declared, documented and
# unread since 2026-08-29, and nothing in the pipeline noticed.
#
# THE THREE VERDICTS.
#   read       reached from R/, run/ or _targets.R -- a live path
#   vestigial  reached ONLY from tests/ -- a path no runner traverses. This is
#              where `smooth_k`'s constructs-off entries sit: the entries exist
#              for a design that only `tests/conditional_ablation.R` builds.
#   dead       reached from nowhere at all
#
# HOW A KEY IS ENUMERATED. Top-level keys, plus each child one level down --
# but ONLY for containers whose children are settings. A container keyed by
# DESIGN VOCABULARY (a signal name, an intervention name, a signal class) is
# a map, and the reachable thing is the map, not its rows: nobody writes
# `cfg$smooth_k$gcs_motor` in code and it would be a defect if they did. The
# test is mechanical -- child name in signals / interventions_extracted /
# signal_classes -- so it cannot drift out of date the way a hand-listed set of
# container keys would.
#
# WHY EVERY HIT IS ADJUDICATED RATHER THAN TRUSTED. `grep_name()` includes a
# bare-symbol arm, so a key named `seed` matches a local variable `seed` that
# has nothing to do with config. That over-matches on purpose: an over-matching
# reachability audit reports a dead key as live only when a same-named local
# exists, and it says so, whereas an under-matching one silently declares live
# keys dead and the audit becomes noise. The report therefore prints BOTH the
# permissive verdict and the strict one (explicit `$key` / `[["key"]]` /
# `"key"` only), and any key where the two disagree is called out.
#
# READS NO DATA. Config text and source text only (hard rule 1).
#
#   Rscript tests/audit_b_config.R          ~5 s, writes out/audit_b_config.md
# ----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  root <- "."
  source(file.path(root, "tests", "audit_common.R"))
}))

cfg <- yaml::read_yaml(file.path(root, "config", "config.yml"))
src <- scan_sources(root)

# `cfg_req()` / `cfg_flag()` reach a key by QUOTED ARGUMENT, not by `$` or
# `[[`, so the accessor matcher below cannot see them. Collected separately and
# counted as accessor evidence. Introduced 2026-09-03 with the fixes for audit
# findings F9, F10 and F11; the first run after those fixes reported three live
# keys as dead, which is the matcher failing safe.
req <- cfg_req_paths(root)

# --- enumerate the declarations ---------------------------------------------

# Names that are DESIGN VOCABULARY rather than settings. A container keyed by
# one of these is a map and is not descended into.
vocab <- unique(c(as.character(unlist(cfg$signals)),
                  as.character(unlist(cfg$interventions_extracted)),
                  as.character(unlist(cfg$signal_classes))))

keys <- list()
for (k in names(cfg)) {
  keys[[length(keys) + 1L]] <- list(key = k, path = k, depth = 1L)
  v <- cfg[[k]]
  if (is.list(v) && !is.null(names(v)) && length(v)) {
    kids <- names(v)
    if (!any(kids %in% vocab)) {
      for (kk in kids) {
        keys[[length(keys) + 1L]] <- list(key = kk, path = paste0(k, ".", kk),
                                          depth = 2L)
      }
    }
  }
}

# --- reachability ------------------------------------------------------------

#' Layers from which `key` is reached by a `cfg_req()`/`cfg_flag()` call.
#' Matched on the LAST path element, which is the key itself, and additionally
#' on the full dotted path so `diagnostics.decoupling_smooth` matches its own
#' children's calls.
req_layers <- function(path, key) {
  hit <- req$key == key | req$path == path | startsWith(req$path, paste0(path, "."))
  unique(req$layer[hit])
}

res <- do.call(rbind, lapply(keys, function(e) {
  arms  <- grep_arms(src, e$key)
  hit_p <- grep_name(src, e$key)          # permissive: includes bare symbol
  hit_s <- grep_name_strict(src, e$key)   # strict: $key / [["key"]] / "key"
  # An ACCESSOR hit is `$key` or `[["key"]]` -- the two forms that can only be
  # a lookup on a list. The quoted arm alone is not an accessor and is where
  # the false positives live.
  hit_a <- arms$dollar | arms$bracket
  lay_req <- req_layers(e$path, e$key)
  lay_p <- unique(c(src$layer[hit_p], lay_req))
  lay_s <- unique(c(src$layer[hit_s], lay_req))
  lay_a <- unique(c(src$layer[hit_a], lay_req))
  data.frame(
    path       = e$path,
    key        = e$key,
    depth      = e$depth,
    verdict    = reach_verdict(lay_a),
    strict     = reach_verdict(lay_s),
    permissive = reach_verdict(lay_p),
    n_accessor = sum(hit_a) + sum(req$key == e$key | req$path == e$path |
                                  startsWith(req$path, paste0(e$path, "."))),
    n_strict   = sum(hit_s),
    n_files    = length(unique(src$path[hit_a])),
    layers     = paste(sort(lay_a), collapse = "/"),
    stringsAsFactors = FALSE)
}))

# The three matchers are nested in strength, so a row where they disagree is a
# row where the evidence is ambiguous and a human must look.
res$disagree <- !(res$verdict == res$strict & res$strict == res$permissive)

# Where a key is reached, name the first live file so an adjudication is cheap.
res$first_live <- vapply(seq_len(nrow(res)), function(i) {
  a <- grep_arms(src, res$key[i])
  hit <- a$dollar | a$bracket
  h <- hit & src$layer != "tests"
  if (!any(h)) {
    rq <- req[(req$key == res$key[i] | req$path == res$path[i] |
               startsWith(req$path, paste0(res$path[i], "."))) & req$layer != "tests", ]
    if (nrow(rq)) return(paste0(sub("^[.]/", "", rq$file[1]), ":", rq$line[1]))
    h <- hit
    if (!any(h)) return("")
  }
  j <- which(h)[1]
  paste0(sub("^[.]/", "", src$path[j]), ":", src$line[j])
}, character(1))

# --- BUNDLE_DESIGN_KEYS, held to its own claim -------------------------------
# R/10_bundle.R says of that vector: "DERIVED, NOT CURATED. This list is
# exactly the set of `cfg$<key>` references in R/, minus `paths` and `pairing`
# which are handled separately." That is a checkable assertion and it belongs
# here, because it is the same failure class the audit exists for: a
# declaration and its consumer drifting apart.
#
# A key frozen as design but referenced by no `cfg$` in R/ is INERT DESIGN --
# `verify_bundle()` compares it between bundle and local config, so an edit
# to it is reported, but nothing anywhere acts on its value. That is not
# harmless: it reads as a control.
suppressWarnings(suppressMessages(source(file.path(root, "R", "10_bundle.R"))))
refs_lib <- cfg_refs(src, "lib")
refs_all <- cfg_refs(src, c("lib", "run", "targets"))
bdk_extra   <- setdiff(BUNDLE_DESIGN_KEYS, c(refs_all, "pairing"))
bdk_missing <- setdiff(refs_lib, c(BUNDLE_DESIGN_KEYS, "paths", "pairing"))

# --- the runner configs ------------------------------------------------------
# Same question, smaller surface. These are per-runner and are consumed by
# `run/internal.R`, `run/test_look.R` and `run/external.R` only, so a dead key
# here is a dead switch on a runner rather than on the design.

# FULL DEPTH, not one level. The runner configs are small and every key in them
# is a setting rather than a map keyed by design vocabulary, so there is no
# reason to stop at depth 2 -- and stopping there had a cost: `test_look.
# severity.apache_counterpart` and `.sofa_counterpart` in `config/internal.yml`
# sat at depth 3, were read by nothing, and this audit did not see them. Their
# twins in `config/external.yml` sat at depth 2 and were reported (finding F7),
# so the same defect was half-found. A depth limit that a declaration can hide
# under is the same failure the audit exists to catch, one level up.
.enumerate_deep <- function(v, prefix = "", vocab = character(0)) {
  out <- list()
  for (k in names(v)) {
    path <- if (nzchar(prefix)) paste0(prefix, ".", k) else k
    out[[length(out) + 1L]] <- list(key = k, path = path)
    x <- v[[k]]
    if (is.list(x) && !is.null(names(x)) && length(x) && !any(names(x) %in% vocab)) {
      out <- c(out, .enumerate_deep(x, path, vocab))
    }
  }
  out
}

# WHICH LAYER IS ALLOWED TO BE A RUNNER CONFIG'S ONLY CONSUMER.
#
# `reach_verdict()` calls a key read only from `tests/` VESTIGIAL, and for
# `config/config.yml` and the three runner configs that is exactly right: those
# configure the pipeline, so a design key nothing outside a test reads is dead
# weight. `config/attribution_eval.yml` is the exception and it is a structural
# one rather than a excuse. The attribution sub-pipeline is a GENERATOR and a
# CONSUMER that both live in `tests/` by the plan's own design
# (`docs/attribution_analysis_plan_20260906.md` section 10.2), because nothing
# in it fits a bundle or produces a reported model. Judging it against `run/`
# would report all 16 of its live keys as vestigial, which is a matcher that
# does not understand the idiom -- the same failure `cfg_req_paths()` was
# written to fix on 2026-09-03.
#
# So the expectation is DECLARED PER FILE rather than assumed, and the four
# files that must be read from a live pipeline layer still are.
RUNNER_CFG <- list(
  "internal.yml"        = c("lib", "run", "targets"),
  "external.yml"        = c("lib", "run", "targets"),
  "attribution.yml"     = c("lib", "run", "targets"),
  # Added 2026-09-06 with the file itself. Its keys declare the replicate
  # counts, the delta policy and the SELECTION RULE, and a dead key in a
  # selection rule is worse than a dead switch on a runner: the rule would still
  # be quoted in a paper as the criterion while nothing read it.
  "attribution_eval.yml" = c("lib", "run", "targets", "tests"))

runner_res <- do.call(rbind, lapply(names(RUNNER_CFG), function(f) {
  live <- RUNNER_CFG[[f]]
  y <- yaml::read_yaml(file.path(root, "config", f))
  ks <- .enumerate_deep(y, "", vocab)
  do.call(rbind, lapply(ks, function(e) {
    a <- grep_arms(src, e$key)
    hit <- a$dollar | a$bracket
    lay <- unique(c(src$layer[hit], req_layers(e$path, e$key)))
    data.frame(file = f, path = e$path, key = e$key,
               verdict = if (!length(lay)) "dead"
                         else if (any(lay %in% live)) "read" else "vestigial",
               n_accessor = sum(hit),
               layers = paste(sort(lay), collapse = "/"),
               stringsAsFactors = FALSE)
  }))
}))

# --- B5: keys the YAML parser turned into booleans ---------------------------
#
# YAML 1.1, which `yaml::read_yaml()` implements, treats the BARE WORDS `y`,
# `n`, `yes`, `no`, `on`, `off`, `true` and `false` as booleans -- in KEY
# position as well as in value position. So `n: 4000` does not declare a
# setting called `n`. It declares one called `FALSE`, and every read of
# `cfg$eval$n` returns NULL.
#
# THIS IS THE AUDIT'S OWN DEFECT CLASS IN A NEW COAT, and it is worse than the
# ones already found: the key is present in the file, spelled the way the
# reader expects, and visible in every diff. `cfg_req()` cannot save you either
# -- it raises on an ABSENT key, and from its point of view the key really is
# absent. FOUND 2026-09-04 in `config/attribution.yml`, where `eval.n: 4000`
# meant the arm silently evaluated on all 10,313 test stays rather than the
# 4,000 declared (finding F18).
#
# The check is a name scan and has no false positives: once mangled, the list
# element is literally named "TRUE" or "FALSE".
.bool_keys <- function(v, prefix = "") {
  out <- character(0)
  for (k in names(v)) {
    path <- if (nzchar(prefix)) paste0(prefix, ".", k) else k
    if (k %in% c("TRUE", "FALSE")) out <- c(out, path)
    if (is.list(v[[k]]) && !is.null(names(v[[k]])))
      out <- c(out, .bool_keys(v[[k]], path))
  }
  out
}

b5 <- do.call(rbind, lapply(list.files(file.path(root, "config"),
                                       pattern = "[.]yml$", full.names = TRUE),
                            function(f) {
  hits <- .bool_keys(yaml::read_yaml(f))
  if (!length(hits)) return(NULL)
  data.frame(file = basename(f), path = hits, stringsAsFactors = FALSE)
}))
if (is.null(b5)) b5 <- data.frame(file = character(0), path = character(0))

# --- report ------------------------------------------------------------------

md <- c(
  "# Audit B — config reachability",
  "",
  "Generated by `tests/audit_b_config.R`. Static analysis over `config/*.yml`",
  "and every `.R` file in `R/`, `run/`, `tests/` and `_targets.R`, with comments",
  "stripped. No data is read.",
  "",
  "Three nested matchers, weakest evidence last:",
  "",
  "- `verdict` — ACCESSOR only: `$key`, `[[\"key\"]]`, or a quoted argument to",
  "  `cfg_req()` / `cfg_flag()`. The three forms that can only ever be a config",
  "  lookup. This is the column to read.",
  "- `strict` — accessor, plus a bare quoted `\"key\"` anywhere.",
  "- `permissive` — the above, plus a bare symbol `key`.",
  "",
  "Rows where the three disagree are listed separately and are the ones a human",
  "must look at. The accessor matcher can miss a key reached through a variable",
  "(`cfg[[nm]]`); the looser two can be fooled by an unrelated string or local of",
  "the same name — `paths.out` is the worked example, scoring reachable on a",
  "`out_dir <- \"out\"` in a test and in fact read by nothing.",
  "")

tally <- as.data.frame(table(res$verdict), stringsAsFactors = FALSE)
names(tally) <- c("verdict", "n")
md <- c(md, "## Tally, config/config.yml", "", md_table(tally), "")

dead <- res[res$verdict == "dead", c("path", "strict", "permissive", "n_strict"), drop = FALSE]
md <- c(md, "## DEAD — declared, reached by no accessor anywhere", "")
md <- c(md, if (nrow(dead)) md_table(dead) else "_none_", "")

vest <- res[res$verdict == "vestigial", c("path", "n_accessor", "first_live"), drop = FALSE]
md <- c(md, "## VESTIGIAL — accessor hits only from tests/", "")
md <- c(md, if (nrow(vest)) md_table(vest) else "_none_", "")

dis <- res[res$disagree, c("path", "verdict", "strict", "permissive",
                           "n_accessor", "n_strict"), drop = FALSE]
md <- c(md, "## Matchers disagree — adjudicate these", "")
md <- c(md, if (nrow(dis)) md_table(dis) else "_none_", "")

md <- c(md, "## BUNDLE_DESIGN_KEYS, against its own derivation claim", "",
        "`R/10_bundle.R` documents that vector as \"exactly the set of `cfg$<key>`",
        "references in R/, minus `paths` and `pairing`\". Two ways it can be wrong:",
        "",
        "- **frozen but inert** — in `BUNDLE_DESIGN_KEYS` and referenced by no",
        "  `cfg$` anywhere in `R/`, `run/` or `_targets.R`. `verify_bundle()` will",
        "  compare it, so an edit is reported; but nothing acts on its value, and it",
        "  reads as a control that is in force.",
        "- **read but unfrozen** — a `cfg$<key>` reference in `R/` that is NOT",
        "  frozen, so an apply site could traverse a different design invisibly.",
        "  This is the direction that breaks hard rule 8.",
        "")
md <- c(md, "**frozen but inert:** ",
        if (length(bdk_extra)) paste0("`", paste(bdk_extra, collapse = "`, `"), "`") else "_none_", "")
md <- c(md, "", "**read but unfrozen:** ",
        if (length(bdk_missing)) paste0("`", paste(bdk_missing, collapse = "`, `"), "`") else "_none_", "")

md <- c(md, "", "## Full table", "",
        md_table(res[order(res$verdict, res$path),
                     c("path", "depth", "verdict", "n_accessor", "n_files",
                       "layers", "first_live")]), "")

md <- c(md, "## B5 - keys the YAML parser turned into booleans", "",
        "YAML 1.1 reads the bare words `y`, `n`, `yes`, `no`, `on`, `off`, `true`",
        "and `false` as booleans in KEY position too, so `n: 4000` declares a key",
        "named `FALSE` and every read of it returns NULL. The key is present,",
        "correctly spelled and visible in the diff, and `cfg_req()` cannot help:",
        "from its point of view the key is absent. Found once, in",
        "`config/attribution.yml` (finding F18).", "")
md <- c(md, if (nrow(b5)) md_table(b5) else "_none_", "")

md <- c(md, "## Runner configs", "",
        "`config/internal.yml`, `config/external.yml`, `config/attribution.yml`",
        "and `config/attribution_eval.yml`. Same matcher, strict only.",
        "")
rt <- as.data.frame(table(runner_res$file, runner_res$verdict), stringsAsFactors = FALSE)
names(rt) <- c("file", "verdict", "n")
md <- c(md, md_table(rt[rt$n > 0, ]), "")
rdead <- runner_res[runner_res$verdict != "read", c("file", "path", "verdict", "n_accessor")]
md <- c(md, "### Not read from a live path", "",
        if (nrow(rdead)) md_table(rdead) else "_none_", "")

audit_write(md, "audit_b_config.md", root)

cat("\n=== AUDIT B: config reachability ===\n")
print(tally, row.names = FALSE)
cat("\nDEAD keys (", nrow(dead), "):\n", sep = "")
if (nrow(dead)) print(dead, row.names = FALSE) else cat("  none\n")
cat("\nVESTIGIAL keys (", nrow(vest), "):\n", sep = "")
if (nrow(vest)) print(vest, row.names = FALSE) else cat("  none\n")
cat("\nStrict/permissive disagreements (", nrow(dis), "):\n", sep = "")
if (nrow(dis)) print(dis, row.names = FALSE) else cat("  none\n")
cat("\nRunner-config keys not read from a live path (", nrow(rdead), "):\n", sep = "")
if (nrow(rdead)) print(rdead, row.names = FALSE) else cat("  none\n")
writeLines("")
cat("B5 keys mangled into booleans by YAML 1.1 (", nrow(b5), "):", sep = "")
writeLines("")
if (nrow(b5)) print(b5, row.names = FALSE) else writeLines("  none")
