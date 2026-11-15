# Shared fixtures + accessors for the CNA test suite.
#
# Engine functions are internal (the NAMESPACE exports nothing user-facing), so
# they are reached with getFromNamespace(). The reference engine (Nestimate /
# cograph) is attached so equivalence tests can call it side by side.

suppressMessages({
  loadNamespace("CNA")
  library(Nestimate)
  library(cograph)
})

cna <- function(name) getFromNamespace(name, "CNA")

# --- repository paths, resolved regardless of cwd ------------------------------
repo_root <- function() {
  root <- getOption("CNA_TEST_ROOT", default = NULL)
  cands <- c(root, ".", "../..")
  hit <- Filter(function(p) !is.null(p) && file.exists(file.path(p, "tools", "generate_jamovi.R")),
                cands)
  if (!length(hit)) testthat::skip("repository root not locatable")
  normalizePath(hit[[1]])
}
read_repo <- function(...) readLines(file.path(repo_root(), ...), warn = FALSE, encoding = "UTF-8")

ANALYSES <- c("CNAEventSingle", "CNAEventCluster", "CNABinarySingle",
              "CNABinaryCluster", "CNADelimited")

regulation_data <- function(n_actors = NULL, seed = 1) {
  d <- utils::read.csv(file.path(repo_root(), "data", "Regulation_long.csv"),
                       stringsAsFactors = FALSE)
  if (!is.null(n_actors)) {
    set.seed(seed)
    d <- d[d$Actor %in% sample(unique(d$Actor), n_actors), , drop = FALSE]
  }
  d
}

# --- deterministic synthetic fixtures -----------------------------------------
# Long-format events: several actions per actor, clustered in time so they
# co-occur within the threshold; two groups.
make_event_data <- function(n_actor = 8, actions = LETTERS[1:5], per = 10, seed = 42) {
  set.seed(seed)
  actors <- sprintf("s%02d", seq_len(n_actor))
  half <- actors[seq_len(ceiling(n_actor / 2))]
  base_time <- as.POSIXct("2021-01-01 00:00:00", tz = "UTC")
  do.call(rbind, lapply(actors, function(a) {
    data.frame(Actor = a, Action = sample(actions, per, replace = TRUE),
               Time = base_time + sort(sample.int(3000, per)),
               Group = if (a %in% half) "G1" else "G2", stringsAsFactors = FALSE)
  }))
}

# One-hot data: 0/1 code columns + actor / session / group.
make_onehot_data <- function(n = 48, codes = paste0("code", 1:4), seed = 7) {
  set.seed(seed)
  m <- as.data.frame(matrix(stats::rbinom(n * length(codes), 1, 0.5), ncol = length(codes)))
  names(m) <- codes
  m$Actor   <- rep(sprintf("a%02d", seq_len(max(1, n %/% 4))), length.out = n)
  m$Session <- rep(1:2, length.out = n)
  m$Group   <- rep(c("G1", "G2"), length.out = n)
  m
}

# Delimited data: keyword lists separated by ";".
make_delimited_data <- function(n = 40, items = c("network", "graph", "matrix",
                                                   "algebra", "model"), seed = 3) {
  set.seed(seed)
  data.frame(id = seq_len(n),
             keywords = vapply(seq_len(n), function(i)
               paste(sample(items, sample(2:4, 1)), collapse = "; "), character(1)),
             stringsAsFactors = FALSE)
}

# --- engine option lists (plain lists, as the engine functions consume) ------
opts_event <- function(...) modifyList(list(
  buildModel_variables_long_action = "Action", buildModel_variables_actor = "Actor",
  buildModel_variables_long_time = "Time", buildModel_variables_long_order = NULL,
  buildModel_variables_group = NULL, buildModel_scaling = "noScaling",
  buildModel_threshold = 3000), list(...))

opts_onehot <- function(...) modifyList(list(
  buildModel_variables_onehot = paste0("code", 1:4), buildModel_variables_actor = "Actor",
  buildModel_variables_session = "Session", buildModel_variables_group = NULL,
  buildModel_window = 1L, buildModel_window_type = "tumbling",
  buildModel_scaling = "noScaling"), list(...))

opts_delimited <- function(...) modifyList(list(
  buildModel_variables_field = "keywords", buildModel_variables_by = NULL,
  buildModel_sep = ";", buildModel_similarity = "none", buildModel_min_occur = 1L,
  buildModel_edge_threshold = 0, buildModel_diagonal = TRUE, buildModel_top_n = 0L),
  list(...))

fixture_event_net <- function(...) cna("cna_build_event_nest")(make_event_data(...), opts_event())
fixture_group_net <- function(...)
  cna("cna_build_event_nest")(make_event_data(...), opts_event(buildModel_variables_group = "Group"))

# --- jmvcore analysis runners ---------------------------------------------------
# Every Variable / Level option has no default in the generated runner, so each
# must be passed (NULL when unused). merge_opts keeps explicit NULL overrides,
# which modifyList() would drop.
merge_opts <- function(defaults, overrides) {
  for (nm in names(overrides)) defaults[nm] <- list(overrides[[nm]])
  defaults
}

runner_defaults <- list(
  CNAEventSingle = list(buildModel_variables_long_action = "Action",
                        buildModel_variables_actor = "Actor",
                        buildModel_variables_long_time = "Time",
                        buildModel_variables_long_order = NULL,
                        buildModel_variables_group = NULL, buildModel_threshold = 3000,
                        compare_group_i = NULL, compare_group_j = NULL),
  CNAEventCluster = list(buildModel_variables_long_action = "Action",
                         buildModel_variables_actor = "Actor",
                         buildModel_variables_long_time = "Time",
                         buildModel_variables_long_order = NULL, buildModel_threshold = 3000),
  CNABinarySingle = list(buildModel_variables_onehot = paste0("code", 1:4),
                         buildModel_variables_actor = "Actor",
                         buildModel_variables_session = "Session",
                         buildModel_variables_group = NULL,
                         compare_group_i = NULL, compare_group_j = NULL),
  CNABinaryCluster = list(buildModel_variables_onehot = paste0("code", 1:4),
                          buildModel_variables_actor = "Actor",
                          buildModel_variables_session = "Session"),
  CNADelimited = list(buildModel_variables_field = "keywords",
                      buildModel_variables_by = NULL))

run_cna <- function(analysis, data, ...) {
  runner <- cna(analysis)
  do.call(runner, merge_opts(c(list(data = data), runner_defaults[[analysis]]), list(...)))
}
run_event_single <- function(data, ...) run_cna("CNAEventSingle", data, ...)

# Every display option of an analysis switched on, with small iteration counts.
all_sections_on <- function(analysis) {
  f <- names(formals(cna(analysis)))
  on <- grep("_show_|^seed_set$|^clustering_run$", f, value = TRUE)
  small <- list(bootstrap_iteration = 30L, permutation_iter = 30L, stability_iter = 50L,
                cluster_choice_kmax = 3L)
  c(stats::setNames(rep(list(TRUE), length(on)), on),
    small[intersect(names(small), f)])
}

# A results table as a data.frame (0 rows if empty or absent).
table_df <- function(results, name) {
  df <- tryCatch(results[[name]]$asDF, error = function(e) NULL)
  if (is.null(df)) data.frame() else df
}

# Section headings that report a failure ("This analysis could not be completed").
section_failures <- function(results) {
  items <- results$itemNames
  texts <- vapply(items, function(nm) {
    it <- results[[nm]]
    if (inherits(it, "Preformatted")) it$content %||% "" else ""
  }, character(1))
  texts[grepl("could not be completed", texts, fixed = TRUE)]
}

# Render a visible image; returns "ok", the error raised, or the failure that
# cna_render() drew into the image instead of raising ("This plot could not
# be drawn" + the message). The drawn text is read back from SVG, so a plot
# whose renderer failed no longer passes as "ok".
render_image <- function(image) {
  testthat::skip_if_not_installed("svglite")
  f <- tempfile(fileext = ".svg")
  svglite::svglite(f, width = image$width / 72, height = image$height / 72)
  status <- tryCatch({ image$.render(); "ok" }, error = function(e) conditionMessage(e),
                     finally = grDevices::dev.off())
  lines <- readLines(f, warn = FALSE)
  unlink(f)
  texts <- sub("<text[^>]*>([^<]*)</text>", "\\1",
               unlist(regmatches(lines, gregexpr("<text[^>]*>[^<]*</text>", lines))))
  if (identical(status, "ok") && any(grepl("could not be drawn", texts, fixed = TRUE)))
    status <- paste(texts, collapse = " ")
  status
}

# Text drawn into an SVG by an expression (cograph labels are <text> nodes).
svg_texts <- function(expr) {
  testthat::skip_if_not_installed("svglite")
  f <- tempfile(fileext = ".svg")
  svglite::svglite(f, 7, 6)
  tryCatch(force(expr), finally = grDevices::dev.off())
  lines <- readLines(f, warn = FALSE)
  unlink(f)
  texts <- unlist(regmatches(lines, gregexpr("<text[^>]*>[^<]*</text>", lines)))
  sub("<text[^>]*>([^<]*)</text>", "\\1", texts)
}

# Parsed jamovi definitions of an analysis ("a" options, "r" results, "u" UI).
read_def <- function(analysis, part = c("a", "r", "u")) {
  part <- match.arg(part)
  yaml::read_yaml(file.path(repo_root(), "jamovi", paste0(analysis, ".", part, ".yaml")))
}
option_names <- function(analysis) vapply(read_def(analysis, "a")$options, `[[`, "", "name")
option_def <- function(analysis, name)
  Filter(function(o) identical(o$name, name), read_def(analysis, "a")$options)[[1]]
