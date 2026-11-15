#!/usr/bin/env Rscript
# tools/generate_jamovi.R
#
# Single source of truth for the CNA module's jamovi definitions. For every
# analysis it writes
#   jamovi/<A>.a.yaml   options
#   jamovi/<A>.u.yaml   user interface
#   jamovi/<A>.r.yaml   results
#   R/<A>.b.R           thin backend (all work is in the shared R/cna_*.R engine)
# and the module manifest jamovi/0000.yaml.
#
# Why a generator: the four (now five) analyses share almost every section, and
# hand-maintaining near-identical copies is what produced the module's earlier
# bugs (an option missing from one analysis's clearWith, controls the engine
# never read). Here each result declares what it depends on and its clearWith is
# derived from that:
#   * model options               -> every result
#   * an analysis's own options   -> that analysis's state holder and outputs
#   * display options             -> only the table or plot they format
# Expensive computations (bootstrap, permutation, stability, communities, ...)
# are cached in hidden state holders whose clearWith excludes display options,
# so changing how a plot looks never reruns the resampling behind it.
#
# Usage, from the repository root (or pass another root to write into):
#   Rscript tools/generate_jamovi.R [root]
#   <jamovi Rscript> -e 'jmvtools::prepare(".")'   # regenerates R/<A>.h.R

VERSION <- "1.1.0"
DATE    <- "2026-11-05"

ROOT <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(ROOT)) ROOT <- "."
stopifnot("the root must hold the package DESCRIPTION" =
            file.exists(file.path(ROOT, "DESCRIPTION")))

# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------

ANALYSES <- list(
  list(name = "CNAEventSingle", kind = "event_single", data = "event",
       cluster = FALSE, title = "Event Data (Single/Group)",
       subgroup = "Event Data", menu = "Single / Group"),
  list(name = "CNAEventCluster", kind = "event_cluster", data = "event",
       cluster = TRUE, title = "Event Data (Clustering)",
       subgroup = "Event Data", menu = "Clustering"),
  list(name = "CNABinarySingle", kind = "binary_single", data = "onehot",
       cluster = FALSE, title = "Binary Data (Single/Group)",
       subgroup = "Binary Data", menu = "Single / Group"),
  list(name = "CNABinaryCluster", kind = "binary_cluster", data = "onehot",
       cluster = TRUE, title = "Binary Data (Clustering)",
       subgroup = "Binary Data", menu = "Clustering"),
  list(name = "CNADelimited", kind = "delimited", data = "delimited",
       cluster = FALSE, title = "Delimited Data (Keywords / Items)",
       subgroup = "Delimited Data", menu = "Keywords / Items")
)

# Capabilities derived from the data type and mode.
describe <- function(a) {
  a$group    <- !a$cluster && a$data != "delimited"  # optional Group variable
  a$multi    <- a$group || a$cluster                 # can hold several networks
  a$resample <- a$data != "delimited"                # has source data to resample
  # transition + co-occurrence view (events and one-hot; not delimited data,
  # which has no order). The module has no sequence analyses: co-occurrence
  # ignores order, so sequence plots and subsequence comparison do not belong.
  a$mixed    <- a$data != "delimited"
  a
}
ANALYSES <- lapply(ANALYSES, describe)

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

compact <- function(x) Filter(Negate(is.null), x)

# A control is a list with a `type`; anything else is a list of controls.
splice <- function(x) {
  x <- compact(x)
  do.call(c, lapply(x, function(el) if (!is.null(el$type)) list(el) else compact(el)))
}

# Mathematical sans-serif bold, as in the existing "━━ Network Analysis ━━"
# dividers.
bold_sans <- function(text) {
  cp <- utf8ToInt(text)
  upper <- cp >= 65L & cp <= 90L
  lower <- cp >= 97L & cp <= 122L
  cp[upper] <- cp[upper] - 65L + 0x1D5D4L
  cp[lower] <- cp[lower] - 97L + 0x1D5EEL
  intToUtf8(cp)
}

# ---------------------------------------------------------------------------
# Option constructors (.a.yaml)
# ---------------------------------------------------------------------------

o_bool <- function(name, title, default = FALSE)
  list(name = name, title = title, type = "Bool", default = default)

o_int <- function(name, title, default, min = NULL, max = NULL)
  compact(list(name = name, title = title, type = "Integer",
               default = as.integer(default),
               min = if (!is.null(min)) as.integer(min),
               max = if (!is.null(max)) as.integer(max)))

o_num <- function(name, title, default, min = NULL, max = NULL)
  compact(list(name = name, title = title, type = "Number",
               default = as.numeric(default),
               min = if (!is.null(min)) as.numeric(min),
               max = if (!is.null(max)) as.numeric(max)))

o_list <- function(name, title, choices, default = names(choices)[[1]])
  list(name = name, title = title, type = "List",
       options = unname(Map(function(n, t) list(name = n, title = t),
                            names(choices), choices)),
       default = default)

o_string <- function(name, title, default)
  list(name = name, title = title, type = "String", default = default)

o_var   <- function(name, title) list(name = name, title = title, type = "Variable")
o_vars  <- function(name, title) list(name = name, title = title, type = "Variables")
o_level <- function(name, title, variable)
  list(name = name, title = title, type = "Level", variable = variable)

# ---------------------------------------------------------------------------
# UI constructors (.u.yaml)
# ---------------------------------------------------------------------------

# Child controls (`...`) are drawn indented under a check box or radio button,
# the way jamovi's own analyses attach a setting to the choice it belongs to;
# they are enabled only while their parent is ticked.
with_children <- function(control, children) {
  children <- splice(children)
  if (!length(children)) return(control)
  control$children <- lapply(children, function(ch) {
    ch$enable <- ch$enable %||% sprintf("(%s)", control$name)
    ch
  })
  control
}

u_check <- function(name, label = NULL, ..., enable = NULL)
  with_children(compact(list(type = "CheckBox", name = name, label = label,
                             enable = enable)), list(...))

u_text <- function(name, label = NULL, format = "number", enable = NULL, suffix = NULL)
  compact(list(type = "TextBox", name = name, format = format,
               label = label, suffix = suffix, enable = enable))

u_combo <- function(name, enable = NULL, label = NULL)
  compact(list(type = "ComboBox", name = name, label = label, enable = enable))

u_radio <- function(option, part, ...)
  with_children(list(type = "RadioButton", name = paste0(option, "_", part),
                     optionName = option, optionPart = part), list(...))

u_radios <- function(option, parts) lapply(parts, function(p) u_radio(option, p))

u_level <- function(name) list(type = "LevelSelector", name = name)

u_label <- function(label, ...)
  list(type = "Label", label = label, children = splice(list(...)))

u_col <- function(...)
  list(type = "LayoutBox", stretchFactor = 1L, children = splice(list(...)))

u_collapse <- function(label, ...)
  list(type = "CollapseBox", label = label, margin = "large",
       stretchFactor = 1L, style = "inline", collapsed = TRUE,
       children = splice(list(...)))

# The label is emitted as an ASCII placeholder and swapped for the bold text
# after serialisation (see emit()), so the YAML writer never escapes or folds it.
u_divider <- function(text)
  list(type = "Label", label = paste0("DIVIDER<", text, ">"),
       style = "list-inline", margin = "large")

u_target <- function(label, name, single = TRUE)
  list(type = "TargetLayoutBox", label = label,
       children = list(compact(list(type = "VariablesListBox", name = name,
                                    maxItemCount = if (single) 1L,
                                    isTarget = TRUE))))

# ---------------------------------------------------------------------------
# Result constructors (.r.yaml)
# ---------------------------------------------------------------------------

# deps: option names, or tokens expanded per analysis:
#   "@model" every option that changes the network(s)
#   "@look"  network-plot appearance options
#   "@seed"  the reproducibility seed
r_item <- function(name, type, title = NULL, visible = FALSE, deps = character(),
                   render = FALSE, width = NULL, height = NULL,
                   columns = NULL, content = NULL)
  list(name = name, type = type, title = title, visible = visible,
       deps = deps, render = render, width = width, height = height,
       columns = columns, content = content)

r_col <- function(name, title, type = "text", format = NULL, visible = NULL)
  compact(list(name = name, title = title, type = type, format = format,
               visible = visible))

r_title <- function(name, title, visible)
  r_item(name, "Preformatted", title = title, visible = visible)

r_state <- function(name, deps)
  r_item(name, "Preformatted", visible = FALSE, deps = deps)

r_image <- function(name, visible, deps, title = NULL, width = 600L, height = 450L)
  r_item(name, "Image", title = title, visible = visible, deps = deps,
         render = TRUE, width = width, height = height)

r_table <- function(name, title, visible, deps, columns)
  r_item(name, "Table", title = title, visible = visible, deps = deps,
         columns = columns)

# "(a || b)" visibility expressions, evaluated by jmvcore against the options.
vis <- function(...) paste0("(", paste(c(...), collapse = " || "), ")")

# The group column: shown only when there is more than one network.
group_col <- function(a) {
  if (a$cluster) return(r_col("group", "Cluster"))
  if (a$group) return(r_col("group", "Group", visible = "(buildModel_variables_group)"))
  r_col("group", "Group", visible = FALSE)
}

# ---------------------------------------------------------------------------
# Shared vocabularies
# ---------------------------------------------------------------------------

# UI option -> cograph::centrality() measure. Co-occurrence networks are
# undirected, so one Strength (and one Closeness) replaces the In/Out pair.
# Listed by family (degree, path, spectral, local structure), which is the
# order of the check boxes and of the table columns.
MEASURES <- list(
  Strength       = list(title = "Strength",                 default = TRUE),
  Degree         = list(title = "Degree",                   default = FALSE),
  Coreness       = list(title = "Coreness (k-core)",        default = FALSE),
  Betweenness    = list(title = "Betweenness",              default = TRUE),
  BetweennessRSP = list(title = "RSP betweenness",          default = FALSE),
  Closeness      = list(title = "Closeness",                default = FALSE),
  Harmonic       = list(title = "Harmonic closeness",       default = FALSE),
  Eccentricity   = list(title = "Eccentricity",             default = FALSE),
  Eigenvector    = list(title = "Eigenvector",              default = FALSE),
  PageRank       = list(title = "PageRank",                 default = FALSE),
  Diffusion      = list(title = "Diffusion",                default = FALSE),
  Clustering     = list(title = "Clustering coefficient",   default = FALSE),
  Constraint     = list(title = "Constraint (structural holes)", default = FALSE),
  EffectiveSize  = list(title = "Effective size",           default = FALSE)
)
MEASURE_OPTS <- paste0("centrality_", names(MEASURES))

P_ADJUST <- c(none = "None", holm = "Holm", hochberg = "Hochberg",
              hommel = "Hommel", bonferroni = "Bonferroni",
              BH = "Benjamini-Hochberg", BY = "Benjamini-Yekutieli",
              fdr = "FDR")

SCALING <- c(noScaling = "No scaling", minmax = "MinMax", max = "Max",
             rank = "Rank")

LOOK_OPTS <- c("plot_layout", "plot_edge_labels", "plot_minimum", "plot_cut",
               "plot_edge_label_size", "plot_node_size", "plot_label_size")

SEED_OPTS <- c("seed_set", "seed_value")

# ---------------------------------------------------------------------------
# Sections. Each returns list(opts, ui, res); ui is a list of UI nodes for the
# top level (sections that belong elsewhere return NULL ui).
# ---------------------------------------------------------------------------

data_vars <- function(a) switch(a$data,
  event = c("buildModel_variables_long_action", "buildModel_variables_actor",
            "buildModel_variables_long_time", "buildModel_variables_long_order",
            if (a$group) "buildModel_variables_group"),
  onehot = c("buildModel_variables_onehot", "buildModel_variables_actor",
             "buildModel_variables_session",
             if (a$group) "buildModel_variables_group"),
  delimited = c("buildModel_variables_field", "buildModel_variables_by"))

model_opts <- function(a) c(
  data_vars(a),
  if (a$data != "delimited") "buildModel_scaling",
  if (a$data == "event") "buildModel_threshold",
  if (a$data == "onehot") c("buildModel_window", "buildModel_window_type"),
  if (a$data == "delimited") c("buildModel_sep", "buildModel_similarity",
                               "buildModel_min_occur", "buildModel_edge_threshold",
                               "buildModel_diagonal", "buildModel_top_n"),
  if (a$cluster) c("clustering_run", "clustering_k", "clustering_dissimilarity",
                   "clustering_method", SEED_OPTS))

sec_inputs <- function(a) {
  actor_title <- if (a$cluster && a$data == "event") "Actor" else "Actor (Optional)"
  opts <- switch(a$data,
    event = list(
      o_var("buildModel_variables_long_action", "Action"),
      o_var("buildModel_variables_actor", actor_title),
      o_var("buildModel_variables_long_time", "Time (Optional)"),
      o_var("buildModel_variables_long_order", "Order (Optional)")),
    onehot = list(
      o_vars("buildModel_variables_onehot", "One-Hot Columns"),
      o_var("buildModel_variables_actor", actor_title),
      o_var("buildModel_variables_session", "Session (Optional)")),
    delimited = list(
      o_vars("buildModel_variables_field", "Item Columns"),
      o_var("buildModel_variables_by", "Transaction ID (Optional)")))
  if (a$group)
    opts <- c(opts, list(o_var("buildModel_variables_group", "Group Variable (Optional)")))

  targets <- lapply(opts, function(o)
    u_target(o$title, o$name, single = identical(o$type, "Variable")))

  list(opts = c(list(list(name = "data", type = "Data")), opts),
       ui = list(list(type = "VariableSupplier", persistentItems = FALSE,
                      stretchFactor = 1L, children = targets)),
       res = list(
         r_item("instructions", "Html", title = "Instructions", visible = TRUE,
                deps = data_vars(a)),
         r_item("errorText", "Preformatted", title = "Error", visible = FALSE,
                deps = "@model"),
         if (a$cluster) r_item("clusteringNote", "Html", visible = FALSE,
                               deps = "clustering_run"),
         r_title("analysisTitle", sprintf("%s — CNA %s", a$title, VERSION),
                 visible = TRUE)))
}

sec_model <- function(a) {
  cluster_opts <- if (a$cluster) list(
    o_bool("clustering_run", "Run clustering", FALSE),
    o_int("clustering_k", "Number of Clusters (k)", 2, 2, 20),
    o_list("clustering_dissimilarity", "Dissimilarity",
           c(hamming = "Hamming", osa = "OSA", lv = "Levenshtein",
             dl = "Damerau-Levenshtein", lcs = "LCS", qgram = "q-gram",
             cosine = "Cosine", jaccard = "Jaccard", jw = "Jaro-Winkler")),
    o_list("clustering_method", "Method",
           c(pam = "PAM", ward.D = "Ward D", ward.D2 = "Ward D2",
             complete = "Complete", average = "Average", single = "Single",
             mcquitty = "McQuitty", median = "Median", centroid = "Centroid")))

  build_opts <- c(
    if (a$data != "delimited") list(o_list("buildModel_scaling", "Scaling", SCALING)),
    if (a$data == "event") list(o_int("buildModel_threshold", "Time threshold (s)", 900)),
    if (a$data == "onehot") list(
      o_int("buildModel_window", "Window size", 1, 1, 100),
      o_list("buildModel_window_type", "Window type",
             c(tumbling = "Tumbling", sliding = "Sliding"))),
    if (a$data == "delimited") list(
      o_string("buildModel_sep", "Separator", ";"),
      o_int("buildModel_min_occur", "Minimum item occurrences", 1, 1),
      o_bool("buildModel_diagonal", "Keep diagonal (item frequency)", TRUE),
      # No "relative" (row-normalised): it is asymmetric, but Nestimate 1.0.7
      # returns it as an undirected network that keeps only one direction.
      o_list("buildModel_similarity", "Similarity",
             c(none = "None (raw counts)", jaccard = "Jaccard",
               cosine = "Cosine (Salton)", inclusion = "Inclusion (Simpson)",
               association = "Association strength", dice = "Dice",
               equivalence = "Equivalence")),
      o_num("buildModel_edge_threshold", "Minimum edge weight", 0, 0),
      o_int("buildModel_top_n", "Top N edges (0 = all)", 0, 0)))

  basic_opts <- list(
    o_bool("buildModel_show_matrix", if (a$data == "delimited") "Matrix" else "Co-occurrence matrix", FALSE),
    o_bool("buildModel_show_summary", "Network summary", FALSE),
    o_bool("seed_set", "Set seed", FALSE),
    o_int("seed_value", "Seed", 1234, 0))

  # Left: what defines the network. Right: the three main outputs, then the
  # seed. Every other output lives in a collapsible section below.
  network_box <- u_label("Network",
    switch(a$data,
      event = list(u_text("buildModel_threshold", "Time threshold", suffix = "s")),
      onehot = list(u_text("buildModel_window", "Window size"),
                    u_combo("buildModel_window_type")),
      delimited = list(
        u_text("buildModel_sep", "Separator", format = "string"),
        u_text("buildModel_min_occur", "Minimum item occurrences"),
        u_combo("buildModel_similarity"),
        u_text("buildModel_edge_threshold", "Minimum edge weight"),
        u_text("buildModel_top_n", "Top edges (0 = all)"),
        u_check("buildModel_diagonal", "Keep diagonal (item frequency)"))),
    if (a$data != "delimited") u_combo("buildModel_scaling"))
  output_box <- u_label("Output",
                        u_check("buildModel_show_plot", "Network plot"),
                        u_check("buildModel_show_matrix",
                                if (a$data == "delimited") "Matrix" else "Co-occurrence matrix"),
                        u_check("buildModel_show_summary", "Network summary"))
  seed_box <- u_label("Reproducibility",
                      u_check("seed_set", "Set seed",
                              u_text("seed_value", "Seed")))
  cluster_box <- if (a$cluster)
    u_label("Clustering",
            u_check("clustering_run", "Run clustering"),
            u_text("clustering_k", "Number of clusters (k)"),
            u_combo("clustering_dissimilarity"), u_combo("clustering_method"))

  columns <- if (a$cluster) {
    list(u_col(cluster_box, network_box), u_col(output_box, seed_box))
  } else {
    list(u_col(network_box), u_col(output_box, seed_box))
  }

  # every model output, so a message written to the heading (a build note or a
  # mixed-network error) is visible whichever output asked for it
  show_any <- c("buildModel_show_matrix", "buildModel_show_summary",
                "buildModel_show_plot", "buildModel_show_histo",
                if (a$data != "delimited") c("buildModel_show_frequencies",
                                             "buildModel_show_mosaic"),
                if (a$mixed) "buildModel_show_mixed")
  list(opts = c(cluster_opts, build_opts, basic_opts),
       ui = list(list(type = "LayoutBox", stretchFactor = 1L, style = "inline",
                      children = columns)),
       res = list(
         r_title("buildModelTitle", "Network", visible = vis(show_any)),
         r_state("buildModelContent", deps = "@model"),
         r_table("buildModel_matrixTable",
                 if (a$data == "delimited") "Co-occurrence / Similarity Matrix"
                 else "Co-occurrence Matrix",
                 visible = vis("buildModel_show_matrix"), deps = "@model",
                 columns = list(group_col(a), r_col("state", "State"))),
         r_table("buildModel_summaryTable", "Network Summary",
                 visible = vis("buildModel_show_summary"), deps = "@model",
                 columns = list(group_col(a), r_col("metric", "Metric"),
                                r_col("value", "Value", "number")))))
}

sec_cluster_diag <- function(a) {
  if (!a$cluster) return(NULL)
  choice_deps <- c(setdiff(model_opts(a), c("clustering_run", "clustering_k")),
                   "cluster_choice_kmin", "cluster_choice_kmax")
  list(
    opts = list(
      o_bool("cluster_show_diagnostics", "Diagnostics table", FALSE),
      o_bool("cluster_show_diagnostics_plot", "Diagnostics plot", FALSE),
      o_bool("cluster_show_choice", "Choosing k table", FALSE),
      o_bool("cluster_show_choice_plot", "Choosing k plot", FALSE),
      o_int("cluster_choice_kmin", "Smallest k", 2, 2, 20),
      o_int("cluster_choice_kmax", "Largest k", 6, 2, 20)),
    ui = list(u_collapse("Cluster Diagnostics",
      u_col(u_label("Output",
                    u_check("cluster_show_diagnostics", "Diagnostics table"),
                    u_check("cluster_show_diagnostics_plot", "Diagnostics plot"),
                    u_check("cluster_show_choice", "Choosing k table"),
                    u_check("cluster_show_choice_plot", "Choosing k plot"))),
      u_col(u_label("Choosing k",
                    u_text("cluster_choice_kmin", "Smallest k"),
                    u_text("cluster_choice_kmax", "Largest k"))))),
    res = list(
      r_title("clusterTitle", "Clustering",
              vis("cluster_show_diagnostics", "cluster_show_diagnostics_plot",
                  "cluster_show_choice", "cluster_show_choice_plot")),
      r_table("clusterDiagTable", "Cluster Diagnostics",
              vis("cluster_show_diagnostics"), deps = "@model",
              columns = list(r_col("cluster", "Cluster"),
                             r_col("size", "Size", "integer"),
                             r_col("pct", "%", "number"),
                             r_col("mean_within_dist", "Mean within distance", "number"),
                             r_col("silhouette", "Mean silhouette", "number"))),
      r_image("clusterDiag_plot", vis("cluster_show_diagnostics_plot"),
              deps = "@model"),
      r_state("clusterChoiceState", deps = choice_deps),
      r_table("clusterChoiceTable", "Choosing the Number of Clusters",
              vis("cluster_show_choice"), deps = choice_deps,
              columns = list(r_col("k", "k", "integer"),
                             r_col("silhouette", "Mean silhouette", "number"),
                             r_col("mean_within_dist", "Mean within distance", "number"),
                             r_col("min_size", "Smallest cluster", "integer"),
                             r_col("max_size", "Largest cluster", "integer"),
                             r_col("size_ratio", "Size ratio", "number"))),
      r_image("clusterChoice_plot", vis("cluster_show_choice_plot"),
              deps = choice_deps, width = 700L)))
}

sec_visual <- function(a) {
  data_plots <- a$data != "delimited"
  opts <- c(
    list(o_bool("buildModel_show_plot", "Network plot", TRUE),
         o_bool("buildModel_show_histo", "Edge weight histogram", FALSE)),
    if (data_plots) list(
      o_bool("buildModel_show_frequencies", "Code frequencies", FALSE),
      o_bool("buildModel_show_mosaic", "Mosaic plot", FALSE)),
    if (a$mixed) list(
      o_bool("buildModel_show_mixed", "Mixed network (transition + co-occurrence)", FALSE),
      o_list("buildModel_mixed_type", "Mixed network display",
             c(group = "Side by side", overlay = "Overlay"))),
    list(
      o_list("plot_layout", "Layout",
             c(oval = "Oval", circle = "Circle", spring = "Spring",
               fr = "Fruchterman-Reingold", gephi_fr = "Gephi (FR)",
               star = "Star", grid = "Grid", random = "Random")),
      o_bool("plot_edge_labels", "Edge labels", TRUE),
      o_num("plot_minimum", "Hide edges below", 0, 0),
      o_num("plot_cut", "Cut value (0 = automatic)", 0, 0),
      o_num("plot_edge_label_size", "Edge label size", 0.8, 0.1, 5),
      o_num("plot_node_size", "Node size (0 = automatic)", 0, 0, 30),
      o_num("plot_label_size", "Node label size (0 = automatic)", 0, 0, 5)))

  # The network plot itself is switched on in the top area; this section holds
  # its appearance and the other plots of the network.
  ui <- list(u_collapse("Plots",
    u_col(u_label("Output",
      u_check("buildModel_show_histo", "Edge weight histogram"),
      if (data_plots) list(
        u_check("buildModel_show_frequencies", "Code frequencies"),
        u_check("buildModel_show_mosaic", "Mosaic plot")),
      if (a$mixed)
        u_check("buildModel_show_mixed", "Mixed network (transition + co-occurrence)",
                u_combo("buildModel_mixed_type", label = "Display")))),
    u_col(u_label("Network Plot Appearance",
      u_combo("plot_layout"),
      u_check("plot_edge_labels", "Edge labels",
              u_text("plot_edge_label_size", "Size")),
      u_text("plot_minimum", "Hide edges below"),
      u_text("plot_cut", "Cut value (0 = automatic)"),
      u_text("plot_node_size", "Node size (0 = automatic)"),
      u_text("plot_label_size", "Node label size (0 = automatic)")))))

  res <- list(
    r_image("buildModel_plot", vis("buildModel_show_plot"), deps = c("@model", "@look")),
    if (a$mixed) r_state("buildModel_mixedState", deps = "@model"),
    if (a$mixed) r_image("buildModel_mixed", vis("buildModel_show_mixed"),
                         deps = c("@model", "@look", "buildModel_mixed_type"),
                         width = 900L, height = 450L),
    r_image("buildModel_histo", vis("buildModel_show_histo"), deps = "@model"),
    if (data_plots) r_image("buildModel_frequencies", vis("buildModel_show_frequencies"),
                            deps = "@model"),
    if (data_plots) r_image("buildModel_mosaic", vis("buildModel_show_mosaic"),
                            deps = "@model"))
  list(opts = opts, ui = ui, res = res)
}

sec_centrality <- function(a) {
  measure_opts <- unname(Map(function(nm, m) o_bool(paste0("centrality_", nm), m$title, m$default),
                             names(MEASURES), MEASURES))
  params <- c("centrality_loops", "centrality_normalize", "centrality_invert")
  deps <- c("@model", params, MEASURE_OPTS)
  list(
    opts = c(list(
      o_bool("centrality_show_table", "Table", FALSE),
      o_bool("centrality_show_plot", "Plot", FALSE),
      o_list("centrality_plot_style", "Plot style",
             c(lollipop = "Lollipop", bar = "Bar", line = "Line", dot = "Dot",
               heatmap = "Heatmap", distribution = "Distribution")),
      o_bool("centrality_loops", "Include loops", FALSE),
      o_bool("centrality_normalize", "Normalize values", FALSE),
      o_bool("centrality_invert", "Invert weights for path-based measures", TRUE)),
      measure_opts),
    ui = list(u_collapse("Centrality",
      u_col(u_label("Output",
                    u_check("centrality_show_table", "Table"),
                    u_check("centrality_show_plot", "Plot",
                            u_combo("centrality_plot_style", label = "Style"))),
            u_label("Settings",
                    u_check("centrality_loops", "Include loops"),
                    u_check("centrality_normalize", "Normalize values"),
                    u_check("centrality_invert", "Invert weights (path-based measures)"))),
      u_col(u_label("Measures",
                    lapply(names(MEASURES), function(nm)
                      u_check(paste0("centrality_", nm), MEASURES[[nm]]$title)))))),
    res = list(
      r_title("centralityTitle", "Centrality",
              vis("centrality_show_table", "centrality_show_plot")),
      r_table("centralityTable", "Centrality Measures", vis("centrality_show_table"),
              deps = deps, columns = list(group_col(a), r_col("state", "State"))),
      r_image("centrality_plot", vis("centrality_show_plot"),
              deps = c(deps, "centrality_plot_style"))))
}

sec_stability <- function(a) {
  if (!a$resample) return(NULL)
  # Normalising divides by the maximum and leaves the CS correlations as they are.
  deps <- c("@model", "centrality_loops",
            "centrality_invert", MEASURE_OPTS, "stability_iter",
            "stability_threshold", "stability_certainty", "stability_method",
            "@seed")
  list(
    opts = list(
      o_bool("stability_show_table", "Table", FALSE),
      o_bool("stability_show_plot", "Plot", FALSE),
      o_int("stability_iter", "Iterations", 1000, 50, 5000),
      o_num("stability_threshold", "Correlation threshold", 0.7, 0, 1),
      o_num("stability_certainty", "Certainty", 0.95, 0, 1),
      o_list("stability_method", "Correlation",
             c(pearson = "Pearson", spearman = "Spearman", kendall = "Kendall"))),
    ui = list(u_collapse("Reliability & Stability",
      u_col(u_label("Output",
                    u_check("stability_show_table", "Tables"),
                    u_check("stability_show_plot", "Plot"))),
      u_col(u_label("Settings",
                    u_text("stability_iter", "Iterations"),
                    u_text("stability_threshold", "Correlation threshold"),
                    u_text("stability_certainty", "Certainty"),
                    u_combo("stability_method"))))),
    res = list(
      r_title("stabilityTitle", "Reliability & Stability",
              vis("stability_show_table", "stability_show_plot")),
      r_state("stabilityState", deps = deps),
      r_table("stabilityTable", "Centrality Stability (CS-coefficients)",
              vis("stability_show_table"), deps = deps,
              columns = list(group_col(a), r_col("measure", "Measure"),
                             r_col("cs", "CS-coefficient", "number"))),
      r_table("reliabilityTable", "Split-half Reliability",
              vis("stability_show_table"), deps = deps,
              columns = list(group_col(a), r_col("metric", "Metric"),
                             r_col("mean", "Mean", "number"),
                             r_col("sd", "SD", "number"))),
      r_image("stability_plot", vis("stability_show_plot"), deps = deps)))
}

COMMUNITY <- c(louvain = "Louvain", leiden = "Leiden", walktrap = "Walktrap",
               fast_greedy = "Fast greedy", infomap = "Infomap",
               label_propagation = "Label propagation",
               edge_betweenness = "Edge betweenness",
               leading_eigenvector = "Leading eigenvector",
               spinglass = "Spinglass", optimal = "Optimal (small networks)")

sec_community <- function(a) {
  deps <- c("@model", "community_methods", "community_gamma", "@seed")
  list(
    opts = list(
      o_list("community_methods", "Methods", COMMUNITY),
      o_num("community_gamma", "Resolution (gamma)", 1, 0, 100),
      o_bool("community_show_table", "Table", FALSE),
      o_bool("community_show_plot", "Plot", FALSE)),
    ui = list(u_collapse("Communities",
      u_col(u_label("Output",
                    u_check("community_show_table", "Table"),
                    u_check("community_show_plot", "Plot"))),
      u_col(u_label("Method",
                    u_combo("community_methods", label = "Algorithm"),
                    u_text("community_gamma", "Resolution γ (Louvain, Leiden, Spinglass)"))))),
    res = list(
      r_title("communityTitle", "Community Detection",
              vis("community_show_table", "community_show_plot")),
      r_state("communityState", deps = deps),
      r_table("communityTable", "Community Assignments", vis("community_show_table"),
              deps = deps, columns = list(group_col(a), r_col("state", "State"),
                                          r_col("community", "Community", "integer"))),
      r_image("community_plot", vis("community_show_plot"), deps = c(deps, "@look"))))
}

sec_edge_betweenness <- function(a) {
  deps <- c("@model", "edgebetweenness_invert")
  list(
    opts = list(
      o_bool("edgebetweenness_show_table", "Table", FALSE),
      o_bool("edgebetweenness_show_plot", "Plot", FALSE),
      o_bool("edgebetweenness_invert", "Invert weights (strong ties are short)", TRUE)),
    ui = list(u_collapse("Edge Betweenness",
      u_col(u_label("Output",
                    u_check("edgebetweenness_show_table", "Table"),
                    u_check("edgebetweenness_show_plot", "Plot"))),
      u_col(u_label("Settings",
                    u_check("edgebetweenness_invert", "Invert weights (strong ties are short)"))))),
    res = list(
      r_title("edgeBetweennessTitle", "Edge Betweenness",
              vis("edgebetweenness_show_table", "edgebetweenness_show_plot")),
      r_table("edgeBetweennessTable", "Edge Betweenness",
              vis("edgebetweenness_show_table"), deps = deps,
              columns = list(group_col(a), r_col("from", "From"), r_col("to", "To"),
                             r_col("weight", "Weight", "number"),
                             r_col("betweenness", "Edge betweenness", "number"))),
      r_image("edgeBetweenness_plot", vis("edgebetweenness_show_plot"),
              deps = c(deps, "@look"))))
}

sec_cliques <- function(a) {
  deps <- c("@model", "cliques_size", "cliques_threshold")
  list(
    opts = list(
      o_bool("cliques_show_table", "Table", FALSE),
      o_bool("cliques_show_plot", "Plot", FALSE),
      o_bool("cliques_show_complex", "Clique complex summary", FALSE),
      o_int("cliques_size", "Size", 3, 2, 10),
      o_num("cliques_threshold", "Threshold", 0, 0),
      o_int("cliques_max_plot", "Cliques to plot", 9, 1, 50)),
    ui = list(u_collapse("Cliques",
      u_col(u_label("Output",
                    u_check("cliques_show_table", "Table"),
                    u_check("cliques_show_complex", "Clique complex summary"),
                    u_check("cliques_show_plot", "Plot",
                            u_text("cliques_max_plot", "Cliques to plot")))),
      u_col(u_label("Settings",
                    u_text("cliques_size", "Size (nodes)"),
                    u_text("cliques_threshold", "Minimum edge weight"))))),
    res = list(
      r_title("cliquesTitle", "Cliques",
              vis("cliques_show_table", "cliques_show_plot", "cliques_show_complex")),
      r_state("cliquesState", deps = deps),
      r_table("cliquesTable", "Cliques", vis("cliques_show_table"), deps = deps,
              columns = list(group_col(a), r_col("clique", "Clique", "integer"),
                             r_col("members", "Members"),
                             r_col("mean_weight", "Mean weight", "number"),
                             r_col("min_weight", "Minimum weight", "number"))),
      r_table("cliquesComplexTable", "Clique Complex", vis("cliques_show_complex"),
              deps = deps, columns = list(group_col(a), r_col("metric", "Metric"),
                                          r_col("value", "Value", "number"))),
      r_image("cliques_plot", vis("cliques_show_plot"),
              deps = c(deps, "cliques_max_plot", "@look"), width = 750L, height = 750L)))
}

sec_backbone <- function(a) {
  deps <- c("@model", "backbone_method", "backbone_level", "backbone_threshold",
            "backbone_lowest")
  list(
    opts = list(
      o_bool("backbone_show_plot", "Plot", FALSE),
      o_bool("backbone_show_table", "Table", FALSE),
      o_list("backbone_method", "Method",
             c(disparity = "Disparity filter", threshold = "Weight threshold",
               lowest = "Lowest weights")),
      o_num("backbone_level", "Significance level", 0.5, 0, 1),
      o_num("backbone_threshold", "Weight threshold", 0, 0),
      o_num("backbone_lowest", "Lowest proportion", 0.05, 0, 1)),
    ui = list(u_collapse("Backbone",
      u_col(u_label("Output",
                    u_check("backbone_show_table", "Table"),
                    u_check("backbone_show_plot", "Plot"))),
      u_col(u_label("Method",
                    u_radio("backbone_method", "disparity",
                            u_text("backbone_level", "Significance level")),
                    u_radio("backbone_method", "threshold",
                            u_text("backbone_threshold", "Remove weights at or below")),
                    u_radio("backbone_method", "lowest",
                            u_text("backbone_lowest", "Remove lowest proportion")))))),
    res = list(
      r_title("backboneTitle", "Backbone Extraction",
              vis("backbone_show_plot", "backbone_show_table")),
      r_state("backboneState", deps = deps),
      r_table("backboneTable", "Backbone Edges", vis("backbone_show_table"),
              deps = deps, columns = list(group_col(a), r_col("from", "From"),
                                          r_col("to", "To"),
                                          r_col("weight", "Weight", "number"),
                                          r_col("status", "Status"))),
      r_image("backbone_plot", vis("backbone_show_plot"), deps = c(deps, "@look"))))
}

LINKS <- c(resource_allocation = "Resource allocation", adamic_adar = "Adamic-Adar",
           common_neighbors = "Common neighbours", jaccard = "Jaccard",
           preferential_attachment = "Preferential attachment", katz = "Katz")

sec_links <- function(a) {
  deps <- c("@model", "links_method", "links_exclude_existing", "links_threshold",
            "links_top_n")
  list(
    opts = list(
      o_bool("links_show_table", "Table", FALSE),
      o_list("links_method", "Method", LINKS),
      o_bool("links_exclude_existing", "Only absent edges", TRUE),
      o_num("links_threshold", "Treat weights at or below as absent", 0, 0),
      o_int("links_top_n", "Top predictions", 20, 1, 1000)),
    ui = list(u_collapse("Link Prediction",
      u_col(u_label("Output",
                    u_check("links_show_table", "Table",
                            u_text("links_top_n", "Top predictions")))),
      u_col(u_label("Settings",
                    u_combo("links_method"),
                    u_check("links_exclude_existing", "Only absent edges"),
                    u_text("links_threshold", "Treat weights at or below as absent"))))),
    res = list(
      r_title("linksTitle", "Link Prediction", vis("links_show_table")),
      r_table("linksTable", "Predicted Links", vis("links_show_table"), deps = deps,
              columns = list(group_col(a), r_col("from", "From"), r_col("to", "To"),
                             r_col("score", "Score", "number"),
                             r_col("rank", "Rank", "integer"),
                             r_col("existing", "Existing edge")))))
}

sec_rules <- function(a) {
  if (!a$resample) return(NULL)
  deps <- c("@model", "rules_min_support", "rules_min_confidence", "rules_min_lift",
            "rules_max_length")
  list(
    opts = list(
      o_bool("rules_show_table", "Table", FALSE),
      o_bool("rules_show_plot", "Plot", FALSE),
      o_num("rules_min_support", "Minimum support", 0.1, 0, 1),
      o_num("rules_min_confidence", "Minimum confidence", 0.5, 0, 1),
      o_num("rules_min_lift", "Minimum lift", 1, 0),
      o_int("rules_max_length", "Maximum rule length", 5, 2, 10),
      o_int("rules_max_rows", "Max rows", 50, 1, 5000)),
    ui = list(u_collapse("Association Rules",
      u_col(u_label("Output",
                    u_check("rules_show_table", "Table",
                            u_text("rules_max_rows", "Max rows")),
                    u_check("rules_show_plot", "Plot"))),
      u_col(u_label("Settings",
                    u_text("rules_min_support", "Minimum support"),
                    u_text("rules_min_confidence", "Minimum confidence"),
                    u_text("rules_min_lift", "Minimum lift"),
                    u_text("rules_max_length", "Maximum rule length"))))),
    res = list(
      r_title("rulesTitle", "Association Rules", vis("rules_show_table", "rules_show_plot")),
      r_state("rulesState", deps = deps),
      r_table("rulesTable", "Association Rules", vis("rules_show_table"),
              deps = c(deps, "rules_max_rows"),
              columns = list(group_col(a), r_col("antecedent", "If (antecedent)"),
                             r_col("consequent", "Then (consequent)"),
                             r_col("support", "Support", "number"),
                             r_col("confidence", "Confidence", "number"),
                             r_col("lift", "Lift", "number"),
                             r_col("conviction", "Conviction", "number"),
                             r_col("count", "Count", "integer"))),
      r_image("rules_plot", vis("rules_show_plot"), deps = deps, width = 700L)))
}

sec_bootstrap <- function(a) {
  if (!a$resample) return(NULL)
  deps <- c("@model", "bootstrap_iteration", "bootstrap_level", "bootstrap_method",
            "bootstrap_range_low", "bootstrap_range_up", "bootstrap_set_threshold",
            "bootstrap_threshold", "bootstrap_ci_method", "@seed")
  tbl_opts <- c("bootstrap_table_significant_only", "bootstrap_table_show_all",
                "bootstrap_table_max_rows")
  list(
    opts = list(
      o_int("bootstrap_iteration", "Iterations", 1000, 10, 10000),
      o_num("bootstrap_level", "Level", 0.05, 0, 1),
      o_list("bootstrap_method", "Method", c(stability = "Stability",
                                             threshold = "Threshold")),
      o_num("bootstrap_range_low", "Lower", 0.75, 0, 10),
      o_num("bootstrap_range_up", "Upper", 1.25, 0, 10),
      o_bool("bootstrap_set_threshold", "Set threshold", FALSE),
      o_num("bootstrap_threshold", "Threshold", 0, 0),
      o_list("bootstrap_ci_method", "Confidence interval",
             c(percentile = "Percentile", basic = "Basic")),
      o_bool("bootstrap_show_plot", "Plot", FALSE),
      o_bool("bootstrap_show_table", "Table", FALSE),
      o_int("bootstrap_table_max_rows", "Max rows", 20, 1, 1000),
      o_bool("bootstrap_table_show_all", "Show all rows", FALSE),
      o_bool("bootstrap_table_significant_only", "Significant only", FALSE)),
    ui = list(u_collapse("Bootstrap",
      u_col(u_label("Output",
                    u_check("bootstrap_show_table", "Table",
                            u_text("bootstrap_table_max_rows", "Max rows"),
                            u_check("bootstrap_table_show_all", "Show all rows"),
                            u_check("bootstrap_table_significant_only", "Significant only")),
                    u_check("bootstrap_show_plot", "Plot")),
            u_label("Settings",
                    u_text("bootstrap_iteration", "Iterations"),
                    u_text("bootstrap_level", "Significance level"),
                    u_combo("bootstrap_ci_method"))),
      u_col(u_label("Method",
                    u_radio("bootstrap_method", "stability",
                            u_text("bootstrap_range_low", "Lower bound"),
                            u_text("bootstrap_range_up", "Upper bound")),
                    u_radio("bootstrap_method", "threshold",
                            u_check("bootstrap_set_threshold", "Set threshold (else automatic)",
                                    u_text("bootstrap_threshold", "Threshold",
                                           enable = "(bootstrap_method_threshold && bootstrap_set_threshold)"))))))),
    res = list(
      r_title("bootstrapTitle", "Bootstrap", vis("bootstrap_show_table", "bootstrap_show_plot")),
      r_state("bootstrapState", deps = deps),
      r_table("bootstrapTable", "Bootstrap Results", vis("bootstrap_show_table"),
              deps = c(deps, tbl_opts),
              columns = list(group_col(a), r_col("from", "From"), r_col("to", "To"),
                             r_col("weight", "Weight", "number"),
                             r_col("p_value", "p-value", "number", format = "zto,pvalue"),
                             r_col("cr_lower", "CR Lower", "number"),
                             r_col("cr_upper", "CR Upper", "number"),
                             r_col("ci_lower", "CI Lower", "number"),
                             r_col("ci_upper", "CI Upper", "number"),
                             r_col("significant", "Significant"))),
      r_image("bootstrap_plot", vis("bootstrap_show_plot"), deps = c(deps, "@look"))))
}

sec_permutation <- function(a) {
  if (!(a$multi && a$resample)) return(NULL)
  deps <- c("@model", "permutation_iter", "permutation_paired", "permutation_level",
            "permutation_adjust", "@seed")
  list(
    opts = list(
      o_bool("permutation_show_table", "Table", FALSE),
      o_bool("permutation_show_plot", "Plot", FALSE),
      o_int("permutation_iter", "Iterations", 1000, 10, 10000),
      o_bool("permutation_paired", "Paired", FALSE),
      o_num("permutation_level", "Level", 0.05, 0, 1),
      o_list("permutation_adjust", "p-value correction", P_ADJUST),
      o_int("permutation_table_max_rows", "Max rows", 50, 1, 1000),
      o_bool("permutation_table_show_all", "Show all rows", FALSE)),
    ui = list(u_collapse("Permutation Test",
      u_col(u_label("Output",
                    u_check("permutation_show_table", "Table",
                            u_text("permutation_table_max_rows", "Max rows"),
                            u_check("permutation_table_show_all", "Show all rows")),
                    u_check("permutation_show_plot", "Plot"))),
      u_col(u_label("Settings",
                    u_text("permutation_iter", "Iterations"),
                    u_text("permutation_level", "Significance level"),
                    u_combo("permutation_adjust"),
                    u_check("permutation_paired", "Paired"))))),
    res = list(
      r_title("permutationTitle", "Permutation Test",
              vis("permutation_show_table", "permutation_show_plot")),
      r_state("permutationState", deps = deps),
      r_table("permutationTable", "Permutation Test Results", vis("permutation_show_table"),
              deps = c(deps, "permutation_table_max_rows", "permutation_table_show_all"),
              columns = list(r_col("group_comparison", "Comparison"),
                             r_col("edge_name", "Edge"),
                             r_col("diff_true", "True Difference", "number"),
                             r_col("effect_size", "Effect Size", "number"),
                             r_col("p_value", "p-value", "number", format = "zto,pvalue"))),
      r_image("permutation_plot", vis("permutation_show_plot"),
              deps = c(deps, "@look"), height = 600L)))
}

sec_compare <- function(a) {
  if (!a$multi) return(NULL)
  pick <- if (a$cluster) c("compare_cluster_i", "compare_cluster_j")
          else c("compare_group_i", "compare_group_j")
  deps <- c("@model", pick, "compare_scaling")
  pick_opts <- if (a$cluster) list(
    o_int("compare_cluster_i", "First Cluster", 1, 1, 20),
    o_int("compare_cluster_j", "Second Cluster", 2, 1, 20))
  else list(
    o_level("compare_group_i", "First group", "(buildModel_variables_group)"),
    o_level("compare_group_j", "Second group", "(buildModel_variables_group)"))
  pick_ui <- if (a$cluster)
    list(u_text("compare_cluster_i", "First cluster"), u_text("compare_cluster_j", "Second cluster"))
  else list(u_level("compare_group_i"), u_level("compare_group_j"))
  list(
    opts = c(list(
      o_bool("compare_show_summary", "Summary Metrics Table", FALSE),
      o_bool("compare_show_network", "Network Properties Table", FALSE),
      o_bool("compare_show_edges", "Edge Differences Table", FALSE),
      o_bool("compare_show_plot", "Comparison Plot", FALSE),
      o_bool("compare_show_diff_plot", "Network Difference Plot", FALSE)),
      pick_opts, list(
      o_list("compare_scaling", "Scaling",
             c(none = "None", minmax = "MinMax", max = "Max", rank = "Rank")),
      o_list("compare_plot_type", "Plot type",
             c(scatter = "Scatter", heatmap = "Heatmap",
               diff_hist = "Difference histogram", weight_dist = "Weight distribution",
               all = "All")),
      o_int("compare_edges_max_rows", "Max rows", 50, 1, 5000))),
    ui = list(u_collapse(if (a$cluster) "Compare Clusters" else "Compare Groups",
      u_col(u_label("Output",
                    u_check("compare_show_summary", "Summary metrics table"),
                    u_check("compare_show_network", "Network properties table"),
                    u_check("compare_show_edges", "Edge differences table",
                            u_text("compare_edges_max_rows", "Max rows")),
                    u_check("compare_show_plot", "Comparison plot",
                            u_combo("compare_plot_type", label = "Type")),
                    u_check("compare_show_diff_plot", "Network difference plot"))),
      u_col(u_label(if (a$cluster) "Clusters" else "Groups", pick_ui),
            u_label("Settings", u_combo("compare_scaling"))))),
    res = list(
      r_item("compareTitle", "Preformatted", title = "Comparison",
             visible = vis("compare_show_summary", "compare_show_network",
                           "compare_show_edges", "compare_show_plot",
                           "compare_show_diff_plot"), deps = deps),
      r_state("compareState", deps = deps),
      r_table("compareSummaryTable", "Summary Metrics", vis("compare_show_summary"),
              deps = deps, columns = list(r_col("category", "Category"),
                                          r_col("metric", "Metric"),
                                          r_col("value", "Value", "number"))),
      r_table("compareNetworkTable", "Network Properties", vis("compare_show_network"),
              deps = deps, columns = list(r_col("metric", "Metric"),
                                          r_col("group_i", "Network 1", "number"),
                                          r_col("group_j", "Network 2", "number"))),
      r_table("compareEdgeTable", "Edge Differences", vis("compare_show_edges"),
              deps = c(deps, "compare_edges_max_rows"),
              columns = list(r_col("from", "From"), r_col("to", "To"),
                             r_col("weight_i", "Network 1", "number"),
                             r_col("weight_j", "Network 2", "number"),
                             r_col("difference", "Difference", "number"),
                             r_col("relative_difference", "Relative difference", "number"))),
      r_image("compare_plot", vis("compare_show_plot"),
              deps = c(deps, "compare_plot_type"), width = 700L),
      r_image("compare_diff_plot", vis("compare_show_diff_plot"),
              deps = c(deps, "@look"))))
}

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

# The order here is the order of the results; the panel follows it too.
sections <- function(a) compact(list(
  inputs        = sec_inputs(a),
  model         = sec_model(a),
  cluster_diag  = sec_cluster_diag(a),
  visual        = sec_visual(a),
  centrality    = sec_centrality(a),
  edge_btw      = sec_edge_betweenness(a),
  community     = sec_community(a),
  cliques       = sec_cliques(a),
  backbone      = sec_backbone(a),
  links         = sec_links(a),
  rules         = sec_rules(a),
  bootstrap     = sec_bootstrap(a),
  stability     = sec_stability(a),
  permutation   = sec_permutation(a),
  compare       = sec_compare(a)))

# Variables and the network settings stay open at the top; everything else is a
# collapsed section under one of (at most) three headings, each section with
# an "Output" column (what to show) and a "Settings"/"Method" column (what to
# compute).
ui_children <- function(a, secs) {
  ui <- function(key) secs[[key]]$ui
  splice(list(
    ui("inputs"), ui("model"), ui("cluster_diag"),
    u_divider("Network Analysis"),
    ui("visual"), ui("centrality"), ui("edge_btw"), ui("community"), ui("cliques"),
    ui("backbone"), ui("links"), ui("rules"),
    if (a$resample) u_divider("Confirmatory Testing"),
    ui("bootstrap"), ui("stability"), ui("permutation"),
    if (a$multi) u_divider(if (a$cluster) "Cluster Comparison" else "Group Comparison"),
    ui("compare")))
}

resolve_deps <- function(deps, a, declared) {
  expanded <- unlist(lapply(deps, function(d) switch(d,
    "@model" = model_opts(a),
    "@look"  = LOOK_OPTS,
    "@seed"  = SEED_OPTS,
    d)), use.names = FALSE)
  expanded <- unique(expanded)
  unknown <- setdiff(expanded, declared)
  if (length(unknown))
    stop(sprintf("%s: clearWith refers to undeclared option(s): %s",
                 a$name, paste(unknown, collapse = ", ")), call. = FALSE)
  expanded
}

result_yaml <- function(item, a, declared) {
  out <- list(name = item$name, title = item$title, type = item$type)
  if (!is.null(item$width))  out$width  <- item$width
  if (!is.null(item$height)) out$height <- item$height
  out$visible <- item$visible
  if (isTRUE(item$render)) out$renderFun <- paste0(".render_", item$name)
  if (!is.null(item$columns)) out$columns <- item$columns
  deps <- resolve_deps(item$deps, a, declared)
  if (length(deps)) out$clearWith <- as.list(deps)
  compact(out)
}

ui_nodes <- function(nodes)
  do.call(c, lapply(nodes, function(n)
    c(list(n), if (!is.null(n$children)) ui_nodes(n$children))))

# Option names bound by UI controls (radio buttons bind their optionName).
ui_bound_options <- function(ui) {
  binders <- c("CheckBox", "TextBox", "ComboBox", "LevelSelector",
               "VariablesListBox", "RadioButton")
  nodes <- Filter(function(n) isTRUE(n$type %in% binders), ui_nodes(ui))
  unique(vapply(nodes, function(n) n$optionName %||% n$name, character(1)))
}

render_items <- function(secs) {
  items <- unname(compact(unlist(lapply(secs, `[[`, "res"), recursive = FALSE)))
  vapply(Filter(function(i) isTRUE(i$render), items), `[[`, character(1), "name")
}

# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

as_logical_yaml <- function(x) {
  r <- ifelse(x, "true", "false")
  class(r) <- "verbatim"
  r
}

# yaml::as.yaml escapes astral-plane characters (the bold dividers) as \UXXXXXXXX
# inside double quotes; turn them back into UTF-8 so the files stay readable.
unescape_astral <- function(txt) {
  hits <- gregexpr("\\\\U[0-9A-Fa-f]{8}", txt)
  regmatches(txt, hits) <- lapply(regmatches(txt, hits), function(m)
    vapply(m, function(code) intToUtf8(strtoi(substring(code, 3), 16L)), character(1)))
  txt
}

expand_dividers <- function(txt) {
  hits <- gregexpr("DIVIDER<[^>]*>", txt)
  regmatches(txt, hits) <- lapply(regmatches(txt, hits), function(m)
    vapply(m, function(tag)
      paste("\u2501\u2501", bold_sans(sub("^DIVIDER<(.*)>$", "\\1", tag)), "\u2501\u2501"),
      character(1)))
  txt
}

emit <- function(x) {
  txt <- yaml::as.yaml(x, indent.mapping.sequence = TRUE,
                       handlers = list(logical = as_logical_yaml))
  expand_dividers(unescape_astral(txt))
}

write_lines <- function(text, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(text, con, useBytes = TRUE)
}

header <- function(what)
  sprintf("# Generated by tools/generate_jamovi.R (CNA %s) -- edit the generator, not this file.\n# %s\n", VERSION, what)

write_analysis <- function(a) {
  secs <- sections(a)
  opts <- unname(unlist(lapply(secs, `[[`, "opts"), recursive = FALSE))
  names_opts <- vapply(opts, `[[`, character(1), "name")
  dup <- names_opts[duplicated(names_opts)]
  if (length(dup)) stop(a$name, ": duplicate options: ", paste(dup, collapse = ", "))
  declared <- names_opts

  # Every UI control binds a declared option, and every option has a control
  # (jamovi's compiler adds a stray box for an option without one).
  ui <- ui_children(a, secs)
  bound <- ui_bound_options(ui)
  stray <- setdiff(bound, declared)
  if (length(stray)) stop(a$name, ": UI binds undeclared option(s): ",
                          paste(stray, collapse = ", "))
  orphan <- setdiff(declared, c("data", bound))
  if (length(orphan)) stop(a$name, ": option(s) without a UI control: ",
                           paste(orphan, collapse = ", "))

  a_yaml <- list(name = a$name, title = a$title, menuGroup = "Co-occurrence Network Analysis",
                 menuSubgroup = a$subgroup, menuTitle = a$menu, version = VERSION,
                 jas = "1.2", options = opts)
  results <- unname(compact(unlist(lapply(secs, `[[`, "res"), recursive = FALSE)))
  res_names <- vapply(results, `[[`, character(1), "name")
  if (anyDuplicated(res_names)) stop(a$name, ": duplicate result names")
  r_yaml <- list(name = a$name, title = a$title, jrs = "1.1",
                 items = lapply(results, result_yaml, a = a, declared = declared))
  u_yaml <- list(title = a$title, name = a$name, jus = "3.0", stage = 0L,
                 compilerMode = "tame", children = ui)

  write_lines(c(header("Options"), "---", emit(a_yaml), "..."),
              file.path(ROOT, "jamovi", paste0(a$name, ".a.yaml")))
  write_lines(c(header("Results"), "---", emit(r_yaml), "..."),
              file.path(ROOT, "jamovi", paste0(a$name, ".r.yaml")))
  write_lines(c(header("User interface"), emit(u_yaml)),
              file.path(ROOT, "jamovi", paste0(a$name, ".u.yaml")))
  write_backend(a, render_items(secs))
  invisible(list(options = declared, results = res_names))
}

write_backend <- function(a, renders) {
  methods <- c(
    sprintf('    .run = function() cna_run_analysis(self, "%s")', a$kind),
    sprintf('    .render_%s = function(image, ...) cna_render(self, "%s", image)',
            renders, renders))
  code <- c(
    sprintf("# Generated by tools/generate_jamovi.R (CNA %s) -- do not edit by hand.", VERSION),
    sprintf("# %s: thin backend. The analysis itself lives in the shared engine", a$title),
    "# (R/cna_*.R); this class only routes .run() and the plot renderers to it.",
    "",
    sprintf("%sClass <- if (requireNamespace(\"jmvcore\", quietly = TRUE)) R6::R6Class(", a$name),
    sprintf('  "%sClass",', a$name),
    sprintf("  inherit = %sBase,", a$name),
    "  private = list(",
    paste(methods, collapse = ",\n"),
    "  )",
    ")")
  write_lines(code, file.path(ROOT, "R", paste0(a$name, ".b.R")))
}

write_manifest <- function() {
  manifest <- list(
    title = "Co-occurrence Network Analysis", name = "CNA", version = VERSION,
    jms = "1.0",
    authors = list("Dylan Girault", "Mohammed Saqr", "Santtu Tikka", "Sonsoles López-Pernas"),
    maintainer = "Sonsoles López-Pernas <sonsoles.lopez@uef.fi>",
    date = DATE, type = "R",
    description = paste(
      "Co-occurrence Network Analysis (CNA) of event (long-format), one-hot",
      "(binary) and delimited (keyword) data, built on Nestimate and cograph.",
      "Single networks, group comparison and clustering; centrality and its",
      "stability, community detection, edge betweenness, cliques and the clique",
      "complex, backbone extraction, link prediction, association rules,",
      "bootstrap validation, permutation tests and network comparison."),
    analyses = lapply(ANALYSES, function(a) list(
      title = a$title, name = a$name, ns = "CNA",
      menuGroup = "Co-occurrence Network Analysis", menuSubgroup = a$subgroup,
      menuTitle = a$menu, category = "analyses")),
    usesNative = TRUE, minApp = "1.8.0",
    datasets = list(
      list(name = "Regulation_long", path = "Regulation_long.csv",
           description = paste("Real. Group regulation actions in collaborative learning,",
                               "long format: Actor, Achievers, Group, Time and Action.",
                               "For the Event Data analyses."),
           tags = list("CNA", "event data", "long format", "real")),
      list(name = "Regulation_sessions", path = "Regulation_sessions.csv",
           description = paste("Real. The Regulation_long sessions as action lists",
                               "(one row per session, actions separated by ';').",
                               "For the Delimited Data analysis (Item Columns: Actions)."),
           tags = list("CNA", "delimited", "real")),
      list(name = "SRL_events", path = "SRL_events.csv",
           description = paste("Simulated (simulab). Self-regulated learning event log:",
                               "8 actions, 200 students in two groups with different",
                               "learning cycles (Strategic: plan-read-summarise-monitor-",
                               "evaluate-reflect; Surface: read-review-practice), courses,",
                               "achievement, Time and Order. For the Event Data analyses,",
                               "group comparison and the mixed network."),
           tags = list("CNA", "event data", "self-regulated learning", "simulated")),
      list(name = "SRL_onehot", path = "SRL_onehot.csv",
           description = paste("Simulated (simulab). The SRL_events log coded 0/1: one row",
                               "per episode of four consecutive events, one column per",
                               "action, with Actor, Group, Course, Achievement and Episode.",
                               "For the Binary Data analyses and the mixed network",
                               "(Actor = Actor, window size 1)."),
           tags = list("CNA", "one-hot", "binary", "self-regulated learning", "simulated")),
      list(name = "LearningActivities", path = "LearningActivities.csv",
           description = paste("Simulated (Nestimate::learning_activities). 0/1 indicators",
                               "of six online learning activities for 200 students over 30",
                               "time points. For the Binary Data analyses",
                               "(Actor = Student) and clustering."),
           tags = list("CNA", "one-hot", "binary", "simulated")),
      list(name = "HumanAI_turns", path = "HumanAI_turns.csv",
           description = paste("Real (Saqr, 2026; Nestimate::human_long and ai_long).",
                               "Human-AI pair-programming turns from 429 sessions: the",
                               "turn's codes as a '; ' list (Codes) and as 17 0/1 columns,",
                               "with Session, Project, Speaker and Turn. For the Delimited",
                               "and the Binary Data analyses (Group = Speaker)."),
           tags = list("CNA", "delimited", "one-hot", "real"))))
  write_lines(c("---", emit(manifest), "..."), file.path(ROOT, "jamovi", "0000.yaml"))
}

if (sys.nframe() == 0L) {
  info <- lapply(ANALYSES, write_analysis)
  write_manifest()
  cat(sprintf("Wrote %d analyses (%s), version %s.\n", length(ANALYSES),
              paste(vapply(ANALYSES, `[[`, character(1), "name"), collapse = ", "),
              VERSION))
  invisible(Map(function(a, i) cat(sprintf("  %-17s %3d options, %3d results\n",
                                           a$name, length(i$options), length(i$results))),
                ANALYSES, info))
}
