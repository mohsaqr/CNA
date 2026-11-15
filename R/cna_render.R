# Plot renderers. Each result image's renderFun calls cna_render(self, key,
# image), which looks up the renderer for that result and draws from the
# analysis's cached state (the network(s) and the analyses' engine objects).
# A renderer returns FALSE when there is nothing to draw yet.

#' Network-plot appearance options -> cograph::splot() arguments.
#'
#' Edge labels are switched differently by input type (verified on rendered
#' SVG, cograph 2.7.12): a network object labels its edges by default and
#' ignores `edge_labels = FALSE`, so hiding uses `edge_label_style = "none"`;
#' a bare matrix is unlabelled by default, so showing needs `edge_labels =
#' TRUE`. Digits are never forced: `weight_digits` turns count labels into
#' "118.00", while cograph's own default prints counts as integers.
#' "0 = automatic" sizes are left to cograph.
#' @param matrix TRUE when the plotted object is a bare weight matrix.
#' @noRd
.cna_look <- function(options, matrix = FALSE) {
  positive <- function(x) if (!is.null(x) && is.finite(x) && x > 0) x
  show_labels <- !isFALSE(options$plot_edge_labels)
  label_args <- if (matrix) list(edge_labels = show_labels)
                else if (!show_labels) list(edge_label_style = "none")
  c(Filter(Negate(is.null), list(
      layout = options$plot_layout %||% "oval",
      edge_label_size = positive(options$plot_edge_label_size),
      minimum = positive(options$plot_minimum),
      cut = positive(options$plot_cut),
      node_size = positive(options$plot_node_size),
      label_size = positive(options$plot_label_size))),
    label_args)
}

#' Appearance options -> cograph::plot_mixed_network() arguments (it names
#' the edge cut-off `threshold` and draws its own edge labels).
#' @noRd
.cna_mixed_look <- function(options) {
  positive <- function(x) if (!is.null(x) && is.finite(x) && x > 0) x
  Filter(Negate(is.null), list(
    layout = options$plot_layout %||% "oval",
    edge_labels = !isFALSE(options$plot_edge_labels),
    edge_label_size = positive(options$plot_edge_label_size),
    node_size = positive(options$plot_node_size),
    threshold = positive(options$plot_minimum)))
}

#' Decimal places for edge labels: none for count (integer) weights, two
#' otherwise (scaled or similarity weights).
#' @noRd
.cna_digits <- function(weights) {
  w <- as.numeric(unlist(weights))
  w <- w[is.finite(w)]
  if (!length(w) || all(abs(w - round(w)) < sqrt(.Machine$double.eps))) 0L else 2L
}

.cna_model_weights <- function(model) lapply(.cna_nets(model), `[[`, "weights")

#' cograph::splot() with the appearance options plus extra arguments.
#' @noRd
.cna_splot <- function(x, options, ...) {
  do.call(cograph::splot, c(list(x), .cna_look(options, matrix = is.matrix(x)), list(...)))
}

#' Rows x columns for n panels.
#' @noRd
.cna_grid <- function(n) {
  ncol <- ceiling(sqrt(n))
  c(ceiling(n / ncol), ncol)
}

#' Draw `n` base-graphics panels with `draw(i)`; par() is restored after.
#' @noRd
.cna_panels <- function(n, draw) {
  op <- graphics::par(mfrow = .cna_grid(n))
  on.exit(graphics::par(op), add = TRUE)
  invisible(lapply(seq_len(n), draw))
}

#' An empty panel carrying a message (used when a plot cannot be drawn).
#' @noRd
.cna_plot_message <- function(title, message) {
  graphics::plot.new()
  graphics::title(main = title)
  graphics::text(0.5, 0.5, paste(strwrap(message, 60), collapse = "\n"), cex = 0.9)
}

#' A bootstrap result with the self co-occurrences (diagonal) set to 0, for
#' drawing only. cograph 2.7.12's splot.net_bootstrap styles the off-diagonal
#' edges of an undirected network but draws its self-loops as well, and stops
#' ("edge_ci_scale must be length 1 or 15, not 10"). The bootstrap table keeps
#' the diagonal rows.
#' @noRd
.cna_boot_without_loops <- function(bs) {
  if (inherits(bs, "net_bootstrap_group")) {
    out <- lapply(bs, .cna_boot_without_loops)
    attributes(out) <- attributes(bs)
    return(out)
  }
  zero <- function(m) {
    if (is.matrix(m)) diag(m) <- 0
    m
  }
  if (!is.null(bs$original$weights)) bs$original$weights <- zero(bs$original$weights)
  if (!is.null(bs$model$weights)) bs$model$weights <- zero(bs$model$weights)
  bs
}

#' Print a ggplot, a list of ggplots (stacked with patchwork) or a gtable.
#' @noRd
.cna_print_gg <- function(p) {
  if (inherits(p, "gtable") || inherits(p, "grob")) {
    grid::grid.newpage()
    grid::grid.draw(p)
  } else if (is.list(p) && !inherits(p, "gg") && !inherits(p, "ggplot")) {
    print(patchwork::wrap_plots(p, ncol = 1L))
  } else {
    print(p)
  }
  TRUE
}

.cna_state <- function(self, name) {
  item <- self$results[[name]]
  if (is.null(item)) NULL else item$state
}

.cna_model_of <- function(self) .cna_state(self, "buildModelContent")

# Each renderer: function(self) -> TRUE (drew) / FALSE (nothing to draw).
.cna_renderers <- list(

  buildModel_plot = function(self) {
    model <- .cna_model_of(self)
    if (is.null(model)) return(FALSE)
    .cna_splot(model, self$options)
    TRUE
  },

  buildModel_mixed = function(self) {
    pairs <- .cna_state(self, "buildModel_mixedState")
    if (is.null(pairs)) return(FALSE)
    o <- self$options
    if (identical(o$buildModel_mixed_type, "overlay")) {
      # one panel per network: transitions (arrows) over co-occurrence
      .cna_panels(length(pairs), function(i)
        do.call(cograph::plot_mixed_network, c(
          list(sym_matrix = pairs[[i]]$cooccurrence$weights,
               asym_matrix = pairs[[i]]$transition$weights,
               title = if (nzchar(names(pairs)[[i]])) names(pairs)[[i]]),
          .cna_mixed_look(o))))
      return(TRUE)
    }
    # side by side: a two-column grid, Transition | Co-occurrence per network
    # (a netobject_group of the pairs, as cograph's splot.wtna_mixed builds)
    # unname: unlist() would otherwise prefix each panel with its group name
    panels <- unlist(unname(Map(function(p, g)
      stats::setNames(list(p$transition, p$cooccurrence),
                      trimws(paste(g, c("Transition", "Co-occurrence")))),
      pairs, names(pairs))), recursive = FALSE)
    grid <- structure(panels, class = "netobject_group")
    # TNA styling named explicitly: cograph 2.7.12 does not list wtna's
    # transition method ("wtna_transition") as TNA-family and would draw that
    # directed network undirected (no arrows, one triangle of each pair)
    do.call(cograph::splot, c(list(grid), .cna_look(o),
                              list(ncol = 2L, tna_styling = TRUE, psych_styling = FALSE)))
    TRUE
  },

  buildModel_histo = function(self) {
    model <- .cna_model_of(self)
    if (is.null(model)) return(FALSE)
    nets <- .cna_nets(model)
    .cna_panels(length(nets), function(i) {
      w <- as.numeric(nets[[i]]$weights)
      w <- w[is.finite(w) & w != 0]
      graphics::hist(w, main = if (nzchar(names(nets)[[i]])) names(nets)[[i]]
                     else "Histogram of edge weights",
                     xlab = "Edge weight", ylab = "Frequency",
                     col = "#56B4E9", border = "white")
    })
    TRUE
  },

  buildModel_frequencies = function(self) {
    model <- .cna_model_of(self)
    if (is.null(model)) return(FALSE)
    # plot_state_frequencies() also prints a frequency summary to the console;
    # the plot is what the result shows, so that text is captured, not shown.
    utils::capture.output(print(Nestimate::plot_state_frequencies(model)))
    TRUE
  },

  buildModel_mosaic = function(self) {
    model <- .cna_model_of(self)
    if (is.null(model)) return(FALSE)
    .cna_print_gg(Nestimate::mosaic_plot(model))
  },

  centrality_plot = function(self) {
    cent <- .cna_state(self, "centralityTable")
    if (is.null(cent) || !nrow(cent)) return(FALSE)
    style <- self$options$centrality_plot_style %||% "lollipop"
    measures <- setdiff(names(cent), c("group", "state"))
    groups <- unique(cent$group)
    as_cograph <- function(df) {
      out <- df[c("state", measures)]
      names(out)[1L] <- "node"
      out
    }
    if (identical(style, "distribution")) {
      .cna_panels(length(groups) * length(measures), function(k) {
        g <- groups[[(k - 1L) %/% length(measures) + 1L]]
        m <- measures[[(k - 1L) %% length(measures) + 1L]]
        cograph::plot_centrality_distribution(
          as_cograph(cent[cent$group == g, , drop = FALSE]), measure = m,
          main = trimws(paste(g, .cna_measure_titles[[m]])))
      })
      return(TRUE)
    }
    plots <- lapply(groups, function(g) {
      df <- as_cograph(cent[cent$group == g, , drop = FALSE])
      p <- if (identical(style, "heatmap")) cograph::plot_centrality_heatmap(df)
           else cograph::plot_centrality(df, style = style)
      if (nzchar(g)) p + ggplot2::ggtitle(g) else p
    })
    .cna_print_gg(if (length(plots) == 1L) plots[[1L]] else plots)
  },

  stability_plot = function(self) {
    st <- .cna_state(self, "stabilityState")
    if (is.null(st)) return(FALSE)
    if (.cna_failed(st$cs)) {
      .cna_plot_message("Centrality stability", cna_stability_failures(st)[["cs"]])
      return(TRUE)
    }
    # Nestimate plots one net_stability; a group result (groups, clusters) has
    # no plot method, so each network is plotted and the panels stacked.
    parts <- .cna_parts(st$cs, inherits(st$cs, "net_stability_group"))
    parts <- Filter(function(cs) length(cs$cs) > 0L, parts)
    if (!length(parts)) {
      .cna_plot_message("Centrality stability", paste(
        "No selected measure varies across nodes, so there is no stability to plot."))
      return(TRUE)
    }
    plots <- Map(function(cs, g) {
      p <- plot(cs)
      if (nzchar(g)) p + ggplot2::ggtitle(g) else p
    }, parts, names(parts))
    .cna_print_gg(if (length(plots) == 1L) plots[[1L]] else unname(plots))
  },

  community_plot = function(self) {
    coms <- .cna_state(self, "communityState")
    if (is.null(coms)) return(FALSE)
    .cna_panels(length(coms), function(i)
      .cna_splot(coms[[i]], self$options,
                 title = if (nzchar(names(coms)[[i]])) names(coms)[[i]]))
    TRUE
  },

  edgeBetweenness_plot = function(self) {
    model <- .cna_model_of(self)
    if (is.null(model)) return(FALSE)
    ebs <- cna_edge_betweenness_compute(model, isTRUE(self$options$edgebetweenness_invert))
    .cna_panels(length(ebs), function(i)
      .cna_splot(ebs[[i]], self$options,
                 title = trimws(paste(names(ebs)[[i]], "Edge betweenness"))))
    TRUE
  },

  cliques_plot = function(self) {
    model <- .cna_model_of(self)
    complexes <- .cna_state(self, "cliquesState")
    if (is.null(model) || is.null(complexes)) return(FALSE)
    size <- self$options$cliques_size
    cliques <- cna_cliques_df(complexes, model, size)
    if (is.null(cliques)) {
      .cna_plot_message("Cliques", sprintf("No cliques of %d nodes at this threshold.", size))
      return(TRUE)
    }
    shown <- do.call(rbind, lapply(split(cliques, factor(cliques$group, unique(cliques$group))),
                                   utils::head, self$options$cliques_max_plot))
    nets <- .cna_nets(model)
    .cna_panels(nrow(shown), function(i) {
      row <- shown[i, ]
      nodes <- strsplit(row$members, ", ", fixed = TRUE)[[1L]]
      # by position: a single network's name is "", which [[ ]] cannot look up
      net <- nets[[match(row$group, names(nets))]]
      w <- unclass(net$weights)[nodes, nodes, drop = FALSE]
      .cna_splot(w, self$options, title = trimws(paste(row$group, "Clique", row$clique)))
    })
    TRUE
  },

  backbone_plot = function(self) {
    backbones <- .cna_state(self, "backboneState")
    if (is.null(backbones)) return(FALSE)
    .cna_panels(length(backbones), function(i)
      .cna_splot(backbones[[i]]$backbone, self$options,
                 title = trimws(paste(names(backbones)[[i]], "Backbone"))))
    TRUE
  },

  rules_plot = function(self) {
    rules <- .cna_state(self, "rulesState")
    if (is.null(rules)) return(FALSE)
    with_rules <- Filter(function(r) isTRUE(r$n_rules > 0), rules)
    if (!length(with_rules)) {
      .cna_plot_message("Association rules", "No rules meet the support, confidence and lift limits.")
      return(TRUE)
    }
    plots <- Map(function(r, g) {
      p <- plot(r)
      if (nzchar(g)) p + ggplot2::ggtitle(g) else p
    }, with_rules, names(with_rules))
    .cna_print_gg(if (length(plots) == 1L) plots[[1L]] else unname(plots))
  },

  bootstrap_plot = function(self) {
    bs <- .cna_state(self, "bootstrapState")
    model <- .cna_model_of(self)
    if (is.null(bs) || is.null(model)) return(FALSE)
    # labels with significance stars; counts printed without ".00"
    labels <- if (!isFALSE(self$options$plot_edge_labels))
      list(edge_label_digits = .cna_digits(.cna_model_weights(model)))
    do.call(cograph::splot, c(list(.cna_boot_without_loops(bs)), .cna_look(self$options),
                              labels))
    TRUE
  },

  permutation_plot = function(self) {
    pt <- .cna_state(self, "permutationState")
    model <- .cna_model_of(self)
    if (is.null(pt) || is.null(model)) return(FALSE)
    # cograph draws nothing at all for a comparison without significant edges
    if (all(vapply(pt, function(cmp) all(cmp$diff_sig == 0), logical(1)))) {
      .cna_plot_message("Permutation test", sprintf(
        "No edge differs significantly between any pair of networks (alpha = %s%s).",
        format(pt[[1]]$alpha),
        if (!identical(pt[[1]]$adjust %||% "none", "none"))
          paste0(", ", pt[[1]]$adjust, "-adjusted") else ""))
      return(TRUE)
    }
    # splot.net_permutation formats its star labels with "%.2f"; the template
    # path with edge_label_digits shows count differences without ".00".
    labels <- if (!isFALSE(self$options$plot_edge_labels))
      list(edge_label_template = "{est}{stars}",
           edge_label_digits = .cna_digits(.cna_model_weights(model)))
    do.call(cograph::splot, c(list(pt), .cna_look(self$options), labels))
    TRUE
  },

  compare_plot = function(self) {
    state <- .cna_state(self, "compareState")
    if (is.null(state)) return(FALSE)
    .cna_print_gg(plot(state$result, type = self$options$compare_plot_type %||% "scatter"))
  },

  compare_diff_plot = function(self) {
    state <- .cna_state(self, "compareState")
    model <- .cna_model_of(self)
    if (is.null(state) || is.null(model)) return(FALSE)
    look <- .cna_look(self$options)
    do.call(cograph::plot_difference, c(list(model, i = state$i, j = state$j,
                                             title = paste(state$group_i, "vs", state$group_j)),
                                        look))
    TRUE
  },

  clusterDiag_plot = function(self) {
    model <- .cna_model_of(self)
    if (is.null(model)) return(FALSE)
    .cna_print_gg(plot(Nestimate::cluster_diagnostics(attr(model, "clustering"))))
  },

  clusterChoice_plot = function(self) {
    cc <- .cna_state(self, "clusterChoiceState")
    if (is.null(cc)) return(FALSE)
    .cna_print_gg(plot(cc))
  }
)

#' Draw result `key` of an analysis (called from each image's renderFun).
#' A failure is drawn into the image rather than raised, so one plot cannot
#' stop the others.
#' @noRd
cna_render <- function(self, key, image = NULL) {
  renderer <- .cna_renderers[[key]]
  if (is.null(renderer))
    .cna_stop(sprintf("no renderer for result '%s'", key), "cna_engine_contract")
  tryCatch(isTRUE(renderer(self)), error = function(e) {
    .cna_plot_message("This plot could not be drawn", conditionMessage(e))
    TRUE
  })
}
