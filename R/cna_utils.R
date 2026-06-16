# Shared utility functions for CNA modules (CNASingle, CNAGroup, CNACluster)

#' Generate dynamic instructions HTML based on module type
#' @param module_type Combined string: "event_single", "event_cluster",
#'   "binary_single", or "binary_cluster"
#' @return HTML string for the instructions
cna_instructions_html <- function(module_type) {
  fields_html <- switch(module_type,
    "event_single"  = '<li><b>Event Data</b>: Assign <b>Action</b> (required), <b>Actor</b>, <b>Time</b>, <b>Order</b> (optional). Add <b>Group Variable</b> to enable group mode.</li>',
    "event_cluster" = '<li><b>Event Data</b>: Assign <b>Action</b> (required), <b>Actor</b> (required for clustering), <b>Time</b>, <b>Order</b> (optional).</li>',
    "binary_single" = '<li><b>Binary Data</b>: Assign <b>One-Hot Columns</b> (at least 2), <b>Actor</b>, <b>Session</b> (optional). Add <b>Group Variable</b> to enable group mode.</li>',
    "binary_cluster" = '<li><b>Binary Data</b>: Assign <b>One-Hot Columns</b> (at least 2), <b>Actor</b>, <b>Session</b> (optional).</li>'
  )

  mode_html <- switch(module_type,
    "event_single"  = '<li><b>Single/Group Network</b>: Builds one network from all data, or separate networks per group if Group Variable is set.</li>',
    "event_cluster" = '<li><b>Clustering</b>: Automatically clusters sequences and builds per-cluster networks. Check "Run Clustering Analysis" to start.</li>',
    "binary_single" = '<li><b>Single/Group Network</b>: Builds one network from all data, or separate networks per group if Group Variable is set.</li>',
    "binary_cluster" = '<li><b>Clustering</b>: Automatically clusters rows and builds per-cluster networks. Check "Run Clustering Analysis" to start.</li>'
  )

  link_html <- switch(module_type,
    "event_single"  = '<li>Learn more: <a href="https://lamethods.org/book2/chapters/ch15-tna/ch15-tna.html" target="_blank">TNA Tutorial</a> | <a href="https://lamethods.org/book2/chapters/ch17-tna-clusters/ch17-tna-clusters.html" target="_blank">Group TNA</a></li>',
    "event_cluster" = '<li>Learn more: <a href="https://lamethods.org/book2/chapters/ch17-tna-clusters/ch17-tna-clusters.html" target="_blank">Group TNA</a></li>',
    "binary_single" = '<li>Learn more: <a href="https://lamethods.org/book2/chapters/ch15-tna/ch15-tna.html" target="_blank">TNA Tutorial</a> | <a href="https://lamethods.org/book2/chapters/ch17-tna-clusters/ch17-tna-clusters.html" target="_blank">Group TNA</a></li>',
    "binary_cluster" = '<li>Learn more: <a href="https://lamethods.org/book2/chapters/ch17-tna-clusters/ch17-tna-clusters.html" target="_blank">Group TNA</a></li>'
  )

  paste0(
    '<div style="border: 2px solid #e6f4fe; border-radius: 15px; padding: 15px; background-color: #e6f4fe; margin-top: 10px;">',
    '<div style="text-align:justify;"><ul>',
    mode_html,
    fields_html,
    link_html,
    '</ul></div></div>'
  )
}

#' Prepare event (long-format) data for CNA
#' @param data The raw data frame from jamovi
#' @param options The jamovi options object
#' @return A list with prepData (from tna::prepare_data) and seq_data
cna_prepare_event_data <- function(data, options) {
  copyData <- data
  copyData[[options$buildModel_variables_long_action]] <-
    as.character(copyData[[options$buildModel_variables_long_action]])

  if (!is.null(options$buildModel_variables_long_time)) {
    copyData[[options$buildModel_variables_long_time]] <-
      as.POSIXct(copyData[[options$buildModel_variables_long_time]])
  }
  if (!is.null(options$buildModel_variables_actor)) {
    copyData[[options$buildModel_variables_actor]] <-
      as.character(copyData[[options$buildModel_variables_actor]])
  }
  if (!is.null(options$buildModel_variables_long_order)) {
    copyData[[options$buildModel_variables_long_order]] <-
      as.character(copyData[[options$buildModel_variables_long_order]])
  }

  threshold <- options$buildModel_threshold

  columnToUseLong <- c(
    options$buildModel_variables_long_time,
    options$buildModel_variables_actor,
    options$buildModel_variables_long_action,
    options$buildModel_variables_long_order
  )

  longData <- copyData[columnToUseLong]

  args_prepare_data <- list(
    data = longData,
    actor = options$buildModel_variables_actor,
    time = options$buildModel_variables_long_time,
    action = options$buildModel_variables_long_action,
    time_threshold = threshold,
    order = options$buildModel_variables_long_order
  )
  args_prepare_data <- args_prepare_data[
    !vapply(args_prepare_data, is.null, logical(1))
  ]

  prepData <- do.call(tna::prepare_data, args_prepare_data)
  seq_data <- as.data.frame(prepData$sequence_data)

  list(prepData = prepData, seq_data = seq_data, longData = longData)
}

#' Prepare binary (one-hot) data for CNA
#' @param data The raw data frame from jamovi
#' @param options The jamovi options object
#' @param group_col Optional group column name (for CNAGroup)
#' @return A list with seq_data and optional group_vector
cna_prepare_binary_data <- function(data, options, group_col = NULL) {
  df <- as.data.frame(data)
  onehot_cols <- options$buildModel_variables_onehot
  actor_col <- options$buildModel_variables_actor
  session_col <- options$buildModel_variables_session
  window_size <- options$buildModel_window

  # Convert one-hot: replace 1 with column name, 0 with NA
  seq_data <- as.data.frame(lapply(onehot_cols, function(col) {
    ifelse(df[[col]] == 1, col, NA)
  }))
  colnames(seq_data) <- onehot_cols

  # Track group variable if provided
  if (!is.null(group_col)) {
    seq_data$..group_var.. <- as.character(df[[group_col]])
  }

  # Aggregate by actor + session/window if provided
  if (!is.null(actor_col) || !is.null(session_col)) {
    if (!is.null(actor_col) && !is.null(session_col)) {
      session_numeric <- as.numeric(factor(df[[session_col]]))
      window_id <- floor((session_numeric - 1) / window_size)
      group_id <- paste(df[[actor_col]], window_id, sep = "_")
    } else if (!is.null(actor_col)) {
      row_window <- floor((seq_len(nrow(df)) - 1) / window_size)
      group_id <- paste(df[[actor_col]], row_window, sep = "_")
    } else {
      session_numeric <- as.numeric(factor(df[[session_col]]))
      group_id <- floor((session_numeric - 1) / window_size)
    }

    seq_data$..group_id.. <- group_id
    cols_to_agg <- if (!is.null(group_col)) {
      c(onehot_cols, "..group_var..")
    } else {
      onehot_cols
    }
    seq_data <- aggregate(
      seq_data[cols_to_agg],
      by = list(..group_id.. = seq_data$..group_id..),
      FUN = function(x) {
        non_na <- na.omit(x)
        if (length(non_na) > 0) non_na[1] else NA
      }
    )
    seq_data$..group_id.. <- NULL
  } else if (window_size > 1) {
    row_window <- floor((seq_len(nrow(df)) - 1) / window_size)
    seq_data$..group_id.. <- row_window
    cols_to_agg <- if (!is.null(group_col)) {
      c(onehot_cols, "..group_var..")
    } else {
      onehot_cols
    }
    seq_data <- aggregate(
      seq_data[cols_to_agg],
      by = list(..group_id.. = seq_data$..group_id..),
      FUN = function(x) {
        non_na <- na.omit(x)
        if (length(non_na) > 0) non_na[1] else NA
      }
    )
    seq_data$..group_id.. <- NULL
  }

  # Extract group vector and remove from seq_data
  group_vector <- NULL
  if (!is.null(group_col)) {
    group_vector <- seq_data$..group_var..
    seq_data$..group_var.. <- NULL
  }

  list(seq_data = seq_data, group_vector = group_vector)
}

#' Set multi-panel layout
#' @param n Number of panels
cna_set_panel_layout <- function(n) {
  if (n == 1) {
    par(mfrow = c(1, 1))
  } else if (n <= 4) {
    par(mfrow = c(2, 2))
  } else if (n <= 6) {
    par(mfrow = c(2, 3))
  } else if (n <= 9) {
    par(mfrow = c(3, 3))
  } else {
    r <- ceiling(sqrt(n))
    par(mfrow = c(r, ceiling(n / r)))
  }
}

#' Run centrality analysis and populate results
#' @param model The TNA model
#' @param options The jamovi options
#' @param results The jamovi results object
cna_run_centrality <- function(model, options, results) {
  if (is.null(model)) return()
  if (!options$centrality_show_table && !options$centrality_show_plot) return()

  centrality_loops <- options$centrality_loops
  centrality_normalize <- options$centrality_normalize

  vectorCharacter <- character(0)
  if (options$centrality_OutStrength)   vectorCharacter <- c(vectorCharacter, "OutStrength")
  if (options$centrality_InStrength)    vectorCharacter <- c(vectorCharacter, "InStrength")
  if (options$centrality_ClosenessIn)   vectorCharacter <- c(vectorCharacter, "ClosenessIn")
  if (options$centrality_ClosenessOut)  vectorCharacter <- c(vectorCharacter, "ClosenessOut")
  if (options$centrality_Closeness)     vectorCharacter <- c(vectorCharacter, "Closeness")
  if (options$centrality_Betweenness)   vectorCharacter <- c(vectorCharacter, "Betweenness")
  if (options$centrality_BetweennessRSP) vectorCharacter <- c(vectorCharacter, "BetweennessRSP")
  if (options$centrality_Diffusion)     vectorCharacter <- c(vectorCharacter, "Diffusion")
  if (options$centrality_Clustering)    vectorCharacter <- c(vectorCharacter, "Clustering")

  cent <- results$centralityTable$state

  if (length(vectorCharacter) > 0 && is.null(cent)) {
    tryCatch({
      cent <- tna::centralities(
        x = model, loops = centrality_loops,
        normalize = centrality_normalize, measures = vectorCharacter
      )
      results$centralityTable$setState(cent)
    }, error = function(e) {
      results$centralityTable$setNote(
        key = "error", note = paste("Error:", e$message)
      )
    })
  }

  # Add columns
  for (measure in vectorCharacter) {
    col_type <- if (measure == "Betweenness") "integer" else "number"
    results$centralityTable$addColumn(name = measure, type = col_type)
  }

  # Populate table
  if (!is.null(cent) && is.data.frame(cent) && nrow(cent) > 0) {
    for (i in 1:nrow(cent)) {
      rowValues <- list()
      if ("group" %in% colnames(cent)) {
        rowValues$group <- as.character(cent[i, "group"])
        rowValues$state <- as.character(cent[i, "state"])
      } else {
        rowValues$group <- ""
        rowValues$state <- as.character(cent[i, "state"])
      }
      for (measure in vectorCharacter) {
        if (measure %in% colnames(cent)) {
          rowValues[[measure]] <- as.numeric(cent[i, measure])
        }
      }
      results$centralityTable$addRow(rowKey = i, values = rowValues)
    }
  }

  results$centralityTitle$setVisible(
    options$centrality_show_table || options$centrality_show_plot
  )
  results$centrality_plot$setVisible(options$centrality_show_plot)
  results$centralityTable$setVisible(options$centrality_show_table)
}

#' Run community detection and populate results
#' @param model The TNA model
#' @param options The jamovi options
#' @param results The jamovi results object
#' @param is_multi Whether this is a multi-network (group/cluster) model
cna_run_community <- function(model, options, results, is_multi) {
  if (is.null(model)) return()
  if (!isTRUE(options$community_show_plot) &&
      !isTRUE(options$community_show_table)) return()

  community_gamma <- as.numeric(options$community_gamma)
  methods <- options$community_methods

  coms <- results$community_plot$state
  if (is.null(coms)) {
    tryCatch({
      coms <- tna::communities(
        x = model, methods = methods, gamma = community_gamma
      )
      results$community_plot$setState(coms)
    }, error = function(e) {
      results$communityErrorText$setContent(
        paste("Community detection error:", e$message)
      )
      results$communityErrorText$setVisible(TRUE)
    })
  }

  # Populate communities table
  if (!is.null(coms) && isTRUE(options$community_show_table)) {
    if (is_multi) {
      row_key <- 1
      method_names <- NULL
      for (group_name in names(coms)) {
        group_coms <- coms[[group_name]]
        if (!is.null(group_coms) && !is.null(group_coms$assignments)) {
          assignments <- group_coms$assignments
          if (is.null(method_names)) {
            method_names <- colnames(assignments)[colnames(assignments) != "state"]
            for (method in method_names) {
              results$communityTable$addColumn(
                name = method, title = method, type = "integer"
              )
            }
          }
          for (i in 1:nrow(assignments)) {
            rowValues <- list(
              group = as.character(group_name),
              state = as.character(assignments[i, "state"])
            )
            for (method in method_names) {
              rowValues[[method]] <- as.integer(assignments[i, method])
            }
            results$communityTable$addRow(rowKey = row_key, values = rowValues)
            row_key <- row_key + 1
          }
        }
      }
    } else {
      if (!is.null(coms$assignments)) {
        assignments <- coms$assignments
        method_names <- colnames(assignments)[-1]
        for (method in method_names) {
          results$communityTable$addColumn(
            name = method, title = method, type = "integer"
          )
        }
        for (i in 1:nrow(assignments)) {
          rowValues <- list(
            group = "",
            state = as.character(assignments[i, "state"])
          )
          for (method in method_names) {
            rowValues[[method]] <- as.integer(assignments[i, method])
          }
          results$communityTable$addRow(rowKey = i, values = rowValues)
        }
      }
    }
  }

  results$community_plot$setVisible(options$community_show_plot)
  results$communityTable$setVisible(options$community_show_table)
  results$communityContent$setVisible(FALSE)
  results$communityTitle$setVisible(
    isTRUE(options$community_show_plot) || isTRUE(options$community_show_table)
  )
}

#' Run cliques analysis and populate results
#' @param model The TNA model
#' @param options The jamovi options
#' @param results The jamovi results object
cna_run_cliques <- function(model, options, results) {
  if (is.null(model)) return()
  if (!isTRUE(options$cliques_show_text) &&
      !isTRUE(options$cliques_show_plot)) return()

  cliques_size <- as.numeric(options$cliques_size)
  cliques_threshold <- as.numeric(options$cliques_threshold)
  if (cliques_threshold == 0) {
    cliques_threshold <- if (cliques_size <= 2) 0.1 else 0.01
  }

  cliques <- results$cliques_multiple_plot$state
  if (is.null(cliques)) {
    cliques <- tna::cliques(
      x = model, size = cliques_size, threshold = cliques_threshold
    )
    results$cliques_multiple_plot$setState(cliques)
    if (isTRUE(options$cliques_show_text)) {
      results$cliquesContent$setContent(cliques)
    }
  }

  results$cliques_multiple_plot$setVisible(options$cliques_show_plot)
  results$cliquesContent$setVisible(options$cliques_show_text)
  results$cliquesTitle$setVisible(
    isTRUE(options$cliques_show_text) || isTRUE(options$cliques_show_plot)
  )
}

#' Run bootstrap analysis and populate results
#' @param model The TNA model
#' @param options The jamovi options
#' @param results The jamovi results object
#' @param is_multi Whether this is a multi-network model
cna_run_bootstrap <- function(model, options, results, is_multi) {
  if (is.null(model)) return()
  if (!isTRUE(options$bootstrap_show_table) &&
      !isTRUE(options$bootstrap_show_plot)) return()

  bs <- results$bootstrap_plot$state
  if (is.null(bs)) {
    bs <- tna::bootstrap(
      x = model,
      iter = options$bootstrap_iteration,
      level = options$bootstrap_level,
      method = options$bootstrap_method,
      threshold = options$bootstrap_threshold,
      consistency_range = c(options$bootstrap_range_low, options$bootstrap_range_up)
    )
    results$bootstrap_plot$setState(bs)
  }

  # Populate bootstrap table
  if (!is.null(bs) && isTRUE(options$bootstrap_show_table)) {
    row_key <- 1
    max_rows <- options$bootstrap_table_max_rows
    show_all <- isTRUE(options$bootstrap_table_show_all)
    significant_only <- isTRUE(options$bootstrap_table_significant_only)

    if (is_multi) {
      tryCatch({
        for (group_name in names(bs)) {
          group_data <- bs[[group_name]]
          if (!is.null(group_data) && !is.null(group_data$summary)) {
            summary_data <- group_data$summary
            if (is.data.frame(summary_data) && nrow(summary_data) > 0) {
              summary_data <- summary_data[
                order(-summary_data$sig, summary_data$p_value),
              ]
              if (significant_only) {
                summary_data <- summary_data[summary_data$sig == TRUE, ]
              }
              if (nrow(summary_data) == 0) next
              for (i in 1:nrow(summary_data)) {
                if (!show_all && row_key > max_rows) break
                results$bootstrapTable$addRow(
                  rowKey = row_key,
                  values = list(
                    group = group_name,
                    from = as.character(summary_data[i, "from"]),
                    to = as.character(summary_data[i, "to"]),
                    weight = as.numeric(summary_data[i, "weight"]),
                    p_value = as.numeric(summary_data[i, "p_value"]),
                    cr_lower = as.numeric(summary_data[i, "cr_lower"]),
                    cr_upper = as.numeric(summary_data[i, "cr_upper"]),
                    ci_lower = as.numeric(summary_data[i, "ci_lower"]),
                    ci_upper = as.numeric(summary_data[i, "ci_upper"]),
                    significant = ifelse(summary_data[i, "sig"], "Yes", "No")
                  )
                )
                row_key <- row_key + 1
              }
              if (!show_all && row_key > max_rows) break
            }
          }
        }
      }, error = function(e) {
        results$bootstrapTable$setNote(
          key = "error", note = paste("Bootstrap table error:", e$message)
        )
      })
    } else {
      if (!is.null(bs$summary) && nrow(bs$summary) > 0) {
        all_edges <- bs$summary
        all_edges <- all_edges[order(-all_edges$sig, all_edges$p_value), ]
        if (significant_only) {
          all_edges <- all_edges[all_edges$sig == TRUE, ]
        }
        if (!show_all && nrow(all_edges) > max_rows) {
          all_edges <- all_edges[1:max_rows, ]
        }
        for (i in 1:nrow(all_edges)) {
          results$bootstrapTable$addRow(
            rowKey = i,
            values = list(
              group = "",
              from = as.character(all_edges[i, "from"]),
              to = as.character(all_edges[i, "to"]),
              weight = all_edges[i, "weight"],
              p_value = all_edges[i, "p_value"],
              cr_lower = all_edges[i, "cr_lower"],
              cr_upper = all_edges[i, "cr_upper"],
              ci_lower = all_edges[i, "ci_lower"],
              ci_upper = all_edges[i, "ci_upper"],
              significant = ifelse(all_edges[i, "sig"], "Yes", "No")
            )
          )
        }
      }
    }
  }

  results$bootstrap_plot$setVisible(options$bootstrap_show_plot)
  results$bootstrapTable$setVisible(options$bootstrap_show_table)
  results$bootstrapTitle$setVisible(
    isTRUE(options$bootstrap_show_plot) || isTRUE(options$bootstrap_show_table)
  )
}

#' Run permutation test and populate results (Group/Cluster only)
#' @param model The TNA model
#' @param options The jamovi options
#' @param results The jamovi results object
cna_run_permutation <- function(model, options, results) {
  if (is.null(model)) return()
  if (!isTRUE(options$permutation_show_table) &&
      !isTRUE(options$permutation_show_plot)) return()

  permTest <- results$permutation_plot$state
  if (is.null(permTest)) {
    tryCatch({
      permTest <- tna::permutation_test(
        x = model,
        iter = options$permutation_iter,
        paired = options$permutation_paired,
        level = options$permutation_level
      )
      results$permutation_plot$setState(permTest)
    }, error = function(e) {
      results$permutationTitle$setVisible(TRUE)
      results$permutationTable$setNote(
        key = "error", note = paste("Error:", e$message)
      )
    })
  }

  # Populate permutation table
  if (!is.null(permTest) && isTRUE(options$permutation_show_table)) {
    row_key <- 1
    max_rows <- options$permutation_table_max_rows
    show_all <- isTRUE(options$permutation_table_show_all)
    for (comp_name in names(permTest)) {
      comp_data <- permTest[[comp_name]]
      if (!is.null(comp_data$edges) && !is.null(comp_data$edges$stats)) {
        stats_df <- comp_data$edges$stats
        stats_df$diff_true <- as.numeric(stats_df$diff_true)
        stats_df$effect_size <- as.numeric(stats_df$effect_size)
        stats_df$p_value <- as.numeric(stats_df$p_value)
        stats_df <- stats_df[
          order(stats_df$p_value, -abs(stats_df$diff_true)),
        ]

        for (i in 1:nrow(stats_df)) {
          if (!show_all && row_key > max_rows) break
          edge_name_value <- stats_df[i, "edge_name"]
          if (is.null(edge_name_value) || is.na(edge_name_value) ||
              edge_name_value == "" || edge_name_value == " -> ") {
            edge_name_value <- rownames(stats_df)[i]
            if (is.null(edge_name_value) || edge_name_value == "") {
              edge_name_value <- paste("Edge", i)
            }
          } else {
            edge_name_value <- trimws(as.character(edge_name_value))
          }

          results$permutationTable$addRow(
            rowKey = row_key,
            values = list(
              group_comparison = comp_name,
              edge_name = edge_name_value,
              diff_true = as.numeric(stats_df[i, "diff_true"]),
              effect_size = as.numeric(stats_df[i, "effect_size"]),
              p_value = as.numeric(stats_df[i, "p_value"])
            )
          )
          row_key <- row_key + 1
        }
        if (!show_all && row_key > max_rows) break
      }
    }
  }

  results$permutation_plot$setVisible(options$permutation_show_plot)
  results$permutationTable$setVisible(options$permutation_show_table)
  results$permutationTitle$setVisible(
    isTRUE(options$permutation_show_table) ||
      isTRUE(options$permutation_show_plot)
  )
}

#' Run compare analysis and populate results (Group/Cluster only)
#' @param model The TNA model
#' @param options The jamovi options
#' @param results The jamovi results object
#' @param mode "group" for LevelSelector, "cluster" for integer selectors
cna_run_compare <- function(model, options, results, mode) {
  if (is.null(model)) return()

  showAnyCompare <- isTRUE(options$compare_show_summary) ||
    isTRUE(options$compare_show_network) ||
    isTRUE(options$compare_show_plot)

  if (!showAnyCompare) return()

  results$compareInstructions$setContent(
    '<div style="border: 2px solid #d4edda; border-radius: 10px; padding: 10px; background-color: #d4edda; margin: 10px 0;">
    <b>Compare Network Properties</b>: Compares general network properties between two groups, including edge weight correlations, distances, and structural metrics.
    </div>'
  )
  results$compareInstructions$setVisible(TRUE)

  available_groups <- names(model)
  num_groups <- length(available_groups)
  group_i_name <- NULL
  group_j_name <- NULL

  if (mode == "group") {
    group_i_name <- options$compare_group_i
    group_j_name <- options$compare_group_j
    if (is.null(group_i_name) || group_i_name == "") {
      group_i_name <- available_groups[1]
    }
    if (is.null(group_j_name) || group_j_name == "") {
      if (num_groups >= 2) {
        group_j_name <- available_groups[2]
      } else {
        group_j_name <- available_groups[1]
      }
    }
  } else {
    # clustering: use integer selectors
    ci <- options$compare_cluster_i
    cj <- options$compare_cluster_j
    if (ci > num_groups || cj > num_groups) {
      results$compareTitle$setContent(
        paste("Invalid cluster index. Available: 1 to", num_groups)
      )
      results$compareTitle$setVisible(TRUE)
      return()
    } else if (ci == cj) {
      results$compareTitle$setContent(
        "Please select two different clusters to compare"
      )
      results$compareTitle$setVisible(TRUE)
      return()
    } else {
      group_i_name <- available_groups[ci]
      group_j_name <- available_groups[cj]
    }
  }

  if (!is.null(group_i_name) && !is.null(group_j_name) &&
      group_i_name != group_j_name) {

    compResult <- results$compare_plot$state
    if (is.null(compResult)) {
      tryCatch({
        compResult <- tna::compare(
          x = model,
          i = group_i_name,
          j = group_j_name,
          scaling = options$compare_scaling,
          network = TRUE
        )
        results$compare_plot$setState(list(
          result = compResult,
          group_i = group_i_name,
          group_j = group_j_name
        ))
      }, error = function(e) {
        results$errorText$setContent(paste("Compare Error:", e$message))
        results$errorText$setVisible(TRUE)
      })
    } else {
      group_i_name <- compResult$group_i
      group_j_name <- compResult$group_j
      compResult <- compResult$result
    }

    results$compareTitle$setContent(
      paste("Comparing:", group_i_name, "vs", group_j_name)
    )

    # Summary metrics table
    if (!is.null(compResult) && isTRUE(options$compare_show_summary)) {
      if (!is.null(compResult$summary_metrics) &&
          nrow(compResult$summary_metrics) > 0) {
        results$compareSummaryTable$setTitle(
          paste("Summary Metrics:", group_i_name, "vs", group_j_name)
        )
        for (i in 1:nrow(compResult$summary_metrics)) {
          results$compareSummaryTable$addRow(
            rowKey = i,
            values = list(
              metric = as.character(compResult$summary_metrics$metric[i]),
              value = as.numeric(compResult$summary_metrics$value[i])
            )
          )
        }
      }
      results$compareSummaryTable$setVisible(TRUE)
    }

    # Network properties table
    if (!is.null(compResult) && isTRUE(options$compare_show_network)) {
      if (!is.null(compResult$network_metrics) &&
          nrow(compResult$network_metrics) > 0) {
        results$compareNetworkTable$setTitle(
          paste("Network Properties:", group_i_name, "vs", group_j_name)
        )
        results$compareNetworkTable$getColumn("group_i")$setTitle(group_i_name)
        results$compareNetworkTable$getColumn("group_j")$setTitle(group_j_name)
        for (i in 1:nrow(compResult$network_metrics)) {
          row_data <- compResult$network_metrics[i, ]
          results$compareNetworkTable$addRow(
            rowKey = i,
            values = list(
              metric = as.character(row_data$metric),
              group_i = as.numeric(row_data[[2]]),
              group_j = as.numeric(row_data[[3]])
            )
          )
        }
      }
      results$compareNetworkTable$setVisible(TRUE)
    }

    results$compare_plot$setVisible(options$compare_show_plot)
    results$compareTitle$setVisible(showAnyCompare)
  }
}

# ── Shared plot functions ──

#' Show build model plot
#' @param plotData The cached model state
#' @param options The jamovi options
#' @param results The jamovi results object
#' @param is_multi Whether this is a multi-network model
cna_plot_model <- function(plotData, options, results, is_multi) {
  if (is.null(plotData)) return(FALSE)

  if (is_multi) {
    cna_set_panel_layout(length(plotData))
  }

  tryCatch({
    plot(
      x = plotData,
      cut = options$buildModel_plot_cut,
      minimum = options$buildModel_plot_min_value,
      edge.label.cex = options$buildModel_plot_edge_label_size,
      node.width = options$buildModel_plot_node_size,
      label.cex = options$buildModel_plot_node_label_size,
      bg = "transparent"
    )
  }, error = function(e) {
    results$errorText$setContent(paste0("Plot error: ", e$message))
    results$errorText$setVisible(TRUE)
  })
  TRUE
}

#' Show histogram plot
#' @param plotData The cached model state
#' @param options The jamovi options
#' @param is_multi Whether this is a multi-network model
cna_plot_histo <- function(plotData, options, is_multi) {
  if (is.null(plotData) || !options$buildModel_show_histo) return(FALSE)

  if (is_multi) {
    n <- length(plotData)
    cna_set_panel_layout(n)
    tryCatch({
      for (i in 1:n) {
        group_name <- names(plotData)[i]
        if (is.null(group_name)) group_name <- paste("Group", i)
        w <- c(plotData[[i]]$weights)
        brks <- seq(0, max(1, max(w, na.rm = TRUE)) + 0.01, length.out = 20)
        hist(
          x = plotData[[i]], breaks = brks,
          main = paste("Histogram -", group_name),
          xlab = "Edge Weights", ylab = "Frequency"
        )
      }
    }, error = function(e) {
      plot(1, type = "n", main = "Histogram Error", sub = e$message)
    })
  } else {
    hist(
      x = plotData, main = "Histogram of Edge Weights",
      xlab = "Edge Weights", ylab = "Frequency"
    )
  }
  TRUE
}

#' Show frequencies plot
cna_plot_frequencies <- function(plotData) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    p <- tna::plot_frequencies(x = plotData)
    if (!is.null(p)) print(p)
  }, error = function(e) {
    hist(plotData, main = "Frequencies Plot",
         xlab = "Edge Weights", ylab = "Frequency")
  })
  TRUE
}

#' Show mosaic plot
cna_plot_mosaic <- function(plotData, options) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    p <- tna::plot_mosaic(x = plotData, digits = options$buildModel_digits)
    print(p)
  }, error = function(e) {
    plot(1, type = "n", main = "Mosaic Plot Error", sub = e$message)
  })
  TRUE
}

#' Show centrality plot
cna_plot_centrality <- function(plotData) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    print(plot(plotData))
  }, error = function(e) {
    plot(1, type = "n", main = "Centrality Plot Error", sub = e$message)
  })
  TRUE
}

#' Show community plot
cna_plot_community <- function(plotData, options, is_multi) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    if (is_multi) {
      cna_set_panel_layout(length(plotData))
    }
    plot(x = plotData, method = options$community_methods, bg = "transparent")
  }, error = function(e) {
    plot(1, type = "n", main = "Community Plot Error", sub = e$message)
  })
  TRUE
}

#' Show cliques plot
cna_plot_cliques <- function(plotData, options) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    n <- lengths(plotData[1])
    if (n == 0) return(FALSE)
    nc <- ceiling(sqrt(n))
    nr <- ceiling(n / nc)
    par(mfrow = c(nr, nc))
    for (i in seq_len(n)) {
      plot(
        x = plotData, ask = FALSE, first = i, n = 1,
        cut = options$cliques_plot_cut,
        minimum = options$cliques_plot_min_value,
        edge.label.cex = options$cliques_plot_edge_label_size,
        node.width = options$cliques_plot_node_size,
        label.cex = options$cliques_plot_node_label_size
      )
    }
  }, error = function(e) {
    plot(1, type = "n", main = "Cliques Plot Error", sub = e$message)
  })
  TRUE
}

#' Show bootstrap plot
cna_plot_bootstrap <- function(plotData, is_multi) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    if (is_multi) {
      cna_set_panel_layout(length(plotData))
    }
    plot(x = plotData, cut = 0.1)
  }, error = function(e) {
    plot(1, type = "n", main = "Bootstrap Plot Error", sub = e$message)
  })
  TRUE
}

#' Show permutation plot
cna_plot_permutation <- function(plotData) {
  if (is.null(plotData)) return(FALSE)
  tryCatch({
    cna_set_panel_layout(length(plotData))
    plot(x = plotData)
  }, error = function(e) {
    plot(1, type = "n", main = "Permutation Plot Error", sub = e$message)
  })
  TRUE
}

#' Show compare plot
cna_plot_compare <- function(stateData, options) {
  if (is.null(stateData)) return(FALSE)
  tryCatch({
    plotData <- stateData$result
    name_x <- stateData$group_i
    name_y <- stateData$group_j
    if (is.null(plotData)) return(FALSE)
    p <- plot(
      x = plotData, type = options$compare_plot_type,
      name_x = name_x, name_y = name_y
    )
    print(p)
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Compare Plot Error", sub = e$message)
  })
  TRUE
}

#' Show network difference plot
cna_plot_network_diff <- function(model, options) {
  if (is.null(model)) return(FALSE)

  n_groups <- length(model)
  n_comparisons <- choose(n_groups, 2)

  if (n_comparisons == 1) {
    par(mfrow = c(1, 1))
  } else if (n_comparisons <= 4) {
    par(mfrow = c(2, 2))
  } else if (n_comparisons <= 6) {
    par(mfrow = c(2, 3))
  } else if (n_comparisons <= 9) {
    par(mfrow = c(3, 3))
  } else {
    rows <- ceiling(sqrt(n_comparisons))
    cols <- ceiling(n_comparisons / rows)
    par(mfrow = c(rows, cols))
  }

  tryCatch({
    group_names <- names(model)
    for (i in 1:(n_groups - 1)) {
      for (j in (i + 1):n_groups) {
        plot_title <- paste(group_names[i], "vs", group_names[j])
        tna::plot_compare(
          x = model, i = i, j = j,
          cut = options$compare_network_diff_plot_cut,
          minimum = options$compare_network_diff_plot_min_value,
          edge.label.cex = options$compare_network_diff_plot_edge_label_size,
          node.width = options$compare_network_diff_plot_node_size,
          label.cex = options$compare_network_diff_plot_node_label_size,
          title = plot_title
        )
      }
    }
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Network Diff Error", sub = e$message)
  })
  TRUE
}

# ============================================================================
#  NESTIMATE ENGINE  (CNA v2)
# ----------------------------------------------------------------------------
#  Co-occurrence networks are built and analysed with the Nestimate package
#  and plotted automatically with cograph::splot. Group networks use
#  Nestimate's native `group=` argument (no manual looping), and splot renders
#  a netobject_group / net_bootstrap_group / net_permutation_group as a grid
#  on its own. Function names carry a `_nest` suffix so the legacy tna helpers
#  above keep working until every analysis is migrated.
#  Verified against Nestimate 0.7.0 / cograph 2.3.9 (see LEARNINGS.md).
# ============================================================================

# Methods cograph::detect_communities actually accepts.
.cna_supported_communities <- c("louvain", "walktrap", "fast_greedy",
                                "label_prop", "infomap", "leiden")

# Pretty titles for cograph::centrality column names.
.cna_centrality_titles <- c(
  strength_all = "Strength", strength = "Strength",
  betweenness = "Betweenness",
  closeness_all = "Closeness", closeness = "Closeness",
  degree_all = "Degree", degree = "Degree",
  eigenvector = "Eigenvector", pagerank = "PageRank"
)

#' Build an event (long-format) co-occurrence network with Nestimate.
#'
#' @return A `netobject` (single mode) or `netobject_group` (group mode,
#'   built via Nestimate's `group=` argument — no manual per-group loop).
cna_build_event_nest <- function(data, options) {
  action <- options$buildModel_variables_long_action
  actor  <- options$buildModel_variables_actor
  time   <- options$buildModel_variables_long_time
  order  <- options$buildModel_variables_long_order
  group  <- options$buildModel_variables_group
  is_multi <- !is.null(group)

  scaling <- options$buildModel_scaling
  scaling <- if (is.null(scaling) || scaling == "noScaling") NULL else scaling

  keep <- Filter(Negate(is.null),
                 list(time, actor, action, order, if (is_multi) group))
  df <- as.data.frame(data)[unlist(keep)]

  df[[action]] <- as.character(df[[action]])
  if (!is.null(actor)) df[[actor]] <- as.character(df[[actor]])
  if (!is.null(time))  df[[time]]  <- as.POSIXct(df[[time]])
  if (!is.null(order)) {
    ord <- suppressWarnings(as.numeric(as.character(df[[order]])))
    if (!anyNA(ord)) df[[order]] <- ord   # keep numeric order (LEARNINGS)
  }

  args <- Filter(Negate(is.null), list(
    data = df, method = "co_occurrence",
    action = action, actor = actor, time = time, order = order,
    time_threshold = options$buildModel_threshold,
    scaling = scaling,
    group = if (is_multi) group else NULL
  ))
  do.call(Nestimate::build_network, args)
}

#' Map the module's centrality checkboxes onto cograph::centrality measure
#' names. Co-occurrence networks are undirected, so In/Out-strength collapse to
#' "strength" and In/Out-closeness collapse to "closeness".
cna_cograph_measures <- function(options) {
  unique(c(
    if (isTRUE(options$centrality_InStrength) ||
        isTRUE(options$centrality_OutStrength)) "strength",
    if (isTRUE(options$centrality_Betweenness)) "betweenness",
    if (isTRUE(options$centrality_Closeness)) "closeness"
  ))
}

#' Compute centralities with cograph::centrality, returning one tidy
#' data.frame: group, state, <cograph measure columns...>. Group mode maps
#' over the netobject_group (cograph::centrality takes a single network).
cna_centrality_compute_nest <- function(model, options, is_multi) {
  measures <- cna_cograph_measures(options)
  if (!length(measures)) measures <- "betweenness"
  one <- function(net, g) {
    ct <- cograph::centrality(net, measures = measures)
    meas <- setdiff(names(ct), "node")
    data.frame(group = g, state = ct$node, ct[meas],
               check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  }
  if (is_multi) do.call(rbind, Map(one, model, names(model)))
  else one(model, "")
}

#' Centrality table (single or group) computed with cograph::centrality.
cna_run_centrality_nest <- function(model, options, results, is_multi) {
  if (is.null(model)) return()
  if (!options$centrality_show_table && !options$centrality_show_plot) return()

  cent <- results$centralityTable$state
  if (is.null(cent)) {
    cent <- tryCatch(
      cna_centrality_compute_nest(model, options, is_multi),
      error = function(e) {
        results$centralityTable$setNote("error", paste("Error:", e$message))
        NULL
      }
    )
    if (!is.null(cent)) results$centralityTable$setState(cent)
  }

  if (!is.null(cent)) {
    meas_cols <- setdiff(names(cent), c("group", "state"))
    invisible(lapply(meas_cols, function(c0) {
      ttl <- if (c0 %in% names(.cna_centrality_titles))
        .cna_centrality_titles[[c0]] else c0
      results$centralityTable$addColumn(name = c0, title = ttl, type = "number")
    }))
    if (isTRUE(options$centrality_show_table) && nrow(cent) > 0) {
      invisible(lapply(seq_len(nrow(cent)), function(i)
        results$centralityTable$addRow(rowKey = i, values = c(
          list(group = cent$group[i], state = cent$state[i]),
          setNames(lapply(meas_cols, function(c0) as.numeric(cent[[c0]][i])),
                   meas_cols)))))
    }
  }

  results$centralityTitle$setVisible(
    options$centrality_show_table || options$centrality_show_plot)
  results$centrality_plot$setVisible(options$centrality_show_plot)
  results$centralityTable$setVisible(options$centrality_show_table)
}

#' Community detection (single or group) via cograph::detect_communities.
cna_run_community_nest <- function(model, options, results, is_multi) {
  if (is.null(model)) return()
  if (!isTRUE(options$community_show_plot) &&
      !isTRUE(options$community_show_table)) return()

  method <- options$community_methods
  if (!(method %in% .cna_supported_communities)) method <- "louvain"

  coms <- results$community_plot$state
  if (is.null(coms)) {
    coms <- tryCatch(
      if (is_multi) lapply(model, cograph::detect_communities, method = method)
      else cograph::detect_communities(model, method = method),
      error = function(e) {
        results$communityErrorText$setContent(
          paste("Community detection error:", e$message))
        results$communityErrorText$setVisible(TRUE)
        NULL
      }
    )
    if (!is.null(coms)) results$community_plot$setState(coms)
  }

  if (!is.null(coms) && isTRUE(options$community_show_table)) {
    results$communityTable$addColumn(name = method, title = method,
                                     type = "integer")
    rows <- if (is_multi)
      do.call(rbind, Map(function(cdf, g)
        data.frame(group = g, state = cdf$node, comm = cdf$community,
                   stringsAsFactors = FALSE), coms, names(coms)))
    else
      data.frame(group = "", state = coms$node, comm = coms$community,
                 stringsAsFactors = FALSE)
    invisible(lapply(seq_len(nrow(rows)), function(i)
      results$communityTable$addRow(rowKey = i, values = setNames(
        list(rows$group[i], rows$state[i], as.integer(rows$comm[i])),
        c("group", "state", method)))))
  }

  results$community_plot$setVisible(options$community_show_plot)
  results$communityTable$setVisible(options$community_show_table)
  results$communityContent$setVisible(FALSE)
  results$communityTitle$setVisible(
    isTRUE(options$community_show_plot) || isTRUE(options$community_show_table))
}

#' Bootstrap (single or group) via Nestimate::bootstrap_network.
cna_run_bootstrap_nest <- function(model, options, results, is_multi) {
  if (is.null(model)) return()
  if (!isTRUE(options$bootstrap_show_table) &&
      !isTRUE(options$bootstrap_show_plot)) return()

  bs <- results$bootstrap_plot$state
  if (is.null(bs)) {
    bs <- tryCatch(Nestimate::bootstrap_network(
      model,
      iter = options$bootstrap_iteration,
      ci_level = options$bootstrap_level,
      inference = options$bootstrap_method,
      consistency_range = c(options$bootstrap_range_low,
                            options$bootstrap_range_up),
      edge_threshold = if (identical(options$bootstrap_method, "threshold"))
        options$bootstrap_threshold else NULL
    ), error = function(e) {
      results$bootstrapTable$setNote("error", paste("Error:", e$message)); NULL
    })
    if (!is.null(bs)) results$bootstrap_plot$setState(bs)
  }

  if (!is.null(bs) && isTRUE(options$bootstrap_show_table)) {
    sig_only <- isTRUE(options$bootstrap_table_significant_only)
    show_all <- isTRUE(options$bootstrap_table_show_all)
    max_rows <- options$bootstrap_table_max_rows

    # one tidy summary, group column included (loop-free across groups)
    summ <- if (is_multi)
      do.call(rbind, Map(function(b, g) {
        s <- b$summary; if (is.null(s) || !nrow(s)) return(NULL)
        s$group <- g; s
      }, bs, names(bs)))
    else { s <- bs$summary; if (!is.null(s) && nrow(s)) s$group <- ""; s }

    if (!is.null(summ) && nrow(summ) > 0) {
      summ <- summ[order(-summ$sig, summ$p_value), ]
      if (sig_only) summ <- summ[summ$sig %in% TRUE, ]
      if (!show_all && nrow(summ) > max_rows) summ <- summ[seq_len(max_rows), ]
      invisible(lapply(seq_len(nrow(summ)), function(i)
        results$bootstrapTable$addRow(rowKey = i, values = list(
          group = as.character(summ$group[i]),
          from = as.character(summ$from[i]), to = as.character(summ$to[i]),
          weight = as.numeric(summ$weight[i]),
          p_value = as.numeric(summ$p_value[i]),
          cr_lower = as.numeric(summ$cr_lower[i]),
          cr_upper = as.numeric(summ$cr_upper[i]),
          ci_lower = as.numeric(summ$ci_lower[i]),
          ci_upper = as.numeric(summ$ci_upper[i]),
          significant = ifelse(isTRUE(summ$sig[i]), "Yes", "No")))))
    }
  }

  results$bootstrap_plot$setVisible(options$bootstrap_show_plot)
  results$bootstrapTable$setVisible(options$bootstrap_show_table)
  results$bootstrapTitle$setVisible(
    isTRUE(options$bootstrap_show_plot) || isTRUE(options$bootstrap_show_table))
}

#' Permutation test (group only) via Nestimate::permutation.
cna_run_permutation_nest <- function(model, options, results) {
  if (is.null(model)) return()
  if (!isTRUE(options$permutation_show_table) &&
      !isTRUE(options$permutation_show_plot)) return()

  pt <- results$permutation_plot$state
  if (is.null(pt)) {
    pt <- tryCatch(Nestimate::permutation(
      model, iter = options$permutation_iter,
      alpha = options$permutation_level,
      paired = isTRUE(options$permutation_paired)
    ), error = function(e) {
      results$permutationTitle$setVisible(TRUE)
      results$permutationTable$setNote("error", paste("Error:", e$message)); NULL
    })
    if (!is.null(pt)) results$permutation_plot$setState(pt)
  }

  if (!is.null(pt) && isTRUE(options$permutation_show_table)) {
    max_rows <- options$permutation_table_max_rows
    show_all <- isTRUE(options$permutation_table_show_all)
    rows <- do.call(rbind, Map(function(comp, nm) {
      s <- comp$summary; if (is.null(s) || !nrow(s)) return(NULL)
      data.frame(comparison = nm,
                 edge = paste(s$from, "->", s$to),
                 diff = s$diff, effect_size = s$effect_size,
                 p_value = s$p_value, stringsAsFactors = FALSE)
    }, pt, names(pt)))
    if (!is.null(rows) && nrow(rows) > 0) {
      rows <- rows[order(rows$p_value, -abs(rows$diff)), ]
      if (!show_all && nrow(rows) > max_rows) rows <- rows[seq_len(max_rows), ]
      invisible(lapply(seq_len(nrow(rows)), function(i)
        results$permutationTable$addRow(rowKey = i, values = list(
          group_comparison = rows$comparison[i],
          edge_name = rows$edge[i],
          diff_true = as.numeric(rows$diff[i]),
          effect_size = as.numeric(rows$effect_size[i]),
          p_value = as.numeric(rows$p_value[i])))))
    }
  }

  results$permutation_plot$setVisible(options$permutation_show_plot)
  results$permutationTable$setVisible(options$permutation_show_table)
  results$permutationTitle$setVisible(
    isTRUE(options$permutation_show_table) ||
      isTRUE(options$permutation_show_plot))
}

#' Network comparison (group only) via Nestimate::compare_model.
cna_run_compare_nest <- function(model, options, results, mode) {
  if (is.null(model)) return()
  show_any <- isTRUE(options$compare_show_summary) ||
    isTRUE(options$compare_show_network) || isTRUE(options$compare_show_plot)
  if (!show_any) return()

  groups <- names(model)
  ng <- length(groups)
  if (mode == "group") {
    gi <- options$compare_group_i; gj <- options$compare_group_j
    if (is.null(gi) || gi == "") gi <- groups[1]
    if (is.null(gj) || gj == "") gj <- groups[min(2, ng)]
    i <- match(gi, groups); j <- match(gj, groups)
  } else {
    i <- options$compare_cluster_i; j <- options$compare_cluster_j
    if (i > ng || j > ng) {
      results$compareTitle$setContent(paste("Invalid index. Available: 1 to", ng))
      results$compareTitle$setVisible(TRUE); return()
    }
    gi <- groups[i]; gj <- groups[j]
  }
  if (is.na(i) || is.na(j) || i == j) {
    results$compareTitle$setContent("Please select two different networks.")
    results$compareTitle$setVisible(TRUE); return()
  }

  state <- results$compare_plot$state
  if (is.null(state)) {
    cmp <- tryCatch(Nestimate::compare_model(
      model, i = i, j = j, scaling = options$compare_scaling),
      error = function(e) {
        results$errorText$setContent(paste("Compare Error:", e$message))
        results$errorText$setVisible(TRUE); NULL
      })
    if (!is.null(cmp)) {
      state <- list(result = cmp, group_i = gi, group_j = gj)
      results$compare_plot$setState(state)
    }
  }
  if (is.null(state)) return()
  cmp <- state$result; gi <- state$group_i; gj <- state$group_j
  results$compareTitle$setContent(paste("Comparing:", gi, "vs", gj))

  if (isTRUE(options$compare_show_summary) && !is.null(cmp$summary_metrics)) {
    sm <- cmp$summary_metrics
    results$compareSummaryTable$setTitle(paste("Summary Metrics:", gi, "vs", gj))
    invisible(lapply(seq_len(nrow(sm)), function(k)
      results$compareSummaryTable$addRow(rowKey = k, values = list(
        metric = as.character(sm$metric[k]), value = as.numeric(sm$value[k])))))
    results$compareSummaryTable$setVisible(TRUE)
  }

  if (isTRUE(options$compare_show_network) && !is.null(cmp$network_metrics)) {
    nm <- cmp$network_metrics
    results$compareNetworkTable$setTitle(paste("Network Properties:", gi, "vs", gj))
    results$compareNetworkTable$getColumn("group_i")$setTitle(gi)
    results$compareNetworkTable$getColumn("group_j")$setTitle(gj)
    invisible(lapply(seq_len(nrow(nm)), function(k)
      results$compareNetworkTable$addRow(rowKey = k, values = list(
        metric = as.character(nm$metric[k]),
        group_i = as.numeric(nm[[2]][k]), group_j = as.numeric(nm[[3]][k])))))
    results$compareNetworkTable$setVisible(TRUE)
  }

  results$compare_plot$setVisible(options$compare_show_plot)
  results$compareTitle$setVisible(show_any)
}

# ── Nestimate plot helpers (cograph::splot does single + group grids) ──

#' Main network plot. splot auto-renders a netobject_group as a grid and uses
#' its own default layout — we never pass a custom layout. Edge-label digits
#' are left to cograph (it suppresses ".00" for the integer co-occurrence
#' counts on its own).
cna_splot_model_nest <- function(model, options) {
  if (is.null(model)) return(FALSE)
  tryCatch(cograph::splot(model, layout = "oval"),
           error = function(e) plot(1, type = "n", main = "Plot error",
                                    sub = e$message))
  TRUE
}

#' Edge-weight histogram (single, or one panel per group).
cna_splot_histo_nest <- function(model, is_multi) {
  if (is.null(model)) return(FALSE)
  draw <- function(net, ttl) {
    w <- as.numeric(net$weights); w <- w[is.finite(w)]
    hist(w, main = ttl, xlab = "Edge weight", ylab = "Frequency",
         col = "grey80", border = "white")
  }
  tryCatch({
    if (is_multi) {
      cna_set_panel_layout(length(model))
      invisible(Map(function(net, g) draw(net, paste("Histogram -", g)),
                    model, names(model)))
    } else draw(model, "Histogram of Edge Weights")
  }, error = function(e) plot(1, type = "n", main = "Histogram error",
                              sub = e$message))
  TRUE
}

#' State-frequency plot via Nestimate::plot_state_frequencies (handles single
#' and group/cluster objects on its own).
cna_splot_freq_nest <- function(model, is_multi = FALSE) {
  if (is.null(model)) return(FALSE)
  tryCatch(
    invisible(utils::capture.output(print(Nestimate::plot_state_frequencies(model)))),
    error = function(e) plot(1, type = "n", main = "Frequencies error",
                             sub = e$message))
  TRUE
}

#' Mosaic plot via Nestimate::mosaic_plot (ggplot; single or group).
cna_splot_mosaic_nest <- function(model) {
  if (is.null(model)) return(FALSE)
  tryCatch(print(Nestimate::mosaic_plot(model)),
           error = function(e) plot(1, type = "n", main = "Mosaic error",
                                    sub = e$message))
  TRUE
}

#' Edge-weight distribution plot via Nestimate::distribution_plot.
cna_splot_distribution_nest <- function(model) {
  if (is.null(model)) return(FALSE)
  tryCatch({
    p <- Nestimate::distribution_plot(model)
    if (inherits(p, c("ggplot", "gg"))) print(p)
  }, error = function(e) plot(1, type = "n", main = "Distribution error",
                              sub = e$message))
  TRUE
}

#' Centrality plot via cograph::plot_centrality (returns a ggplot). Group mode
#' stacks one plot per group with patchwork.
cna_splot_centrality_nest <- function(model, options, is_multi) {
  if (is.null(model)) return(FALSE)
  measures <- cna_cograph_measures(options)
  if (!length(measures)) measures <- "betweenness"
  style <- options$centrality_plot_style
  if (is.null(style) || !nzchar(style)) style <- "lollipop"
  nets <- if (is_multi) model else stats::setNames(list(model), "")
  tryCatch({
    if (style == "distribution") {
      # base-graphics distribution of centrality values (one panel per group)
      if (is_multi) cna_set_panel_layout(length(nets))
      invisible(Map(function(net, g) cograph::plot_centrality_distribution(
        net, main = if (nzchar(g)) g else "Centrality distribution"),
        nets, names(nets)))
    } else {
      # ggplot styles: lollipop / bar / line / dot (plot_centrality) or
      # heatmap (plot_centrality_heatmap); groups stacked with patchwork
      one <- function(net, g) {
        p <- if (style == "heatmap")
          cograph::plot_centrality_heatmap(net, measures = measures)
        else
          cograph::plot_centrality(net, measures = measures, style = style)
        if (nzchar(g)) p <- p + ggplot2::ggtitle(g)
        p
      }
      if (is_multi)
        print(patchwork::wrap_plots(Map(one, nets, names(nets)), ncol = 1))
      else print(one(nets[[1]], ""))
    }
  }, error = function(e) plot(1, type = "n", main = "Centrality plot error",
                              sub = e$message))
  TRUE
}

#' Community plot. splot.cograph_communities colours nodes by community using
#' its own default layout (no custom layout passed). Group = one panel each.
cna_splot_community_nest <- function(coms, options, is_multi) {
  if (is.null(coms)) return(FALSE)
  tryCatch({
    if (is_multi) {
      cna_set_panel_layout(length(coms))
      invisible(lapply(coms, cograph::splot, layout = "oval"))
    } else {
      cograph::splot(coms, layout = "oval")
    }
  }, error = function(e) plot(1, type = "n", main = "Community error",
                              sub = e$message))
  TRUE
}

#' Bootstrap plot. splot renders net_bootstrap / net_bootstrap_group (grid).
#' edge_label_digits = 0: co-occurrence weights are integer counts, so the
#' default 2-decimal labels would show a redundant ".00" (significance stars
#' are kept).
cna_splot_bootstrap_nest <- function(bs) {
  if (is.null(bs)) return(FALSE)
  tryCatch(cograph::splot(bs, layout = "oval", edge_label_digits = 0),
           error = function(e) plot(1, type = "n", main = "Bootstrap error",
                                    sub = e$message))
  TRUE
}

#' Permutation plot. splot renders net_permutation_group as a grid.
#' splot.net_permutation hardcodes "%.2f" for its star labels, so we go through
#' the template path ({est}{stars}) with edge_label_digits = 0 to show the
#' integer diff with stars and no redundant ".00".
cna_splot_permutation_nest <- function(pt) {
  if (is.null(pt)) return(FALSE)
  tryCatch(cograph::splot(pt, layout = "oval",
                          edge_label_template = "{est}{stars}",
                          edge_label_digits = 0),
           error = function(e) plot(1, type = "n", main = "Permutation error",
                                    sub = e$message))
  TRUE
}

#' Comparison plot. plot.net_comparison returns a ggplot, so print() it.
cna_splot_compare_nest <- function(state, options) {
  if (is.null(state) || is.null(state$result)) return(FALSE)
  tryCatch({
    p <- plot(state$result)
    if (inherits(p, "ggplot")) print(p)
  }, error = function(e) plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
                              main = "Compare error", sub = e$message))
  TRUE
}

# ----------------------------------------------------------------------------
#  One-hot MIXED networks: wtna(method = "both") -> wtna_mixed
#  ($transition + $cooccurrence). splot renders the pair side-by-side.
# ----------------------------------------------------------------------------

#' Build a one-hot mixed network (transition + co-occurrence) with Nestimate.
#' @return A `wtna_mixed` object, or a named list of them in group mode.
cna_build_onehot_mixed_nest <- function(data, options, group_col = NULL) {
  codes  <- options$buildModel_variables_onehot
  window <- options$buildModel_window
  if (is.null(window) || window < 1) window <- 3L
  df <- as.data.frame(data)

  build_one <- function(sub)
    Nestimate::wtna(sub[codes], method = "both", codes = codes,
                    window_size = as.integer(window))

  if (is.null(group_col)) return(build_one(df))
  # group mode without a manual loop: split + lapply
  parts <- split(df, factor(as.character(df[[group_col]])))
  lapply(parts, build_one)
}

#' Build a one-hot (binary) co-occurrence network with Nestimate. Single or
#' group (native `group=`). Window / session aggregation is handled by
#' build_network's one-hot path.
cna_build_onehot_nest <- function(data, options) {
  codes   <- options$buildModel_variables_onehot
  actor   <- options$buildModel_variables_actor
  session <- options$buildModel_variables_session
  group   <- options$buildModel_variables_group
  is_multi <- !is.null(group)
  window <- options$buildModel_window
  if (is.null(window) || window < 1) window <- 1L
  mode <- if (identical(options$buildModel_window_type, "sliding"))
    "overlapping" else "non-overlapping"
  scaling <- options$buildModel_scaling
  scaling <- if (is.null(scaling) || scaling == "noScaling") NULL else scaling

  keep <- Filter(Negate(is.null), c(codes, actor, session, if (is_multi) group))
  df <- as.data.frame(data)[unlist(keep)]

  args <- Filter(Negate(is.null), list(
    data = df, method = "co_occurrence", codes = codes, format = "onehot",
    actor = actor, session = session,
    window_size = as.integer(window), mode = mode,
    scaling = scaling, group = if (is_multi) group else NULL
  ))
  do.call(Nestimate::build_network, args)
}

#' Cluster event sequences and build per-cluster co-occurrence networks.
#' Returns a netobject_group keyed by cluster (so every group _nest helper and
#' splot's grid apply unchanged), or NULL until the Run-clustering box is on.
cna_build_event_cluster_nest <- function(data, options) {
  if (!isTRUE(options$clustering_run)) return(NULL)
  action <- options$buildModel_variables_long_action
  actor  <- options$buildModel_variables_actor
  time   <- options$buildModel_variables_long_time
  order  <- options$buildModel_variables_long_order

  keep <- Filter(Negate(is.null), c(actor, time, action, order))
  df <- as.data.frame(data)[unlist(keep)]
  df[[action]] <- as.character(df[[action]])
  if (!is.null(actor)) df[[actor]] <- as.character(df[[actor]])
  if (!is.null(time))  df[[time]]  <- as.POSIXct(df[[time]])
  if (!is.null(order)) {
    ord <- suppressWarnings(as.numeric(as.character(df[[order]])))
    if (!anyNA(ord)) df[[order]] <- ord
  }
  pargs <- Filter(Negate(is.null), list(
    data = df, actor = actor, time = time, action = action,
    order = order, time_threshold = options$buildModel_threshold))
  seqs <- do.call(Nestimate::prepare, pargs)$sequence_data

  Nestimate::cluster_network(
    seqs, k = options$clustering_k, method = "co_occurrence",
    cluster_by = options$clustering_method,
    dissimilarity = options$clustering_dissimilarity)
}

#' Cluster one-hot windows and build per-cluster co-occurrence networks.
cna_build_onehot_cluster_nest <- function(data, options) {
  if (!isTRUE(options$clustering_run)) return(NULL)
  codes   <- options$buildModel_variables_onehot
  actor   <- options$buildModel_variables_actor
  session <- options$buildModel_variables_session
  window <- options$buildModel_window
  if (is.null(window) || window < 1) window <- 1L
  wtype <- if (identical(options$buildModel_window_type, "sliding"))
    "overlapping" else "non-overlapping"

  keep <- Filter(Negate(is.null), c(codes, actor, session))
  df <- as.data.frame(data)[unlist(keep)]
  seqs <- Nestimate::prepare_onehot(
    df, cols = codes, actor = actor, session = session,
    window_size = as.integer(window), window_type = wtype)

  Nestimate::cluster_network(
    seqs, k = options$clustering_k, method = "co_occurrence",
    cluster_by = options$clustering_method,
    dissimilarity = options$clustering_dissimilarity)
}

#' Mixed-network plot (wtna both: transition + co-occurrence side by side).
cna_splot_mixed_nest <- function(mixed) {
  if (is.null(mixed)) return(FALSE)
  tryCatch(cograph::splot(mixed, layout = "oval"),
           error = function(e) plot(1, type = "n", main = "Mixed plot error",
                                    sub = e$message))
  TRUE
}


# ----------------------------------------------------------------------------
#  Reliability & Stability (Nestimate): centrality_stability (CS per measure,
#  plotted), casedrop_reliability (edge-weight CS), network_reliability
#  (split-half). centrality_stability uses Nestimate net_centrality measure
#  names (InStrength/OutStrength/Betweenness/Closeness), not cograph's.
# ----------------------------------------------------------------------------

#' Map the centrality checkboxes onto Nestimate net_centrality measure names.
cna_nestimate_measures <- function(options) {
  unique(c(
    if (isTRUE(options$centrality_InStrength))  "InStrength",
    if (isTRUE(options$centrality_OutStrength)) "OutStrength",
    if (isTRUE(options$centrality_Betweenness)) "Betweenness",
    if (isTRUE(options$centrality_Closeness))   "Closeness"
  ))
}

#' Tidy CS-coefficients: one row per measure (+ an "Edge weights" row from the
#' case-drop reliability), with a group column. Handles single + group objects.
cna_stability_rows <- function(cs, cd, is_multi) {
  one_stab <- function(st, g) data.frame(
    group = g, measure = names(st$cs), cs = as.numeric(st$cs),
    row.names = NULL, stringsAsFactors = FALSE)
  one_edge <- function(c0, g) data.frame(
    group = g, measure = "Edge weights", cs = as.numeric(c0$cs),
    row.names = NULL, stringsAsFactors = FALSE)
  if (is_multi) {
    stab <- do.call(rbind, Map(one_stab, cs, names(cs)))
    edge <- if (!is.null(cd)) do.call(rbind, Map(one_edge, cd, names(cd))) else NULL
  } else {
    stab <- one_stab(cs, "")
    edge <- if (!is.null(cd)) one_edge(cd, "") else NULL
  }
  rbind(stab, edge)
}

#' Reliability & Stability tables + plot (single, group, or cluster).
cna_run_stability_nest <- function(model, options, results, is_multi) {
  if (is.null(model)) return()
  if (!isTRUE(options$stability_show_table) &&
      !isTRUE(options$stability_show_plot)) return()

  measures <- cna_nestimate_measures(options)
  if (!length(measures)) measures <- c("Betweenness", "Closeness")
  iter <- options$stability_iter
  thr  <- options$stability_threshold

  cs <- results$stability_plot$state
  if (is.null(cs)) {
    cs <- tryCatch(Nestimate::centrality_stability(
      model, measures = measures, iter = iter, threshold = thr),
      error = function(e) {
        results$stabilityTable$setNote("error", paste("Error:", e$message)); NULL
      })
    if (!is.null(cs)) results$stability_plot$setState(cs)
  }

  if (isTRUE(options$stability_show_table) && !is.null(cs)) {
    cd <- tryCatch(Nestimate::casedrop_reliability(model, iter = iter,
                   threshold = thr), error = function(e) NULL)
    rows <- cna_stability_rows(cs, cd, is_multi)
    invisible(lapply(seq_len(nrow(rows)), function(i)
      results$stabilityTable$addRow(rowKey = i, values = list(
        group = rows$group[i], measure = rows$measure[i],
        cs = as.numeric(rows$cs[i])))))

    nr <- tryCatch(Nestimate::network_reliability(model, iter = iter),
                   error = function(e) NULL)
    if (!is.null(nr) && !is.null(nr$summary) && nrow(nr$summary) > 0) {
      s <- nr$summary
      gcol <- if ("model" %in% names(s)) as.character(s$model) else rep("", nrow(s))
      invisible(lapply(seq_len(nrow(s)), function(i)
        results$reliabilityTable$addRow(rowKey = i, values = list(
          group = gcol[i], metric = as.character(s$metric[i]),
          mean = as.numeric(s$mean[i]), sd = as.numeric(s$sd[i])))))
    }
  }

  show <- isTRUE(options$stability_show_table) || isTRUE(options$stability_show_plot)
  results$stabilityTitle$setVisible(show)
  results$stabilityTable$setVisible(isTRUE(options$stability_show_table))
  results$reliabilityTable$setVisible(isTRUE(options$stability_show_table))
  results$stability_plot$setVisible(isTRUE(options$stability_show_plot))
}

#' Centrality-stability plot (splot renders net_stability / _group).
cna_splot_stability_nest <- function(cs) {
  if (is.null(cs)) return(FALSE)
  tryCatch(cograph::splot(cs),
           error = function(e) plot(1, type = "n", main = "Stability error",
                                    sub = e$message))
  TRUE
}
