# Plot output, read back from SVG text (cograph labels are <text> nodes).
# Pins the edge-label behaviour verified for cograph 2.7.12: counts print as
# integers, the "Edge labels" option shows and hides labels on both network
# objects and bare matrices, and appearance options reach splot().

numeric_labels <- function(texts) grep("^-?[0-9]*\\.?[0-9]+\\**$", texts, value = TRUE)
look <- function(...) modifyList(list(plot_layout = "oval", plot_edge_labels = TRUE,
                                      plot_minimum = 0, plot_cut = 0,
                                      plot_edge_label_size = 0.8, plot_node_size = 0,
                                      plot_label_size = 0), list(...))

test_that("count weights are labelled as integers, without '.00'", {
  net <- fixture_event_net()
  labels <- numeric_labels(svg_texts(cna(".cna_splot")(net, look())))
  expect_gt(length(labels), 0)
  expect_false(any(grepl("\\.00$", labels)))
  w <- net$weights
  expect_true(all(unique(as.numeric(labels)) %in% unique(as.vector(w))))
})

test_that("'Edge labels' off hides the labels of a network object and a matrix", {
  net <- fixture_event_net()
  off <- look(plot_edge_labels = FALSE)
  expect_length(numeric_labels(svg_texts(cna(".cna_splot")(net, off))), 0)
  m <- unclass(net$weights)
  expect_length(numeric_labels(svg_texts(cna(".cna_splot")(m, off))), 0)
  # a bare matrix is unlabelled by default in cograph; "on" must label it
  expect_gt(length(numeric_labels(svg_texts(cna(".cna_splot")(m, look())))), 0)
})

test_that("'Hide edges below' removes the weaker edges' labels", {
  net <- fixture_event_net()
  w <- net$weights
  cut <- stats::quantile(w[upper.tri(w)], 0.75)
  labels <- as.numeric(numeric_labels(svg_texts(cna(".cna_splot")(net, look(plot_minimum = cut)))))
  expect_gt(length(labels), 0)
  expect_true(all(labels >= cut))
})

test_that("scaled weights keep their decimals", {
  dat <- make_event_data()
  net <- cna("cna_build_event_nest")(dat, opts_event(buildModel_scaling = "minmax"))
  labels <- numeric_labels(svg_texts(cna(".cna_splot")(net, look())))
  expect_true(any(grepl("\\.", labels)))
  expect_equal(cna(".cna_digits")(net$weights), 2L)
  expect_equal(cna(".cna_digits")(fixture_event_net()$weights), 0L)
})

test_that("a renderer failure is drawn into the image, not raised", {
  self <- list(options = list(), results = list())
  expect_error(cna("cna_render")(self, "no_such_plot"), class = "cna_engine_contract")
  broken <- list(options = list(plot_layout = "oval"),
                 results = list(buildModelContent = list(state = structure(list(), class = "netobject"))))
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_true(cna("cna_render")(broken, "buildModel_plot"))
})

test_that("panel grids cover n plots", {
  g <- cna(".cna_grid")
  invisible(lapply(1:12, function(n) expect_gte(prod(g(n)), n)))
  expect_equal(g(1), c(1, 1))
  expect_equal(g(4), c(2, 2))
})
