#!/usr/bin/env Rscript
# tools/preview_panels.R
#
# Draws the options panel of every analysis as HTML, from the jamovi
# definitions (<A>.u.yaml for the layout, <A>.a.yaml for titles and defaults),
# so a layout change can be reviewed without opening jamovi. Controls show
# their default state; a control whose `enable` is false at the defaults is
# greyed, as jamovi shows it. Sections start collapsed, as in jamovi.
#
# Usage (default output: tmp/cna_panels.html):
#   Rscript tools/preview_panels.R [out.html] [before-dir]
# With a before-dir (a folder holding older <A>.u.yaml / <A>.a.yaml files) the
# old and new panels are drawn side by side.

# the repository root is the parent of this script's folder
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg), "..")) else "."
stopifnot("the repository root must hold the DESCRIPTION" =
            file.exists(file.path(ROOT, "DESCRIPTION")))

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1L) args[[1]] else file.path(ROOT, "tmp", "cna_panels.html")
before_dir <- if (length(args) >= 2L) args[[2]] else NULL

ANALYSES <- c("CNAEventSingle", "CNAEventCluster", "CNABinarySingle",
              "CNABinaryCluster", "CNADelimited")

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Option lookup and the value of every name an `enable` expression can use:
# Bool options and RadioButton names (<option>_<part>) at their defaults.
read_options <- function(a_file) {
  opts <- yaml::read_yaml(a_file)$options
  stats::setNames(opts, vapply(opts, `[[`, "", "name"))
}

enable_values <- function(opts) {
  bools <- Filter(\(o) identical(o$type, "Bool"), opts)
  lists <- Filter(\(o) identical(o$type, "List"), opts)
  radio <- unlist(unname(lapply(lists, \(o) {
    parts <- vapply(o$options, `[[`, "", "name")
    stats::setNames(parts == o$default, paste0(o$name, "_", parts))
  })))
  c(vapply(bools, \(o) isTRUE(o$default), logical(1)), radio)
}

is_enabled <- function(expr, values) {
  if (is.null(expr)) return(TRUE)
  tokens <- unique(regmatches(expr, gregexpr("[A-Za-z_][A-Za-z0-9_]*", expr))[[1]])
  env <- list2env(as.list(stats::setNames(
    lapply(tokens, \(t) if (t %in% names(values)) values[[t]] else TRUE), tokens)))
  isTRUE(eval(parse(text = expr), env))
}

list_title <- function(o, part) {
  hit <- Filter(\(x) identical(x$name, part), o$options)
  if (length(hit)) hit[[1]]$title else part
}

# Divider labels are drawn by the generator as "━━ <bold text> ━━".
plain_divider <- function(label) {
  cp <- utf8ToInt(gsub("━", "", label))
  upper <- cp >= 0x1D5D4L & cp <= 0x1D5EDL
  lower <- cp >= 0x1D5EEL & cp <= 0x1D607L
  cp[upper] <- cp[upper] - 0x1D5D4L + 65L
  cp[lower] <- cp[lower] - 0x1D5EEL + 97L
  trimws(intToUtf8(cp))
}

render_node <- function(n, opts, values, parent_on = TRUE) {
  on <- parent_on && is_enabled(n$enable, values)
  cls <- if (on) "" else " off"
  kids <- function(parent_on = on)
    paste(vapply(n$children %||% list(), render_node, "", opts = opts, values = values,
                 parent_on = parent_on), collapse = "")
  o <- opts[[n$optionName %||% n$name %||% ""]]
  switch(n$type,
    VariableSupplier = sprintf(
      '<div class="supplier"><div class="varlist">Variables</div><div class="targets">%s</div></div>',
      kids()),
    TargetLayoutBox = sprintf('<div class="target"><div class="tlabel">%s</div><div class="tbox"></div></div>',
                              esc(n$label)),
    LayoutBox = if (length(n$children) > 1L && all(vapply(n$children, \(c) identical(c$type, "LayoutBox"), logical(1))))
      sprintf('<div class="columns">%s</div>', kids())
      else sprintf('<div class="col">%s</div>', kids()),
    Label = if (grepl("━", n$label))
      sprintf('<div class="divider">%s</div>', esc(plain_divider(n$label)))
      else sprintf('<div class="group"><div class="glabel">%s</div>%s</div>', esc(n$label), kids()),
    CollapseBox = sprintf('<details class="collapse"><summary>%s</summary><div class="cbody">%s</div></details>',
                          esc(n$label), kids()),
    CheckBox = {
      checked <- isTRUE(o$default)
      sprintf('<div class="ctl%s"><span class="box">%s</span> %s</div>%s', cls,
              if (checked) "&#9745;" else "&#9744;", esc(n$label %||% o$title),
              if (length(n$children)) sprintf('<div class="sub">%s</div>', kids()) else "")
    },
    RadioButton = {
      chosen <- identical(o$default, n$optionPart)
      sprintf('<div class="ctl%s"><span class="box">%s</span> %s</div>%s', cls,
              if (chosen) "&#9673;" else "&#9711;", esc(list_title(o, n$optionPart)),
              if (length(n$children)) sprintf('<div class="sub">%s</div>', kids()) else "")
    },
    ComboBox = sprintf('<div class="ctl%s">%s <span class="combo">%s &#9662;</span></div>', cls,
                       esc(n$label %||% o$title), esc(list_title(o, o$default))),
    TextBox = sprintf('<div class="ctl%s">%s <span class="text">%s</span>%s</div>', cls,
                      esc(n$label %||% o$title), esc(format(o$default %||% "")),
                      if (!is.null(n$suffix)) paste0(" ", esc(n$suffix)) else ""),
    LevelSelector = sprintf('<div class="ctl%s">%s <span class="combo">&mdash; &#9662;</span></div>',
                            cls, esc(o$title)),
    VariablesListBox = "",
    sprintf('<div class="ctl">[%s %s]</div>', esc(n$type), esc(n$name %||% "")))
}

render_panel <- function(dir, analysis) {
  u_file <- file.path(dir, paste0(analysis, ".u.yaml"))
  a_file <- file.path(dir, paste0(analysis, ".a.yaml"))
  if (!file.exists(u_file) || !file.exists(a_file)) return('<p class="none">(not present)</p>')
  opts <- read_options(a_file)
  values <- enable_values(opts)
  u <- yaml::read_yaml(u_file)
  sprintf('<div class="panel"><div class="ptitle">%s</div>%s</div>', esc(u$title),
          paste(vapply(u$children, render_node, "", opts = opts, values = values),
                collapse = ""))
}

css <- '
:root { --bg:#f3f3f3; --panel:#ffffff; --ink:#333; --muted:#888; --line:#d9d9d9;
        --accent:#3e6da9; --head:#e9eef5; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#1e1e1e; --panel:#2a2a2a; --ink:#ddd; --muted:#888; --line:#444; --accent:#8fb3e6; --head:#33393f; } }
:root[data-theme="dark"] { --bg:#1e1e1e; --panel:#2a2a2a; --ink:#ddd; --muted:#888; --line:#444;
  --accent:#8fb3e6; --head:#33393f; }
body { background:var(--bg); color:var(--ink); font:13px/1.45 -apple-system, "Segoe UI", Roboto, sans-serif;
       margin:0; padding:16px; }
h1 { font-size:18px; margin:0 0 4px; } p.lead { margin:0 0 12px; color:var(--muted); max-width:70ch; }
.tabs { display:flex; flex-wrap:wrap; gap:6px; margin-bottom:12px; }
.tabs button { border:1px solid var(--line); background:var(--panel); color:var(--ink); padding:5px 10px;
               border-radius:4px; cursor:pointer; font:inherit; }
.tabs button[aria-selected="true"] { border-color:var(--accent); color:var(--accent); font-weight:600; }
.tools { margin-bottom:12px; } .tools button { font:inherit; padding:4px 10px; cursor:pointer; }
.compare { display:flex; gap:16px; align-items:flex-start; flex-wrap:wrap; }
.side { flex:1 1 440px; min-width:0; max-width:560px; } .side h2 { font-size:14px; margin:0 0 6px; }
.panel { background:var(--panel); border:1px solid var(--line); border-radius:6px; padding:12px 14px; }
.ptitle { font-size:16px; font-weight:600; margin-bottom:10px; }
.supplier { display:flex; gap:10px; margin-bottom:12px; }
.varlist { flex:0 0 32%; border:1px solid var(--line); min-height:120px; padding:4px; color:var(--muted); }
.targets { flex:1; } .target { margin-bottom:6px; } .tlabel { font-weight:600; font-size:12px; }
.tbox { border:1px solid var(--line); height:20px; }
.columns { display:flex; gap:14px; margin:6px 0 10px; } .col { flex:1 1 0; min-width:0; }
.group { margin-bottom:8px; } .glabel { font-weight:600; margin-bottom:2px; }
.ctl { padding:1px 0 1px 2px; } .ctl.off { color:var(--muted); opacity:.6; }
.sub { margin-left:20px; } .box { display:inline-block; width:14px; }
.combo, .text { border:1px solid var(--line); padding:0 6px; border-radius:3px; white-space:nowrap; }
.divider { margin:14px 0 6px; font-weight:700; letter-spacing:.02em; color:var(--accent);
           border-bottom:1px solid var(--line); padding-bottom:2px; }
details.collapse { border-top:1px solid var(--line); }
details.collapse > summary { cursor:pointer; padding:7px 4px; font-weight:600; background:var(--head); margin:0 -4px; }
.cbody { padding:6px 0 4px; }
.none { color:var(--muted); }
@media (max-width:640px) { .columns { flex-direction:column; gap:0; } .supplier { flex-direction:column; } }
'

js <- '
function show(id) {
  document.querySelectorAll(".analysis").forEach(el => el.hidden = el.id !== id);
  document.querySelectorAll(".tabs button").forEach(b => b.setAttribute("aria-selected", b.dataset.id === id));
}
function expandAll(open) { document.querySelectorAll("details.collapse").forEach(d => d.open = open); }
'

sections <- vapply(ANALYSES, function(a) {
  after <- render_panel(file.path(ROOT, "jamovi"), a)
  body <- if (is.null(before_dir)) after else
    sprintf('<div class="compare"><div class="side"><h2>Before</h2>%s</div><div class="side"><h2>After</h2>%s</div></div>',
            render_panel(before_dir, a), after)
  sprintf('<section class="analysis" id="%s"%s>%s</section>', a,
          if (a == ANALYSES[[1]]) "" else " hidden", body)
}, "")
tabs <- paste(sprintf('<button data-id="%s" aria-selected="%s" onclick="show(\'%s\')">%s</button>',
                      ANALYSES, ANALYSES == ANALYSES[[1]], ANALYSES, ANALYSES), collapse = "")

html <- paste0(
  '<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>CNA Options Panels</title><style>', css, '</style></head><body>',
  '<h1>CNA options panels</h1>',
  '<p class="lead">Drawn from the jamovi definition files. Controls show their defaults; ',
  'greyed controls are disabled until their parent option is ticked. Sections start ',
  'collapsed, as in jamovi.</p>',
  '<div class="tabs">', tabs, '</div>',
  '<div class="tools"><button onclick="expandAll(true)">Expand all</button> ',
  '<button onclick="expandAll(false)">Collapse all</button></div>',
  paste(sections, collapse = ""), '<script>', js, '</script></body></html>')

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
con <- file(out_file, open = "w", encoding = "UTF-8")
writeLines(html, con, useBytes = TRUE)
close(con)
cat("Wrote", out_file, "\n")
