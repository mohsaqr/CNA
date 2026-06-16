# CNA — Co-occurrence Network Analysis (jamovi module)

A jamovi module for **Co-occurrence Network Analysis** of sequence (long-format)
and one-hot (binary) data. Networks are estimated with the **Nestimate** package
and plotted with **cograph** (`splot`). Supports single networks, group
comparison, automatic clustering, and one-hot **mixed** (transition +
co-occurrence) networks.

> Extracted from the [JTNA](https://github.com/mohsaqr/JTNA1.2) workspace as a
> standalone module and rebuilt on the Nestimate + cograph engine (2026-06).

## Analyses

Menu group in jamovi: **Co-occurrence Network Analysis**.

| Analysis | Data type | Mode |
|----------|-----------|------|
| `CNAEventSingle`   | Sequence (long format) | Single / Group (auto when a Group variable is set) |
| `CNAEventCluster`  | Sequence (long format) | Clustering |
| `CNABinarySingle`  | One-hot (binary)       | Single / Group + **mixed** network |
| `CNABinaryCluster` | One-hot (binary)       | Clustering |

Single vs Group is **auto-detected** (presence of a Group variable) — there is
no separate "group" menu. Clustering is gated behind a *Run clustering* checkbox.

## Features (every analysis)

- **Network** — co-occurrence network, plotted with `cograph::splot` (oval layout).
- **Visualization** — network plot, edge-weight histogram, **frequencies**
  (`plot_state_frequencies`), **mosaic** (`mosaic_plot`), **distribution**
  (`distribution_plot`).
- **Centrality** — table (`cograph::centrality`) + plot with **6 styles**
  (lollipop / bar / line / dot / heatmap / distribution).
- **Community detection** — `cograph::detect_communities`
  (louvain / walktrap / fast_greedy / label_prop / infomap / leiden).
- **Reliability & Stability** — centrality-stability CS-coefficients
  (`centrality_stability`), edge-weight case-drop stability
  (`casedrop_reliability`), split-half reliability (`network_reliability`),
  plus a stability plot.
- **Bootstrap** — edge bootstrap (`bootstrap_network`).
- **Permutation** — between-group/cluster edge permutation test (`permutation`).
- **Compare** — descriptive comparison of two groups/clusters (`compare_model`).
- **One-hot mixed** (binary single) — `wtna(method = "both")`: transition +
  co-occurrence networks side by side.

Group networks and clusters are built **without manual loops** — via Nestimate's
native `group=` argument and `cluster_network()`, both returning a
`netobject_group` that flows through every downstream analysis and `splot` grid.

## Engine

```
data ─▶ Nestimate::build_network(method="co_occurrence")   # single / group
      └ Nestimate::cluster_network(method="co_occurrence")  # clustering
      └ Nestimate::wtna(method="both")                      # one-hot mixed
            │  -> netobject / netobject_group / wtna_mixed
            ▼
analysis:  cograph::centrality · cograph::detect_communities ·
           Nestimate::bootstrap_network · permutation · compare_model ·
           centrality_stability / casedrop_reliability / network_reliability
            ▼
plots:     cograph::splot (+ plot_centrality / mosaic_plot /
           distribution_plot / plot_state_frequencies)
```

See [`docs/NESTIMATE_COOCCURRENCE_STUDY.md`](docs/NESTIMATE_COOCCURRENCE_STUDY.md)
for a deep dive on how Nestimate computes co-occurrence.

## Structure

```
DESCRIPTION          # package metadata (Imports: Nestimate, cograph, + closure)
NAMESPACE
R/                   # jamovi *.h.R (generated) + *.b.R (backends) + cna_utils.R
  cna_utils.R        # the Nestimate engine: builders + analyses + splot helpers
jamovi/              # *.a.yaml (options) *.u.yaml (UI) *.r.yaml (results) 0000.yaml
data/                # bundled example dataset (Regulation_long.csv)
docs/                # NESTIMATE_COOCCURRENCE_STUDY.md
```

The UI (`*.u.yaml`) follows the JTNA arrangement: inline build options at the
top, `━━ section ━━` dividers, and two-column `Label`-grouped CollapseBoxes.

## Building / installing

**Build with jamovi's own R** (it already carries the dependency tree, so the
bundled packages match jamovi's ABI):

```bash
"/Applications/jamovi.app/Contents/Frameworks/R.framework/Versions/<ver>-<arch>/Resources/bin/Rscript" \
  -e 'jmvtools::install(".")'
```

This compiles the module, installs it into
`~/Library/Application Support/jamovi/modules/CNA`, and writes `CNA_1.0.0.jmo`.
Then **restart jamovi**.

> Do **not** build with the system R if its package library differs from
> jamovi's — jmvtools would skip bundling deps that jamovi lacks, producing an
> unusable `.jmo`. See `LEARNINGS.md`.

## Dependencies

`Nestimate` (GitHub: `mohsaqr/Nestimate`), `cograph` (GitHub: `mohsaqr/cograph`),
`jmvcore`, `R6`, `ggplot2`, `patchwork`, plus Nestimate's tree
(`glasso`, `data.table`, `cluster`, `brglm2`, …) and `igraph` / `stringdist`
(used by cograph communities / Nestimate clustering dissimilarities).

## License

MIT — see [LICENSE](LICENSE).
