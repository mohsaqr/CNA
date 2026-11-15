# CNA — Co-occurrence Network Analysis (jamovi module)

A jamovi module for **Co-occurrence Network Analysis** of event (long-format),
one-hot (binary) and delimited (keyword) data. Networks are estimated with
**Nestimate** and drawn with **cograph**. Single networks, group comparison and
clustering, with centrality and its stability, communities, cliques, backbones,
link prediction, association rules, bootstrap and permutation inference.

> Extracted from the [JTNA](https://github.com/mohsaqr/JTNA1.2) workspace and
> rebuilt on the Nestimate + cograph engine. Version 1.1.0 (2026-11) restores
> every function the earlier tna-based CNA had and adds what Nestimate 1.0 and
> cograph 2.7 offer for co-occurrence networks.

## Analyses

Menu group in jamovi: **Co-occurrence Network Analysis**.

| Analysis | Data | Mode |
|----------|------|------|
| `CNAEventSingle`   | Events, long format (Action, Actor, Time, Order) | Single, or one network per group |
| `CNAEventCluster`  | Events, long format | Clusters of sessions (actor × time-gap session), one network per cluster |
| `CNABinarySingle`  | One-hot (binary) columns | Single, or one network per group |
| `CNABinaryCluster` | One-hot (binary) columns | Clusters, one network per cluster |
| `CNADelimited`     | Item lists (keywords "a; b; c"), or one item per row with a transaction ID | Single, with similarity normalisation |

Group mode switches on when a Group variable is assigned. Clustering waits for
*Run clustering*; *Choosing k* works before that.

The module has **no sequence analyses**: co-occurrence ignores order, so there
are no sequence index/heatmap plots, no state-distribution-over-time plot and
no subsequence comparison. The one order-based view is the optional mixed
network, which shows the transition network beside the co-occurrence network.

### Options panel

Every analysis uses the same layout:

- **Variables**, then the open top area: **Network** (what defines the network:
  time threshold, window, scaling or the delimited-data settings; *Clustering*
  first in the cluster analyses), **Output** (network plot, co-occurrence
  matrix, network summary) and **Reproducibility** (seed).
- Collapsed sections under three headings: *Network Analysis* (Plots,
  Centrality, Edge Betweenness, Communities, Cliques, Backbone, Link
  Prediction, Association Rules), *Confirmatory Testing* (Bootstrap,
  Reliability & Stability, Permutation Test) and *Group / Cluster Comparison*.
- Each section has an **Output** column (what to show) and a **Settings** or
  **Method** column (what to compute). Options that belong to one choice sit
  indented under it (a table's row limit under *Table*, a backbone method's
  parameter under that method) and are greyed until it is chosen.

`Rscript tools/preview_panels.R [out.html] [before-dir]` draws the panels as
HTML from the definition files, for review without jamovi.

## What each analysis offers

| Section | What it does | Engine |
|---|---|---|
| **Network** | weight matrix table, network summary table, network plot | `build_network`, `cooccurrence`, `summary`, `splot` |
| **Plots** | edge-weight histogram, code frequencies, mosaic plot; network plot appearance | `plot_state_frequencies`, `mosaic_plot`, `splot` |
| **Mixed network** (event and one-hot analyses) | the transition network next to the co-occurrence network, side by side (Transition \| Co-occurrence per network) or overlaid. Transitions are counted within the same units and settings as the co-occurrences (events: the same sessions; one-hot: the same actor × session and windows) | `build_network(method = "frequency")`, `splot`, `plot_mixed_network` |
| **Plot appearance** | layout, edge labels, hide edges below, cut value, edge-label / node / node-label size — shared by every network plot | `cograph::splot` |
| **Centrality** | 14 measures by family — degree: strength, degree, coreness; paths: betweenness, RSP betweenness, closeness, harmonic, eccentricity; spectral: eigenvector, PageRank, diffusion (tna's, range-normalised to 0–1); local structure: clustering (weighted, Zhang & Horvath, 0–1), constraint, effective size; loops, normalisation, inverted weights; 6 plot styles | `cograph::centrality`; Diffusion `net_centrality` |
| **Reliability & stability** | CS-coefficients for every selected measure, edge-weight case-drop, split-half reliability; correlation method and certainty | `centrality_stability`, `casedrop_reliability`, `network_reliability` |
| **Communities** | Louvain, Leiden, walktrap, fast greedy, Infomap, label propagation, edge betweenness, leading eigenvector, spinglass, optimal; resolution γ | `cograph::communities` |
| **Edge betweenness** | table and network plot | `net_edge_betweenness` |
| **Cliques** | cliques of a chosen size above a weight threshold (table + grid plot) and the clique-complex summary (clique counts, Euler characteristic, Betti numbers, Q-analysis) | `build_simplicial`, `betti_numbers`, `q_analysis` |
| **Backbone** | disparity filter, weight threshold or lowest proportion; kept/removed table and backbone plot | `net_prune` |
| **Link prediction** | resource allocation, Adamic-Adar, common neighbours, Jaccard, preferential attachment, Katz | `predict_links` |
| **Association rules** | Apriori rules with support, confidence, lift and conviction; table and plot | `association_rules` |
| **Bootstrap** | edge stability or threshold inference; percentile or basic intervals | `bootstrap_network` |
| **Permutation** (groups/clusters) | edge differences between every pair of networks; p-value correction over the unique edges | `permutation` |
| **Compare** (groups/clusters) | summary metrics, network properties, edge differences, comparison plot (5 types), network difference plot | `compare_model`, `plot_difference` |
| **Clustering** (cluster analyses) | 9 distances, 9 linkages, cluster diagnostics, choosing k | `build_clusters`, `cluster_diagnostics`, `cluster_choice` |
| **Reproducibility** | *Set seed* fixes every resampling and randomised step | — |

### Notes on the statistics

- **Resampling needs units.** Bootstrap, permutation and stability resample
  sessions for event data (each actor's events split by the time threshold;
  sessions of one actor are drawn independently) and actor × session units
  for one-hot data. Without an Actor (or a Session) the data are one unit:
  intervals would have zero width, CS coefficients would be 0 and every
  p-value 1. The module therefore does not run those analyses and says why.
- **Permutation p-value correction** runs over the unique edges of each
  comparison. Nestimate 1.0.7 corrects over all n × n cells, which counts every
  undirected edge twice; the module requests unadjusted p-values and applies
  the correction itself, and the plot's stars use the same values.
- **Path-based centralities invert the weights** by default (a frequent
  co-occurrence is a short distance), as Nestimate does — eccentricity
  included (cograph ignores the switch for it, so the module gives it the
  lengths 1/w). Effective size never counts self-loops (Burt, 1992).
- **Diffusion** is tna's diffusion centrality (row sums of W + W² + … + Wⁿ),
  range-normalised to 0–1 as Nestimate does by default: on raw co-occurrence
  counts the power series runs into the billions.
- **Clustering** is weighted: tna's coefficient (Zhang & Horvath, 2005) on the
  weights divided by their maximum, so it lies in 0–1 (tna's value on raw
  counts divided by 2 × the largest weight). The unweighted coefficient is 1
  for every node of a complete network, which most co-occurrence networks
  are.
- **Stability assesses the table's own values**: each resample is scored with
  the same definition as the centrality table. A centrality that does
  not vary across nodes (for example degree in a complete network) has no
  stability coefficient; the table says so.
- **Link prediction** ignores the diagonal (a code co-occurring with itself is
  not a link). Jaccard compares neighbour sets; the other scores use weights.
- **Groups** are ordered by the Group's levels (jamovi's level order),
  otherwise numerically or alphabetically, never by row order; comparisons
  read "first vs second" in that order.
- **Transaction ID mode** (Delimited): every row holds one item per Item
  column; items are pooled over the columns, trimmed, and blank items and rows
  without an ID are left out (counted in a note). The "relative"
  (row-normalised) similarity is not offered: it is directed, and Nestimate
  1.0.7 returns it as an undirected network with half its edges.
- **One-hot data**: code columns must hold 0/1 (numbers, logicals or 0/1
  labels); rows must be in time order within each actor/session. Sliding
  windows count as in tna's `import_onehot(window_type = "sliding")`. Binary
  Clustering clusters one windowed sequence per actor × session and builds
  each cluster's network from its own rows (so scaling applies). Stability,
  bootstrap and permutation resample whole actors/sessions and keep their
  windows. See [`docs/ONEHOT_AUDIT.md`](docs/ONEHOT_AUDIT.md) for the audit and
  what Nestimate `a1754ab` fixed.
- **Scaled models are never resampled**: bootstrap, permutation and stability
  are defined on raw counts; set Scaling to "No scaling" for them.
- **Equivalence.** Every result is checked against independent references
  (base-R formulas, igraph, arules, stringdist, cluster, tna) in
  `tests/testthat/test-equiv-*.R`; the findings, the defects fixed and the
  upstream defects still open are in
  [`docs/EQUIVALENCE_AUDIT.md`](docs/EQUIVALENCE_AUDIT.md).

## Example datasets

Bundled in `data/`, shown in jamovi's Data Library, and rebuilt by
`Rscript tools/make_datasets.R` (seeded; the simulated data come from simulab and
Nestimate):

| Dataset | Source | Format | Use it with |
|---|---|---|---|
| `Regulation_long` | real | events (Actor, Group, Time, Action) | Event Data |
| `Regulation_sessions` | real | `;`-separated action lists per session | Delimited (Item Columns = Actions) |
| `SRL_events` | simulated (simulab) | self-regulated learning event log: 8 actions, 200 students in two groups with different learning cycles (Strategic: plan → read → summarise → monitor → evaluate → reflect; Surface: read → review → practice), Course, Achievement, Time, Order | Event Data, group comparison, mixed network |
| `SRL_onehot` | simulated (simulab) | the same log coded **0/1** per episode of 4 events (8 action columns, Actor, Group, Episode) | Binary Data (Actor = Actor, window 1), mixed network |
| `LearningActivities` | simulated (Nestimate) | **0/1** indicators of 6 online activities, 200 students × 30 time points | Binary Data (Actor = Student), clustering |
| `HumanAI_turns` | real (Saqr, 2026) | human–AI pair-programming turns: codes as a `;` list (`Codes`) **and** as 17 **0/1** columns, with Session, Project, Speaker, Turn | Delimited (Codes) and Binary Data (Actor = Session, Group = Speaker) |

The test suite checks the planted SRL structure is recovered (Strategic: Reflect →
Plan, Summarize → Monitor; Surface: Read → Review) and that every dataset runs in
its intended analyses.

## Re-run and caching

Each result declares the options that change it; its `clearWith` is derived from
that by the generator (below). Expensive computations are cached in hidden state
results whose `clearWith` excludes display options:

| You change… | What recomputes |
|---|---|
| A model option (variables, scaling, threshold, window, clustering, seed for clusters) | the network(s) and every analysis |
| An analysis option (measures, method, iterations, …) | that analysis only — the network is reused |
| A display option (plot appearance, plot type, row limits) | only the affected plot or table — nothing is recomputed |

## Structure

```
DESCRIPTION, NAMESPACE
R/
  cna_build.R      data + options -> Nestimate network(s)
  cna_compute.R    analyses as pure functions returning tidy data frames
  cna_render.R     plot renderers (one per result image)
  cna_run.R        orchestration: input checks, caching, filling the results
  <Analysis>.b.R   generated thin backends (route .run and renderers to the engine)
  <Analysis>.h.R   generated by jmvtools from the YAML
jamovi/            generated *.a.yaml, *.r.yaml, *.u.yaml and 0000.yaml
tools/
  generate_jamovi.R   single source of the five analyses' definitions
  build_module.R      build the .jmo (optionally bundling a local engine)
  make_datasets.R     derives data/Regulation_sessions.csv
data/              Regulation_long.csv, Regulation_sessions.csv
tests/             testthat suite (run_tests.R + testthat/)
docs/              NESTIMATE_COOCCURRENCE_STUDY.md
```

**Do not edit `jamovi/*.yaml` or `R/*.b.R` by hand.** Change
`tools/generate_jamovi.R`, then regenerate:

```bash
Rscript tools/generate_jamovi.R
"<jamovi Rscript>" -e 'jmvtools::prepare(".")'   # regenerates R/*.h.R
```

The generator refuses to write definitions in which a `clearWith` names an
undeclared option, a control binds an unknown option, or an option has no
control. The test suite checks the committed files are exactly its output.

## Building and installing

Use **jamovi's own R** (`/Applications/jamovi.app/Contents/Frameworks/R.framework/Versions/<ver>-<arch>/Resources/bin/Rscript`):

```bash
"<jamovi Rscript>" tools/build_module.R
# or bundle an engine that is newer locally than on GitHub:
"<jamovi Rscript>" tools/build_module.R --engine ~/Documents/Github/Nestimate
```

This writes `CNA_<version>.jmo` without touching jamovi. `jmc` installs the
DESCRIPTION Remotes (Nestimate, cograph) from GitHub HEAD; `--engine <dir>`
installs that source into the build library and replaces the GitHub copy inside
the `.jmo`. The script lists the bundled package versions at the end.

Install by sideloading the `.jmo` in jamovi (*Modules ▸ Sideload*), or unzip it
into the modules directory (the archive layout matches it) and restart jamovi:

```bash
unzip -o CNA_1.1.0.jmo -d "$HOME/Library/Application Support/jamovi/modules/"
```

## Testing

```bash
"<jamovi Rscript>" tests/run_tests.R
```

`run_tests.R` puts `build/R*/` (CNA, Nestimate, cograph) on the library path and
reports every failure. The suite (748 expectations) covers:

- **Equivalence** with the wrapped engine call for every builder and analysis,
  and with an independent computation where one exists (igraph cliques and edge
  betweenness, `p.adjust` over the unique edges, raw weight differences,
  strength as row sums, Euler characteristic = alternating f-vector sum =
  alternating Betti sum).
- **Integration**: all five analyses through jmvcore with every section on; no
  section may fail, every visible table has rows or a note saying why not, every
  visible image renders.
- **Regressions**: ISO time stamps keep their time of day, cluster scaling,
  transaction IDs with a separator, spinglass γ, stale section messages, plots
  that would be blank, resampling without units.
- **Manifest invariants**: generator idempotence, no ghost options, no undeclared
  option reads, display options never clear a cached computation, each output
  clears with the state it draws, every renderer exists.
- **Plot output** read back from SVG: integer labels for counts, label on/off,
  hidden weak edges.

## Dependencies

`Nestimate` (>= 1.0.5, GitHub `mohsaqr/Nestimate`) and `cograph` (>= 2.7.11,
GitHub `mohsaqr/cograph`), `jmvcore`, `R6`, `ggplot2`, `patchwork`, and their
closure (`idiographic`, `psychnets`, `igraph`, `stringdist`, `brglm2`, …),
declared in full because jamovi does not resolve Remotes transitively.

## License

MIT — see [LICENSE](LICENSE).
