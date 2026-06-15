# How Nestimate Computes Co-occurrence — A Detailed Study

> Reference notes for the standalone CNA module. Studies the co-occurrence
> algorithm(s) in the sibling `Nestimate` package
> (`/Users/mohammedsaqr/Documents/Github/Nestimate`) and contrasts them with
> how this CNA jamovi module produces the same kind of network.

## TL;DR

Nestimate has **two independent co-occurrence systems**:

1. **`cooccurrence()`** (`R/co_occurrence.R`) — a general-purpose,
   transaction/bibliometric engine. Items that appear in the same *transaction*
   (document, basket, session) are linked. Base-R only. 6 input formats, 8
   similarity normalizations.
2. **`build_network(method = "co_occurrence")`** (`R/estimators.R`,
   `R/build_network.R`) — a *sequence-positional* engine. States that appear
   within the same sliding *window* of a sequence are linked.

Both return the **same `netobject`** S3 type, so all downstream analysis
(centrality, bootstrap, plotting, comparison) is method-agnostic.

By contrast, **this CNA jamovi module computes nothing itself** — it delegates
to the `tna` package (`tna::build_model`, `tna::group_model`). See the last
section.

---

## System A — `cooccurrence()` (transaction-based)

File: `Nestimate/R/co_occurrence.R` (~395 lines). Argument names follow the
"citenets" convention. The whole thing is a clean **5-stage pipeline**:

```
detect format → parse to transactions → filter → count → normalize → threshold/top_n → wrap
```

### 1. Format auto-detection — `.co_detect_format(data, field, by, sep)`

Dispatches on the *combination* of arguments and data class:

| Format | Trigger |
|--------|---------|
| `list` | `data` is a bare list |
| `multi_delimited` | `sep` + `field` (length > 1) |
| `delimited` | `sep` + `field` (single col) |
| `long` | `field` + `by` |
| `binary` | matrix/df, all values ∈ {0,1} |
| `wide` | matrix/df, non-binary |

Each format has **one dedicated parser** (`.co_parse_*`) and they all return the
**same intermediate representation**: a *list of character vectors*, one per
transaction, de-duplicated and trimmed. This is the key design move — every
input shape collapses to one canonical structure before any math happens.

### 2. Filtering

- Drop empty transactions.
- `min_occur`: drop entities appearing in fewer than *k* transactions (computed
  via `table(unlist(transactions))`), then drop transactions emptied by that.

### 3. Counting — the algorithmic heart

```r
.co_transactions_to_matrix(transactions)  # B: rows = transactions, cols = items (logical)
.co_compute_matrix(B)                      # C = crossprod(B * 1L)   -> items x items
```

The entire co-occurrence count is **one `crossprod`**:
`C = Bᵀ B`. `C[i, j]` = number of transactions in which items *i* and *j* both
appear. The **diagonal `C[i, i]` = item frequency `fᵢ`** (number of transactions
containing *i*). Fully vectorized, no loops.

### 4. Normalization — `.co_normalize(C, freq, N, method)`

Eight measures, all built from `C` and the frequency vector `freq`, using
`outer()` for the pairwise denominators (again, no loops). With `f` = freq,
`C` = co-occurrence count:

| `similarity` | Formula | Name |
|--------------|---------|------|
| `none` | `C` | raw counts |
| `jaccard` | `C / (fᵢ + fⱼ − C)` | Jaccard |
| `cosine` | `C / √(fᵢ·fⱼ)` | Salton cosine |
| `inclusion` | `C / min(fᵢ, fⱼ)` | Simpson |
| `association` | `C / (fᵢ·fⱼ)` | assoc. strength (van Eck & Waltman 2009) |
| `dice` | `2C / (fᵢ + fⱼ)` | Dice |
| `equivalence` | `C² / (fᵢ·fⱼ)` | cosine² |
| `relative` | row-normalized (`C / rowSums`) | conditional prob. |

Every denominator guards `denom[denom == 0] <- 1` to avoid divide-by-zero.

### 5. Post-filters + wrap

- `threshold`: zero out weights below a cutoff (applied *after* normalization).
- `top_n`: keep only the top-N edges by weight (cutoff from the sorted upper
  triangle, applied symmetrically to both triangles).
- `.wrap_netobject(W, method = "co_occurrence_fn", directed = FALSE)` →
  `netobject`. `$params` records similarity, threshold, min_occur, diagonal,
  top_n, n_transactions, n_items.

---

## System B — `build_network(method = "co_occurrence")` (sequence-positional)

File: `Nestimate/R/estimators.R`, fn `.estimator_co_occurrence`. Registered in
`R/estimator_registry.R` under names `co_occurrence` (aliases `cna`/`wcna`) and
`wtna_cooccurrence`.

This path answers a *different* question: in **sequence** data (a trajectory of
states over time), which states occur near each other? Routing:

- **matrix input** → `.transition_matrix_input(data, type = "co_occurrence")`.
- **wide + `codes`** or **auto-detected binary one-hot** →
  `.estimator_wtna_core(..., wtna_method = "cooccurrence")`.
- **wide sequence** → `.count_cooccurrence_wide()`.
- **long sequence** → `.count_cooccurrence_long()`.

### Windowed counting — `.count_cooccurrence_wide`

The sequence is sliced into **non-overlapping windows** of
`effective_window = window_size × window_span` columns. Within each window, it
counts *all* ordered position pairs `(j, k)` (including self-pairs), incrementing
`C[from_state, to_state]`. Because both `(j,k)` and `(k,j)` are counted, the
result is symmetric → undirected. A `weighted = TRUE` option down-weights each
contribution by `1 / sequence_length`.

Note: this counter uses nested `for` loops over windows × positions × positions
(inside a private helper). It is the positional analogue of System A's single
`crossprod`, trading vectorization for explicit window/position control.

---

## Shared infrastructure (the part worth copying)

- **`.wrap_netobject(mat, data, method, directed, inits)`** (`R/mcml.R`) — the
  single return type. Validates squareness/numeric, builds `edges` (via
  `.extract_edges_from_matrix`) and a tidy `nodes_df`
  (`id, label, name, x, y`), stores `weights`, `method`, `directed`. **Every
  estimator funnels through this**, which is exactly why centrality / bootstrap /
  plot / compare don't care how the matrix was produced.
- **`.void_markers <- c("%", "*", "", "NA", "NaN")`** (`R/estimators.R`) —
  values treated as missing in sequences (used by the wide parser).
- **`estimator_registry`** — `register_estimator(name, fn, ...)`;
  `build_network()` is a thin dispatcher keyed by `method`. New methods plug in
  without editing the dispatcher; aliases map friendly names onto canonical ones.

---

## Architecture lessons (what "how Nestimate does this" teaches)

1. **One canonical intermediate, many front doors.** Six input formats all
   reduce to a *list of character-vector transactions* before any math. Detection
   and parsing are isolated, one tiny function each.
2. **The math is two base-R calls.** `crossprod(B)` for counts, `outer()` for
   normalization denominators. No package dependencies, no loops, trivially
   testable.
3. **Normalization is a separable switchboard.** Adding a similarity measure =
   one `else if` branch operating on `(C, freq, N)`. Nothing else changes.
4. **One return type unifies everything.** The `netobject` S3 object is the
   contract; the registry/dispatch pattern keeps `build_network` open for
   extension. This is the tidy-output / clean-delegation principle in practice.
5. **Two engines, one vocabulary.** Transaction co-occurrence and
   sequence-positional co-occurrence are genuinely different algorithms, but they
   share names, return type, and downstream code.

---

## Contrast with this CNA jamovi module

The CNA module (the package in *this* repo) is a **presentation layer**, not an
algorithm:

- Backends (`R/CNA*.b.R`) call `tna::prepare_data()`, then `tna::build_model()`
  / `tna::group_model()`. The co-occurrence matrix is computed **inside `tna`**.
- `R/cna_utils.R` only does (a) data shaping (`cna_prepare_event_data`,
  `cna_prepare_binary_data` — turning event/one-hot tables into the sequence
  form `tna` expects), (b) thin wrappers around `tna` for centrality / community
  / cliques / bootstrap / permutation / comparison, and (c) plotting.

**Implication for this standalone repo:** if you ever want CNA to be
self-contained (not tied to the `tna` runtime), `Nestimate/R/co_occurrence.R` is
the reference port — it is dependency-light base R and supports *more* input
formats and *more* similarity measures than the `tna` build path the module
currently leans on.
