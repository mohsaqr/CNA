# CNA — Co-occurrence Network Analysis (jamovi module)

A jamovi module for **Co-occurrence Network Analysis**, supporting long-format
and one-hot (binary) encoded data, with optional group comparison and automatic
clustering. Includes centrality measures, community detection, clique analysis,
bootstrap validation, and permutation testing.

> Extracted from the [JTNA](https://github.com/mohsaqr/JTNA1.2) workspace as a
> standalone module (2026-06-15).

## Analyses

| Analysis | Data type | Mode |
|----------|-----------|------|
| `CNAEventSingle`   | Event (long format)  | Single / Group |
| `CNAEventCluster`  | Event (long format)  | Clustering |
| `CNABinarySingle`  | Binary (one-hot)     | Single / Group |
| `CNABinaryCluster` | Binary (one-hot)     | Clustering |

Menu group in jamovi: **Co-occurrence Network Analysis**.

## Structure

```
DESCRIPTION          # package metadata (depends on jmvcore, R6, ggplot2, tna)
NAMESPACE
R/                   # jamovi-generated *.h.R (header) + hand-written *.b.R (backend)
  cna_utils.R        # shared co-occurrence / network helpers
jamovi/              # *.a.yaml (analysis), *.u.yaml (UI), *.r.yaml (results), 0000.yaml (manifest)
data/                # bundled example dataset (Regulation_long.csv)
```

## Building

This is a jamovi module. Build/install with the `jmvtools` package:

```r
jmvtools::install()   # build and install into a local jamovi
```

## License

MIT — see [LICENSE](LICENSE).
