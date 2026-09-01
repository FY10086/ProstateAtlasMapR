# ProstateAtlasMapR

Map new single-cell datasets onto a prostate atlas, and run spatial
deconvolution (RCTD) with the same H1/H2/H3 labels.

- **Reference mapping**: Seurat `FindTransferAnchors()` + `MapQuery()`,
  projecting query cells into atlas **`umap.rpca`**
- **RCTD**: `spacexr` deconvolution on Visium / spot Seurat objects
- **Plot helpers** for prediction scores and RCTD weights / composition
- Code only in this repo; lean reference (`.qs`, ~0.9 GB) via
  [GitHub Releases](https://github.com/FY10086/ProstateAtlasMapR/releases)

## Installation

```r
# install.packages(c("remotes", "openssl"))  # openssl: checksum for get_reference()
remotes::install_github("FY10086/ProstateAtlasMapR")
# for RCTD:
# install.packages("spacexr")
# for pie charts:
# install.packages("scatterpie")
# for H3 side-legend layouts:
# install.packages("cowplot")
```

## Quick start

### 1. Get the lean reference

```r
library(ProstateAtlasMapR)

## First use downloads ~0.9 GB (SHA-256 checked), then reuses cache:
ref <- get_reference()
## Or use a local copy:
# ref <- get_reference(path = "prostate_atlas_reference.qs")
```

Release asset:
https://github.com/FY10086/ProstateAtlasMapR/releases/download/v0.1.0/prostate_atlas_reference.qs  
SHA-256: `c83d94b5326aad495cb854f98fbb31ca10a0b1143eaf45bab055f368a8071d5f`

### 2. Reference mapping (scRNA)

```r
query <- map_query(
  reference = ref,
  query = query,
  annotation_level = "H2"
  # query_assay = "sketch"   # if needed
)

plot_rm_score_hist(query, annotation_level = "H2")
plot_rm_score(query, annotation_level = "H2")
plot_rm_prediction(query, annotation_level = "H2", score_threshold = 0.6)
plot_rm_prediction(query, annotation_level = "H3", score_threshold = 0.6)
# Sets options(repr.plot.*) automatically: H1/H2 = 12x10, H3 = 22x10
```

### 3. Spatial RCTD

```r
spatial <- run_rctd(
  reference = ref,
  spatial = spatial,
  annotation_level = "H2",
  coords = "spatial",          # DimReduc name, or omit for GetTissueCoordinates()
  doublet_mode = "full"        # or "doublet" / "multi"
)

# annotation_level = NULL → last run_rctd() level (else "H2")
plot_rctd_first_type(spatial, annotation_level = "H2")
plot_rctd_first_type(spatial, annotation_level = "H3")  # single-col full-name legend
plot_rctd_composition(spatial, annotation_level = "H2")
plot_rctd_weights(spatial, annotation_level = "H2", n_types = 11)
plot_rctd_weight_vs_marker(spatial, annotation_level = "H2", n_types = 6)
plot_rctd_weight_vs_marker(spatial, annotation_level = "H3", n_types = 6)  # H3 gene-prefix markers
plot_rctd_pie(spatial, annotation_level = "H2", top_n = 12)   # needs scatterpie
plot_rctd_pie(spatial, annotation_level = "H3", top_n = 12)   # single-col full-name legend
plot_rctd_qc(spatial, annotation_level = "H2")
```

### 4. (Maintainers) Build the lean reference

```r
ref <- build_reference(
  qc_path = "/path/to/QC_fin.qs",
  save_path = "prostate_atlas_reference.qs"
)
```

## Main functions

| Function | Role |
|----------|------|
| `build_reference()` | Slim atlas → distributable `patlas_reference` |
| `get_reference()` | Download+cache Release `.qs` (or load local `path`) |
| `map_query()` | Transfer H1/H2/H3 labels + project to `umap.rpca` |
| `run_rctd()` | spacexr RCTD; writes level-prefixed weights + meta onto spatial Seurat |
| `plot_rm_score_hist()` | Prediction-score histogram |
| `plot_rm_score()` | Score FeaturePlot on `umap.rpca` |
| `plot_rm_prediction()` | Predicted labels (H3: cowplot UMAP \| single-col full-name legend, 22×10) |
| `plot_rctd_first_type()` | Spatial dominant cell type (H3: large map \| single-col full-name legend) |
| `plot_rctd_composition()` | Mean weight + dominant-count summary |
| `plot_rctd_weights()` | Faceted weight maps |
| `plot_rctd_weight_vs_marker()` | Weight vs marker gene (H3 resolves marker from gene prefix) |
| `plot_rctd_pie()` | Per-spot composition pies (H3: same single-col side legend as first_type) |
| `plot_rctd_qc()` | Max weight + Shannon entropy maps |

## Plot sizes (`options(repr.plot.*)`)

Set automatically by each plot helper (Jupyter / IRkernel):

| Plot | H1 / H2 | H3 |
|------|---------|-----|
| `plot_rm_prediction()` | 12 × 10 | 22 × 10 |
| `plot_rctd_first_type()` | 12 × 10 | ~map + legend (often ≥ 17 × 13) |
| `plot_rctd_pie()` | 13 × 11 | ~map + legend (tall single-col) |
| `plot_rctd_composition()` | 13 × 7 | 16 × 8 |
| `plot_rctd_weights()` | width 13; height ≈ 4.3 × rows + 0.8 | same |
| `plot_rctd_weight_vs_marker()` | width 13; height ≈ 6.2 × n_types + 0.6 | same |
| `plot_rctd_qc()` | 11 × 5.5 | same |

**H3 spatial legends** (`first_type`, `pie`): always **one column** of full subtype names (no 2-col split); need **`cowplot`**.

## Notes

- Query prep for RM: `NormalizeData` is enough; split layers are joined inside `map_query()`.
- Reference download: `get_reference()` defaults to Release
  [`v0.1.0`](https://github.com/FY10086/ProstateAtlasMapR/releases/tag/v0.1.0);
  install **openssl** for SHA-256 checks. Override with
  `options(ProstateAtlasMapR.reference_url=..., ProstateAtlasMapR.reference_sha256=...)`
  or `get_reference(path=...)`.
- RCTD: pass `coords` as a reduction name (e.g. `"spatial"` / `"physical"`) when you store array coords as a DimReduc.
- `run_rctd()` results (per level, e.g. `H2`):
  - `@misc$rctd_weights_<level>` (and unprefixed `@misc$rctd_weights` for the last run)
  - meta: `rctd_<level>_first_type`, `rctd_<level>_spot_class`, …
  - optional assay `"RCTD_<level>"` (and `"RCTD"` for the last run)
  - `@misc$rctd_annotation_level` records the last run
- Plot helpers look up `rctd_weights_<level>` first; re-run `run_rctd(..., annotation_level = ...)` if that level is missing.
- `plot_rctd_weight_vs_marker()` skips RCTD assays for expression, auto-`NormalizeData` if needed, and for H3 names like `POSTN+FAP+ …` uses the first gene present in the assay.
