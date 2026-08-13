# ProstateAtlasMapR

Map new single-cell datasets onto a prostate single-cell atlas with Seurat's
`FindTransferAnchors()` + `MapQuery()` workflow.

- Choose annotation level: `H1` / `H2` / `H3` (`celltype_manual_*`)
- Project query cells into the atlas **`umap.rpca`** space without moving
  existing reference coordinates
- Ship **code only** in this repo; the lean reference (`.qs`, ~0.9 GB) goes to
  GitHub Releases / Zenodo and is fetched by `get_reference()`

> Spatial deconvolution is planned but **not implemented yet**.

## Installation

```r
# install.packages("remotes")
remotes::install_github("FY10086/ProstateAtlasMapR")
```

## Quick start

### 1. Get the lean reference

**Local file (recommended while developing):**

```r
library(ProstateAtlasMapR)

ref <- get_reference(
  path = "prostate_atlas_reference.qs"
)
```

**After you publish a Release / Zenodo asset:**

```r
options(
  ProstateAtlasMapR.reference_url = "https://github.com/FY10086/ProstateAtlasMapR/releases/download/v0.1.0/prostate_atlas_reference.qs",
  ProstateAtlasMapR.reference_sha256 = "<sha256>"
)
ref <- get_reference()  # download once + cache
```

### 2. Map a query Seurat object

```r
library(Seurat)

# query_assay defaults to "RNA" if present, else DefaultAssay(query).
# For Seurat sketch objects, set query_assay = "sketch".
# Split layers are joined automatically inside map_query().
query <- map_query(
  reference = ref,
  query = query,
  annotation_level = "H2",
  query_assay = "sketch"   # omit if your assay is RNA
)

table(query$predicted.celltype_manual_H2)
DimPlot(query, reduction = "umap.rpca",
        group.by = "predicted.celltype_manual_H2")
```

Optional side-by-side plot:

```r
plot_query_mapping(ref, query, annotation_level = "H2")
```

### 3. (Maintainers) Build the lean reference from the full atlas

```r
ref <- build_reference(
  qc_path = "/path/to/QC_fin.qs",
  save_path = "prostate_atlas_reference.qs"
)
```

Then upload `prostate_atlas_reference.qs` to **GitHub Releases** (or Zenodo),
compute SHA-256, and set the defaults in `R/zzz.R` / `options()`.

**Do not commit the `.qs` file into this git repository** (≈900 MB; GitHub
file limit is 100 MB).

## Package API

| Function | Role |
|----------|------|
| `build_reference()` | Maintainers: slim the full atlas + attach UMAP projection model |
| `get_reference()` | Users: load local path or download/cache remote lean reference |
| `map_query()` | Users: transfer labels + project into `umap.rpca` |
| `plot_query_mapping()` | Optional visualization helper |

## How it works (short)

1. Query expression × reference `rpca` loadings → shared PC/`rpca` space  
2. Anchor-based label transfer in that space  
3. Project into published `umap.rpca` via a `uwot` model attached with
   `n_epochs = 0` so reference coordinates stay unchanged  

## License

MIT
