# ProstateAtlasMapR 0.0.0.9000

* Core: `build_reference()`, `get_reference()`,
  `map_query()`, `run_rctd()`.
* RM plots: `plot_rm_score_hist()`, `plot_rm_score()`,
  `plot_rm_prediction()` (Unknown via `score_threshold`).
* RCTD plots: `plot_rctd_first_type()`, `plot_rctd_composition()`,
  `plot_rctd_weights()`, `plot_rctd_weight_vs_marker()`, `plot_rctd_pie()`,
  `plot_rctd_qc()`. Spatial panels use `coord_fixed()` (true tissue aspect,
  matching `10-spatial.ipynb`); pie radii use median NN pitch.
* Removed `plot_query_mapping()`.
* Lean reference data is distributed as a GitHub Release asset
  (`prostate_atlas_reference.qs`, v0.1.0).
* `get_reference()` defaults to GitHub Release
  `v0.1.0` (`prostate_atlas_reference.qs`) with baked-in SHA-256; downloads
  follow redirects, use an atomic `.part` file, and cache under
  `tools::R_user_dir("ProstateAtlasMapR", "cache")`. Removed redundant
  `download_reference()` alias.