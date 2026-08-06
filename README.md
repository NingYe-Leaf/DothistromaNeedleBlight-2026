# Ordinal classification of Dothistroma needle blight severity in *Pinus radiata* using UAV multispectral imagery

Analysis code accompanying the manuscript. The pipeline takes UAV multispectral
orthomosaics and field-scored crown polygons, and produces the crown-level
severity models, tables and figures reported in the paper.

## Overview

| Stage | What happens |
|---|---|
| 1 | Load orthomosaics and crown polygons; convert field scores to numeric, flag dead crowns |
| 2 | Convert DN to reflectance, derive a NIR shadow threshold, mask shadow, drop small shadow patches |
| 3 | Extract sunlit pixels within crowns |
| 4 | Compute 82 vegetation indices, aggregate to crowns with six distributional statistics |
| 5 | Split **plots** (not crowns) into train and test |
| 6 | Rank features on training crowns only, then prune collinear features greedily |
| 7 | Quantify selection stability with a plot-level bootstrap |
| 8 | Fit Ordinal Forest and fabOF at four class granularities |
| 9 | Produce figures, site-level metrics and per-trial models |

The train/test split is defined **before** feature selection, so ranking,
collinearity pruning and evaluation all share the same held-out plots. This
avoids selection bias in the reported performance.

## Data availability

This repository contains code only. The datasets used during the current study
are available from the corresponding author on reasonable request.

The pipeline expects two UAV orthomosaics (11-band MicaSense RedEdge-P Dual,
16-bit DN) and two crown polygon shapefiles carrying a field-scored DNB severity
value and a plot identifier. Point `RASTER_SITE1/2` and `CROWNS_SITE1/2` at them
in the configuration block.

## Requirements

R 4.2 or later.

```r
install.packages(c("sf", "terra", "raster", "dplyr", "tidyr", "tibble",
                   "ggplot2", "patchwork", "openxlsx", "ranger", "caret",
                   "psych", "pracma", "LaplacesDemon", "ordinalForest",
                   "scales", "forcats"))

# not on CRAN
devtools::install_github("phibuc/fabOF")
devtools::install_github("ModelOriented/treeshap")
```

`treeshap` is optional. If it is absent the SHAP analysis is skipped with a
message and the rest of the pipeline still runs.

## Running

Edit the paths at the top of `RPBC_DNB_Pipeline.R`, then:

```r
source("RPBC_DNB_Pipeline.R")
```

Everything that needs changing lives in the `CONFIGURATION` block. Nothing below
it should require editing.

### Entry point

`START_FROM` controls where the pipeline begins. Both paths give identical
results; the second just avoids repeating the slowest preprocessing step.

| Value | Starts from |
|---|---|
| `"rasters"` (default) | the orthomosaics and crown shapefiles, caching extracted pixels to `pixel_extraction.csv.gz` |
| `"pixels"` | that cached file, skipping shadow masking and pixel extraction |

Run once with `"rasters"`, then switch to `"pixels"` while iterating on the
modelling.

### Run-control flags

| Flag | Default | Effect when `FALSE` |
|---|---|---|
| `RUN_QC_PLOTS` | `TRUE` | Skips crown/raster overlap and NIR threshold diagnostics |
| `RUN_STABILITY` | `TRUE` | Skips the bootstrap stability analysis (the slowest optional step) |
| `RUN_TRIAL_MODELS` | `TRUE` | Skips the per-trial models, leaving only the combined-trial models |
| `EXPORT_GENETICS` | `TRUE` | Skips the crown feature export used by the companion BLUP analysis |

Set `RUN_STABILITY` and `RUN_TRIAL_MODELS` to `FALSE` for a quick end-to-end
check of the core pipeline.

## Repository layout

```
RPBC_DNB_Pipeline.R                 the whole pipeline
README.md
.gitignore
```

All outputs, including the cached pixel table, are written to `DIR_ANALYSIS`,
which sits outside the repository by default. Add it to `.gitignore` if you
point it inside.

## Runtime

Several hours on a desktop workstation for a full run. The dominant costs are
the ordinal forest hyperparameter grid (60 combinations per granularity) and the
bootstrap stability analysis (`BOOTSTRAP_B = 100` full selection re-runs). All
forests use `num.threads = 0`, so they use every available core.

## Severity class definitions

Class boundaries are defined once, in the configuration block, and reused by the
combined and per-trial pipelines. Changing `CUT_3` changes it everywhere.

| Granularity | Boundaries (% severity) |
|---|---|
| 3-class | 20, 40 |
| 5-class | 20, 40, 60, 80 |
| 10-class | 10, 20, ..., 90 |
| 19-class | the observed 5-point field grades, read from the data |

Tolerance metrics (MAE, W±5%, W±10%) compare the **continuous** field score
against the predicted **class midpoint**, identically at overall, site and
per-trial level. At 3- and 5-class granularity these approach exact class
agreement because the midpoints are widely spaced, so read them alongside QWK
rather than on their own.

## Key outputs

All written to `DIR_ANALYSIS`.

**Tables**

- `results_combined_trials.csv` — headline performance, both algorithms, four granularities
- `results_site_level.csv` — the same metrics split by trial
- `results_three_pipeline_comparison.csv` — combined versus per-trial models
- `feature_ranking_full.csv`, `selected_features.txt`, `selected_features_ranking.csv`
- `feature_selection_stability.csv` — bootstrap selection frequencies
- `feature_overlap_three_pipelines.csv` / `.xlsx`
- `ordfor_<granularity>_predictions.csv`, `_varimp.csv`, `_shap_values.csv`
- `sessionInfo.txt`

**Figures**

- `confusion_matrices_19class.png` — OF and fabOF side by side
- `spectral_profile_combined.png`, `spectral_profile_per_trial.png`
- `class_distribution_train_test.png`
- `split_map_s1.png`, `split_map_s2.png`
- `feature_selection_stability.png`
- `shap_importance_OF_19class.png`, `shap_importance_fabOF_19class.png`

fabOF SHAP values describe the latent continuous score the model thresholds
internally, not the final ordinal class. They are supplementary and labelled as
such in the figure.

## Reproducibility

All random operations are seeded from `SEED` in the configuration block. A
session snapshot is written to `sessionInfo.txt` on every run.

