# Ordinal classification of Dothistroma needle blight severity in *Pinus radiata* using UAV multispectral imagery

Analysis code accompanying the manuscript. The main pipeline takes UAV
multispectral orthomosaics and field-scored crown polygons, and produces the
crown-level severity models, tables and figures reported in the paper. Two
further scripts carry out the validation analyses: one tests how far the
reported performance survives spatial separation between training and test
crowns, the other tests how far it depends on where the crown boundaries were
drawn.

```
RPBC_DNB_Pipeline.R             main pipeline: preprocessing to models
spatial_validation.R            buffered and spatially blocked validation
crown_boundary_sensitivity.R    boundary perturbation analysis
README.md
.gitignore
```

Run the main pipeline first. Both validation scripts read artefacts it writes.

## Main pipeline

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

## Validation analyses

Both scripts are self-contained apart from the artefacts listed under
[Inputs](#inputs-the-validation-scripts-expect). Neither re-tunes
hyperparameters or repeats feature selection, so each isolates one source of
variation.

### `spatial_validation.R`

A plot-level split prevents trees from the same plot appearing on both sides,
but it does not make training and test crowns spatially independent: in these
trials the median test crown lies about 4.6 m from the nearest training crown.
The script increases that separation in two ways and reports what happens.

| Check | Separation achieved by | Model refitted |
|---|---|---|
| Buffered evaluation | discarding test crowns close to training data | no, saved predictions are re-scored |
| Blocked cross-validation | holding out whole blocks of plots, then deleting training crowns within an exclusion radius | yes, fabOF only |
| Random-removal control | removing the same number of training crowns at random | yes, fabOF only |

The control exists because the exclusion radius removes the nearest, and so the
most informative, training crowns. Without it the resulting loss cannot be
separated from the cost of a smaller training set.

Set `RUN_BLOCKED_CV <- FALSE` for the buffered check alone, which needs no
model fitting and finishes in under a minute.

### `crown_boundary_sensitivity.R`

The models draw predominantly on within-crown distributional statistics, which
depend on which pixels fall inside the boundary. Every crown is eroded and
dilated by one and two pixels at each trial's native ground sampling distance,
the statistics are recomputed, and the existing model is applied without
refitting. The script reports the realised change in pixel count, the
correlation of each selected feature with its original value, and the resulting
movement in predictions.

This tests sensitivity to boundary **position**. It does not reproduce the
structured, shape-dependent errors that automated segmentation would introduce,
and feature selection is not repeated, so it says nothing about whether the
same features would be chosen again.

## Which script produces what

| Manuscript item | Script |
|---|---|
| Table 4, Table 5 | `RPBC_DNB_Pipeline.R` |
| Figures 3 to 8 | `RPBC_DNB_Pipeline.R` |
| Appendix Tables 1, 2, 3, 5 | `RPBC_DNB_Pipeline.R` |
| Appendix Figures 2 to 6 | `RPBC_DNB_Pipeline.R` |
| Appendix Table 4, Sections 3.2 and 4 spatial figures | `spatial_validation.R` |
| Appendix Table 6, Sections 2.6 and 4 boundary figures | `crown_boundary_sensitivity.R` |

## Data availability

This repository contains code only. The datasets used during the current study
are available from the corresponding author on reasonable request.

The pipeline expects two UAV orthomosaics (11-band MicaSense RedEdge-P Dual,
16-bit DN) and two crown polygon shapefiles carrying a field-scored DNB severity
value and a plot identifier. Point `RASTER_SITE1/2` and `CROWNS_SITE1/2` at them
in the configuration block.

### Inputs the validation scripts expect

Both read these from `DIR_ANALYSIS`.

| Artefact | Filename | Used by |
|---|---|---|
| Selected-feature list | `selected_features.txt` | both |
| Crown feature table | `crown_features_for_genetic_analysis.csv` | both |
| 19-class OF test predictions | `ordfor_19-class_predictions.csv` | both |
| 19-class fabOF test predictions | `fabof_19-class_predictions.csv` | `spatial_validation.R`, optional |
| Reported performance table | `results_combined_trials.csv` | `spatial_validation.R`, to verify itself |
| Shadow-masked reflectance rasters | `Trial1_reflectance_masked.tif`, `Trial2_…` | `crown_boundary_sensitivity.R` |

Output names and column headings have changed between pipeline versions, so
each input is resolved from a short list of candidates given in the
`CONFIGURATION` block of each script, and the first that exists is used. If
yours differ, edit the relevant `NAMES_*` vector rather than hunting through the
code. A missing required input stops the script and names what it looked for.

The fabOF predictions are optional, because `EXPORT_FABOF_PREDICTIONS` can be
turned off in the main pipeline. Without them `spatial_validation.R` says so and
reports the OF arm of the buffered check only, which loses the fabOF row of
Appendix Table 4.

Each script also re-scores the unperturbed case first and checks it against the
reported figures, printing `PASS` or `FAIL`. A failure there means an input has
changed and nothing downstream should be trusted.

## Requirements

R 4.2 or later.

```r
install.packages(c("sf", "terra", "raster", "dplyr", "tidyr", "tibble",
                   "rlang", "ggplot2", "patchwork", "openxlsx", "ranger",
                   "caret", "psych", "pracma", "LaplacesDemon",
                   "ordinalForest", "scales", "forcats"))

# not on CRAN
devtools::install_github("phibuc/fabOF")
devtools::install_github("ModelOriented/treeshap")
```

`treeshap` is optional. If it is absent the SHAP analysis is skipped with a
message and the rest of the pipeline still runs.

## Running

Edit the paths at the top of each script, then:

```r
source("RPBC_DNB_Pipeline.R")           # first; writes everything below depends on
source("spatial_validation.R")
source("crown_boundary_sensitivity.R")
```

Everything that needs changing lives in the `CONFIGURATION` block of each
script. Nothing below it should require editing.

### Entry point

`START_FROM` controls where the main pipeline begins. Both paths give identical
results; the second just avoids repeating the slowest preprocessing step.

| Value | Starts from |
|---|---|
| `"rasters"` (default) | the orthomosaics and crown shapefiles, caching extracted pixels to `pixel_extraction.csv.gz` |
| `"pixels"` | that cached file, skipping shadow masking and pixel extraction |

Run once with `"rasters"`, then switch to `"pixels"` while iterating on the
modelling.

### Run-control flags

| Script | Flag | Default | Effect when `FALSE` |
|---|---|---|---|
| main | `RUN_QC_PLOTS` | `TRUE` | Skips crown/raster overlap and NIR threshold diagnostics |
| main | `RUN_STABILITY` | `TRUE` | Skips the bootstrap stability analysis (the slowest optional step) |
| main | `RUN_TRIAL_MODELS` | `TRUE` | Skips the per-trial models, leaving only the combined-trial models |
| main | `EXPORT_GENETICS` | `TRUE` | Skips the crown feature export used by the companion BLUP analysis |
| main | `EXPORT_FABOF_PREDICTIONS` | `TRUE` | Skips saving per-crown fabOF predictions, which `spatial_validation.R` needs for its fabOF arm |
| spatial | `RUN_BLOCKED_CV` | `TRUE` | Leaves only the buffered check, which refits nothing |
| spatial | `RUN_CONTROL` | `TRUE` | Skips the random-removal arm, losing the decomposition |
| boundary | `RUN_PREDICT` | `TRUE` | Reports feature stability only, without refitting or scoring |

Set `RUN_STABILITY` and `RUN_TRIAL_MODELS` to `FALSE` for a quick end-to-end
check of the core pipeline.

## Runtime

| Script | Typical | Dominated by |
|---|---|---|
| `RPBC_DNB_Pipeline.R` | several hours | the ordinal forest hyperparameter grid (60 combinations per granularity) and the bootstrap stability analysis (`BOOTSTRAP_B = 100` full selection re-runs) |
| `spatial_validation.R` | ~20 minutes | 25 fabOF fits; under a minute with `RUN_BLOCKED_CV <- FALSE` |
| `crown_boundary_sensitivity.R` | ~1 hour | pixel extraction, five times over |

All forests use `num.threads = 0`, so they use every available core. The
boundary script re-extracts roughly 23 million pixel rows per perturbation; if
the rasters sit on a network share, copy them to a local disk first, as that
will dominate everything else.

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

Quadratic weighted kappa in the validation scripts is computed on the full
5 to 95% scale rather than on the grades present in each subset. A grade missing
from a subset would otherwise shorten the ordinal distance between its
neighbours and inflate the statistic, which matters because the subsets shrink
as separation increases.

## Key outputs

All written to `DIR_ANALYSIS`.

**Main pipeline, tables**

- `results_combined_trials.csv` — headline performance, both algorithms, four granularities
- `results_site_level.csv` — the same metrics split by trial
- `results_three_pipeline_comparison.csv` — combined versus per-trial models
- `feature_ranking_full.csv`, `selected_features.txt`, `selected_features_ranking.csv`
- `crown_features_for_genetic_analysis.csv` — crown feature table, also read by both validation scripts
- `ordfor_<granularity>_predictions.csv`, `fabof_<granularity>_predictions.csv`, `_varimp.csv`, `_shap_values.csv`
- `feature_selection_stability.csv` — bootstrap selection frequencies
- `feature_overlap_three_pipelines.csv` / `.xlsx`
- `sessionInfo.txt`

**Main pipeline, figures**

- `confusion_matrices_19class.png` — OF and fabOF side by side
- `spectral_profile_combined.png`, `spectral_profile_per_trial.png`
- `class_distribution_train_test.png`
- `split_map_s1.png`, `split_map_s2.png`
- `feature_selection_stability.png`
- `shap_importance_OF_19class.png`, `shap_importance_fabOF_19class.png`

**Spatial validation**

- `spatial_buffered_evaluation.csv` — performance against buffer distance
- `spatial_blocked_cv_folds.csv`, `spatial_blocked_cv_summary.csv`
- `spatial_decomposition.csv` — the loss split into training-set size and spatial dependence

**Boundary sensitivity**

- `boundary_perturbation.csv` — realised change in crown pixel count per trial
- `boundary_feature_stability.csv`, `boundary_feature_stability_summary.csv`
- `boundary_performance.csv`, `boundary_agreement.csv`

fabOF SHAP values describe the latent continuous score the model thresholds
internally, not the final ordinal class. They are supplementary and labelled as
such in the figure.

## Reproducibility

All random operations are seeded from `SEED` in the configuration block of each
script. A session snapshot is written to `sessionInfo.txt` on every run of the
main pipeline.

Multithreaded `ranger` is not bit-reproducible across differing thread counts,
so a refitted model may differ from the original on a small number of crowns.
Both validation scripts therefore verify the unperturbed case against the
reported figures rather than assuming it.
