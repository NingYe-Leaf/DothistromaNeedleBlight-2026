# ==============================================================================
# Ordinal classification of Dothistroma needle blight severity in Pinus radiata
# using UAV multispectral imagery
#
# Analysis pipeline: shadow masking -> vegetation indices -> crown-level features
#                    -> leakage-free feature selection -> ordinal forest models
#
# Sensor    MicaSense RedEdge-P Dual (11 bands), DJI Matrice 350 RTK
# Models    Ordinal Forest (OF) and frequency-adjusted borders OF (fabOF)
# Outputs   All tables and figures reported in the manuscript
#
# Reproducibility notes
#   - Edit the CONFIGURATION block below; nothing else needs changing.
#   - The train/test split is defined at plot level BEFORE feature selection,
#     so ranking, collinearity pruning and evaluation share the same held-out
#     plots (avoids selection bias, Ambroise & McLachlan 2002).
#   - All random operations are seeded via SEED.
# ==============================================================================


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

## ---- Paths -----------------------------------------------------------------
DIR_BASE     <- "J:/RPBC_Dothi_Analysis"          # <- edit this
DIR_REF      <- file.path(DIR_BASE, "Reference")
DIR_MS       <- file.path(DIR_BASE, "Multispectral")
DIR_SHADOW   <- file.path(DIR_BASE, "Shadow_Masks")
DIR_ANALYSIS <- file.path(DIR_BASE, "Analysis")

## Raw imagery and crown polygons. Not distributed with this repository; see the
## data availability statement in the README.
RASTER_SITE1 <- file.path(DIR_MS, "Mutispectral_16122025_SunNBoard_georeference2.tif")
RASTER_SITE2 <- file.path(DIR_MS, "Mutispectral_15122025_SunNBoard.tif")
CROWNS_SITE1 <- file.path(DIR_REF, "Crown_Health_Reference_2026.shp")
CROWNS_SITE2 <- file.path(DIR_REF, "Crown_Health_Reference_2026_site2_clean.shp")

## ---- Entry point -----------------------------------------------------------
## "rasters" runs the full chain from the orthomosaics and crown polygons, and
##           caches the extracted pixels to PIXEL_TABLE on the way through.
## "pixels"  skips straight to index calculation by reading that cache, which
##           saves the masking and extraction step on repeat runs.
START_FROM   <- "rasters"
PIXEL_TABLE  <- file.path(DIR_ANALYSIS, "pixel_extraction.csv.gz")

## Columns the cached pixel table carries. b5 (panchromatic) is unused by the
## index formulas but is kept so the table matches the raw sensor. Dead crowns
## and shadowed pixels are already excluded, so no status flag is needed.
PIXEL_TABLE_COLS <- c(paste0("b", 1:11),
                      "Label_12", "DOTHI___IN", "plot_id", "site")

## ---- Run control -----------------------------------------------------------
SEED              <- 42
RUN_QC_PLOTS      <- TRUE    # crown/raster overlap and NIR threshold diagnostics
RUN_STABILITY     <- TRUE    # bootstrap feature-selection stability
RUN_TRIAL_MODELS  <- TRUE    # per-trial models in addition to the combined model
EXPORT_GENETICS   <- TRUE    # crown feature table for the companion BLUP analysis
EXPORT_FABOF_PREDICTIONS <- TRUE   # per-crown fabOF predictions, read by spatial_validation.R
CM_SHOW_50PCT_REF <- FALSE   # dotted reference line at 50% severity on the CM figure

## ---- Preprocessing ---------------------------------------------------------
DN_SCALE       <- 65535        # 16-bit DN -> reflectance
NIR_BAND       <- 11           # 842 nm
MIN_SHADOW_M2  <- 0.02         # minimum shadow patch area retained (m2)
NA_THRESHOLD   <- 0.10         # drop features with >10% NA among training crowns

## ---- Feature selection -----------------------------------------------------
COR_CUTOFF     <- 0.7          # greedy collinearity prune, |Pearson r|
RF_TREES_UNIV  <- 200          # trees per univariate ranking forest
BOOTSTRAP_B    <- 100          # plot-level bootstrap replicates for stability
BOOTSTRAP_RF   <- 100          # trees per univariate forest inside the bootstrap

## ---- Train/test split ------------------------------------------------------
TRAIN_FRACTION <- 0.7          # proportion of PLOTS assigned to training

## ---- Ordinal forest hyperparameters ----------------------------------------
OF_TUNE  <- list(nsets = 1000, nbest =  5, ntreeperdiv = 100, ntreefinal = 1000)
OF_FINAL <- list(nsets = 1000, nbest = 10, ntreeperdiv = 100, ntreefinal = 5000)
OF_MTRY_GRID     <- c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
OF_NODESIZE_GRID <- c(5, 10, 20, 50, 100, 150)

## ---- fabOF hyperparameters -------------------------------------------------
FABOF_NTREES <- 5000
FABOF_MTRY   <- 10

## ---- Severity class definitions --------------------------------------------
## Defined ONCE here and reused by every model, so the combined and per-trial
## pipelines can never drift apart.
CUT_3  <- c(20, 40)                                     # 3-class breakpoints (%)
CUT_5  <- c(20, 40, 60, 80)                             # 5-class breakpoints (%)
CUT_10 <- seq(10, 90, by = 10)                          # 10-class breakpoints (%)
## 19-class levels are the observed 5-point field grades and are read from data.

MID_3  <- c("1_Low" = 10, "2_Moderate" = 30, "3_High" = 67.5)
MID_5  <- c("1_Low" = 10, "2_Moderate" = 30, "3_High" = 50,
            "4_VeryHigh" = 70, "5_Severe" = 90)
MID_10 <- setNames(seq(5, 95, by = 10),
                   c("1_0to10","2_10to20","3_20to30","4_30to40","5_40to50",
                     "6_50to60","7_60to70","8_70to80","9_80to90","10_90plus"))


# ==============================================================================
# 2. PACKAGES
# ==============================================================================
required <- c("sf", "terra", "raster", "dplyr", "tidyr", "tibble", "ggplot2",
              "patchwork", "openxlsx", "ranger", "caret", "psych", "pracma",
              "LaplacesDemon", "ordinalForest", "scales", "forcats")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")

suppressPackageStartupMessages({
  library(sf); library(terra); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(openxlsx); library(ranger)
  library(caret); library(psych); library(ordinalForest)
})

# fabOF is not on CRAN
if (!requireNamespace("fabOF", quietly = TRUE))
  stop("fabOF not installed. Run: devtools::install_github('phibuc/fabOF')")
library(fabOF)

HAS_TREESHAP <- requireNamespace("treeshap", quietly = TRUE)
if (!HAS_TREESHAP)
  message("treeshap not installed - SHAP analysis will be skipped.")

for (d in c(DIR_SHADOW, DIR_ANALYSIS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 3. UTILITIES
# ==============================================================================

clean_names_fn <- function(x) gsub("[^a-zA-Z0-9_]", "_", x)

out_path <- function(...) file.path(DIR_ANALYSIS, ...)

save_fig <- function(plot_obj, filename, width = 190, height = 120, dpi = 300) {
  ggsave(out_path(filename), plot_obj, width = width, height = height,
         units = "mm", dpi = dpi, limitsize = FALSE)
  cat(sprintf("  Figure saved: %s\n", filename))
}

save_table <- function(df, filename) {
  write.csv(df, out_path(filename), row.names = FALSE)
  cat(sprintf("  Table saved: %s\n", filename))
}

## Coefficient of variation, guarded against near-zero means and outliers
safe_cv <- function(x) {
  v <- x[!is.na(x) & is.finite(x)]
  if (length(v) < 2) return(NA_real_)
  m <- mean(v); s <- sd(v)
  if (is.na(m) || is.na(s) || m == 0) return(NA_real_)
  cv <- s / abs(m)
  if (cv > 10) NA_real_ else cv
}

## Quadratic weighted kappa; returns NA rather than halting the run
compute_qwk <- function(true_cat, pred_cat, class_levels, label = "") {
  unexpected <- setdiff(as.character(pred_cat), class_levels)
  if (length(unexpected))
    cat(sprintf("  Warning [%s]: unexpected predicted values: %s\n",
                label, paste(unexpected, collapse = ", ")))
  t_num <- as.numeric(factor(as.character(true_cat), levels = class_levels))
  p_num <- as.numeric(factor(as.character(pred_cat), levels = class_levels))
  tryCatch(psych::cohen.kappa(cbind(t_num, p_num))$weighted.kappa,
           error = function(e) NA_real_)
}

## Classification and tolerance metrics shared by OF and fabOF
compute_metrics <- function(cm, pred_numeric, true_numeric, naive_acc, step = 5) {
  list(
    accuracy     = unname(cm$overall["Accuracy"]),
    balanced_acc = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
    kappa        = unname(cm$overall["Kappa"]),
    naive        = naive_acc,
    vs_naive     = unname(cm$overall["Accuracy"]) - naive_acc,
    mae          = mean(abs(pred_numeric - true_numeric), na.rm = TRUE),
    rmse         = sqrt(mean((pred_numeric - true_numeric)^2, na.rm = TRUE)),
    r2           = 1 - sum((true_numeric - pred_numeric)^2, na.rm = TRUE) /
                       sum((true_numeric - mean(true_numeric))^2),
    within_1     = mean(abs(pred_numeric - true_numeric) <= 1 * step, na.rm = TRUE),
    within_2     = mean(abs(pred_numeric - true_numeric) <= 2 * step, na.rm = TRUE)
  )
}


# ==============================================================================
# 4. SEVERITY CLASS CONSTRUCTION
# One definition per granularity, used by every model in the study.
# ==============================================================================

cut_to_classes <- function(x, breakpoints, labels) {
  factor(cut(x, breaks = c(-Inf, breakpoints, Inf), labels = labels,
             right = FALSE, include.lowest = TRUE),
         levels = labels, ordered = TRUE)
}

## Returns the four granularity specifications for a given dataset.
## `levels_19` is data-driven (the observed 5-point field grades).
make_granularities <- function(y_numeric) {
  levels_19 <- sort(unique(y_numeric))
  list(
    `3-class` = list(
      fn   = function(d) cut_to_classes(d, CUT_3, names(MID_3)),
      mids = MID_3),
    `5-class` = list(
      fn   = function(d) cut_to_classes(d, CUT_5, names(MID_5)),
      mids = MID_5),
    `10-class` = list(
      fn   = function(d) cut_to_classes(d, CUT_10, names(MID_10)),
      mids = MID_10),
    `19-class` = list(
      fn   = function(d) factor(d, levels = levels_19, ordered = TRUE),
      mids = setNames(levels_19, as.character(levels_19)))
  )
}

GRAN_ORDER <- c("3-class", "5-class", "10-class", "19-class")


# ==============================================================================
# 5. PREPROCESSING FUNCTIONS
# ==============================================================================

## Field scores are stored as text and include "dead" crowns; convert to numeric
## and flag dead crowns so they can be excluded from modelling.
fix_health <- function(crowns) {
  crowns %>%
    mutate(raw     = trimws(as.character(DOTHI___IN)),
           is_dead = !is.na(raw) & tolower(raw) == "dead",
           DOTHI___IN = ifelse(is_dead, NA_real_,
                               suppressWarnings(as.numeric(raw)))) %>%
    select(-raw)
}

otsu_threshold <- function(values, n_bins = 256) {
  breaks <- seq(min(values), max(values), length.out = n_bins + 1)
  probs  <- hist(values, breaks = breaks, plot = FALSE)$counts
  probs  <- probs / sum(probs)
  best_t <- 0; best_v <- 0
  for (i in seq_len(n_bins - 1)) {
    w0 <- sum(probs[1:i]); w1 <- sum(probs[(i + 1):n_bins])
    if (w0 == 0 || w1 == 0) next
    mu0 <- sum(probs[1:i] * breaks[1:i]) / w0
    mu1 <- sum(probs[(i + 1):n_bins] * breaks[(i + 1):n_bins]) / w1
    v <- w0 * w1 * (mu0 - mu1)^2
    if (v > best_v) { best_v <- v; best_t <- breaks[i + 1] }
  }
  best_t
}

## Zero crossings of the density gradient, used to locate the shadow/sunlit
## trough in a bimodal NIR histogram.
find_local_min <- function(values) {
  n <- length(values)
  res <- list(); k <- 1
  for (i in 2:n) {
    if (values[i - 1] < 0 && values[i] >= 0 && i > 15) {
      dp <- values[i:min(i + 15, n)]
      dn <- values[max(i - 15, 1):i]
      res[[k]] <- data.frame(Index = i - 1,
                             definition = sum(abs(dn), abs(dp)),
                             pos_def = sum(dp[dp > 0]))
      k <- k + 1
    }
  }
  if (!length(res)) return(data.frame(Index = numeric(0), definition = numeric(0),
                                      pos_def = numeric(0)))
  do.call(rbind, res)
}

## Local-minimum threshold when the NIR distribution is multimodal, Otsu otherwise.
find_nir_threshold <- function(raster_obj, crowns_vect, site_label) {
  nir <- as.vector(values(crop(raster_obj[[NIR_BAND]], crowns_vect)))
  nir <- nir[!is.na(nir) & nir > 0 & nir < 0.99]
  dens <- density(nir, n = 1024)
  method <- "Otsu"; thr <- otsu_threshold(nir)
  if (LaplacesDemon::is.multimodal(nir)) {
    z <- find_local_min(pracma::gradient(dens$y))
    z <- z[z$pos_def > 0, , drop = FALSE]
    if (nrow(z)) {
      cand <- dens$x[z$Index[which.max(z$definition)]]
      if (cand > quantile(nir, 0.01) && cand < quantile(nir, 0.99)) {
        thr <- cand; method <- "LocalMin"
      }
    }
  }
  cat(sprintf("  [%s] NIR threshold = %.4f (%s)\n", site_label, thr, method))
  list(threshold = thr, method = method, density = dens, nir = nir)
}

plot_nir_distribution <- function(nir_res, site_label) {
  df <- data.frame(x = nir_res$density$x, y = nir_res$density$y)
  thr <- nir_res$threshold
  p <- ggplot(df, aes(x, y)) +
    geom_area(data = subset(df, x <= thr), fill = "#2166ac", alpha = 0.35) +
    geom_area(data = subset(df, x >= thr), fill = "#d6604d", alpha = 0.35) +
    geom_line(colour = "grey20", linewidth = 0.8) +
    geom_vline(xintercept = thr, linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = thr + 0.01, y = max(df$y) * 0.97, hjust = 0, vjust = 1,
             size = 3.5, colour = "grey15",
             label = sprintf("Threshold = %.4f\n(%s)", thr, nir_res$method)) +
    labs(x = "NIR reflectance", y = "Density",
         title = sprintf("NIR reflectance distribution - %s", site_label),
         subtitle = sprintf("n = %s pixels", format(length(nir_res$nir), big.mark = ","))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          panel.grid.minor = element_blank())
  ggsave(file.path(DIR_SHADOW, sprintf("%s_NIR_distribution.png", site_label)),
         p, width = 180, height = 110, units = "mm", dpi = 300)
}

## Threshold NIR, remove shadow patches smaller than MIN_SHADOW_M2, and return
## both the raw and the size-filtered mask.
build_shadow_mask <- function(raster_obj, crowns_vect, threshold) {
  sm <- crop(raster_obj[[NIR_BAND]], crowns_vect)
  sm[sm >  threshold] <- NA
  sm[sm <= threshold] <- 1
  clumps <- rast(raster::clump(raster::raster(sm), directions = 4))
  px_area <- (res(raster_obj)[1] * 111320) * (res(raster_obj)[2] * 111320)
  min_px  <- MIN_SHADOW_M2 / px_area
  keep    <- data.frame(freq(clumps))
  keep    <- keep$value[keep$count > min_px]
  mask_clean <- clumps %in% keep
  mask_clean[mask_clean == 0] <- NA
  list(raw   = sm,
       clean = terra::project(mask_clean, raster_obj, method = "near", align = TRUE))
}

extract_pixels <- function(raster_obj, crowns, site_label) {
  ed <- terra::extract(raster_obj, vect(crowns), ID = TRUE)
  colnames(ed)[2:12] <- paste0("b", 1:11)
  ed$Label_12   <- crowns$Label_12[ed$ID]
  ed$DOTHI___IN <- crowns$DOTHI___IN[ed$ID]
  ed$is_dead    <- crowns$is_dead[ed$ID]
  ed$plot_id    <- crowns$plot_id[ed$ID]
  ed$site       <- site_label
  ed
}


# ==============================================================================
# 6. VEGETATION INDICES
# MicaSense RedEdge-P Dual band centres:
#   b1  Coastal blue 444 nm    b2  Blue      475 nm
#   b3  Green 1      531 nm    b4  Green 2   560 nm
#   b5  Panchromatic (excluded from index calculation)
#   b6  Red 1        650 nm    b7  Red 2     668 nm
#   b8  Red edge 1   705 nm    b9  Red edge 2 717 nm
#   b10 Red edge 3   740 nm    b11 NIR       842 nm
# ==============================================================================

RAW_BANDS <- paste0("b", c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11))

calc_vegetation_indices <- function(px) {
  px %>%
    mutate(
      NDVI       = (b11 - b7) / (b11 + b7),
      NDRE       = (b11 - b9) / (b11 + b9),
      SAVI       = 1.5 * (b11 - b7) / (b11 + b7 + 0.5),
      DVI        = b11 - b7,
      EVI        = 2.5 * (b11 - b7) / (b11 + 6 * b7 - 7.5 * b2 + 1),
      EVI2       = 2.5 * (b11 - b7) / (b11 + 2.4 * b7 + 1),
      RDVI       = (b11 - b7) / sqrt(b11 + b7),
      MSR        = (b11 / b7 - 1) / (sqrt(b11 / b7) + 1),
      MSAVI      = (2 * b11 + 1 - sqrt((2 * b11 + 1)^2 - 8 * (b11 - b7))) / 2,
      SR         = b11 / b7,
      MTVI1      = 1.2 * (1.2 * (b11 - b4) - 2.5 * (b7 - b4)),
      MTVI2      = (1.5 * (1.2 * (b11 - b4) - 2.5 * (b7 - b4))) /
                     sqrt((2 * b11 + 1)^2 - (6 * b11 - 5 * sqrt(b7)) - 0.5),
      MTVI2_den  = sqrt((2 * b11 + 1)^2 - (6 * b11 - 5 * sqrt(b7)) - 0.5),
      MCARI1     = 1.2 * (2.5 * (b11 - b7) - 1.3 * (b11 - b4)),
      SR528_574  = b3 / b4,
      SR800_635  = b11 / b6,
      NDVI2      = (b10 - b8) / (b10 + b8),
      NDRE3      = (b11 - b10) / (b11 + b10),
      NDRE_705   = (b11 - b8) / (b11 + b8),
      NDVI_RE    = (b11 - b6) / (b11 + b6),
      RE_Slope   = (b10 - b8) / 35,
      VOG1       = b10 / b9,
      GM1        = b10 / b4,
      MCARI      = ((b8 - b7) - 0.2 * (b8 - b4)) * (b8 / b7),
      MCARI3     = ((b10 - b8) - 0.2 * (b10 - b4)) * (b10 / b8),
      CARI       = (b8 - b7) - 0.2 * (b8 - b4),
      TCARI      = 3 * ((b8 - b7) - 0.2 * (b8 - b4) * (b8 / b7)),
      OSAVI      = (1 + 0.16) * (b11 - b7) / (b11 + b7 + 0.16),
      TCA_OSA    = TCARI / OSAVI,
      TVI        = 0.5 * (120 * (b10 - b4) - 200 * (b7 - b4)),
      PSSR_c     = b11 / b2,
      PSNDc      = (b11 - b2) / (b11 + b2),
      CAR        = b8 / b10,
      RARSa      = b7 / b8,
      RARS       = b10 / b3,
      RGI        = b8 / b4,
      RR         = 1 / b8,
      Carter1    = b8 / b1,
      ARI        = (1 / b4) - (1 / b8),
      ARI2       = (1 / b4) - (1 / b9),
      CRI550_515 = (1 / b3) - (1 / b4),
      CRI700_515 = (1 / b3) - (1 / b8),
      DCabCxc    = b7 / (b4 * 3 * b8),
      DNIRCabCxc = b11 / (b4 * b8),
      SIPI       = (b11 - b1) / (b11 + b7),
      PSRI       = (b7 - b2) / b10,
      NPCI       = (b7 - b1) / (b7 + b1),
      PBI        = b11 / b4,
      CVI        = (b11 * b8) / (b4^2),
      CI         = b10 / b8,
      PRI        = (b3 - b4) / (b3 + b4),
      PRIm3      = (b7 - b3) / (b7 + b3),
      PRIm4      = (b4 - b3 - b7) / (b4 + b3 + b7),
      PRIn       = ((b4 - b3) / (b4 + b3)) /
                     (((b11 - b7) / sqrt(b11 + b7)) * (b8 / b7)),
      PRI_CI     = ((b4 - b3) / (b4 + b3)) * ((b10 / b8) - 1),
      CCI        = (b3 - b6) / (b3 + b6),
      LIC3       = b1 / b10,
      GCC        = (b3 + b4) / ((b6 + b7) + (b3 + b4) + (b1 + b2)),
      GNDVI      = (b11 - b4) / (b11 + b4),
      GNDVI2     = (b11 - b3) / (b11 + b3),
      GRVI       = b11 / b3,
      NGRDI      = (b4 - b7) / (b4 + b7),
      VARI       = (b4 - b7) / (b4 + b7 - b2),
      GLI        = (2 * b4 - b7 - b2) / (2 * b4 + b7 + b2),
      DWSI4      = b4 / b7,
      R_index    = b8 / b7,
      B_index    = b1 / b2,
      BGI2       = b1 / b4,
      BRI        = b1 / b8,
      CN1        = (b2 - b7) / (b2 + b7),
      CN2        = b2 / b7,
      SRPI       = b1 / b7,
      LIC4       = b1 / b3,
      LIC5       = b8 / b3,
      NDWI       = (b3 - b11) / (b3 + b11),
      HI         = (b3 - b8) / (b3 + b8) - 0.5 * b8,
      CLS        = (b8 - b4) / (b8 + b4) - b10,
      PMI        = (b3 - b4) / (b3 + b4) + b9,
      SBRI       = (b4 - b3) / (b4 + b3) + 0.5 * b8,
      VIN_B8     = b11 / b8,
      VIN_B9     = b11 / b9,
      VIN_B10    = b11 / b10
    ) %>%
    mutate(across(NDVI:VIN_B10, ~ ifelse(is.infinite(.x) | is.nan(.x), NA, .x)))
}

## Six distributional statistics per index and per raw band, aggregated to crowns.
summarise_crowns <- function(vi_data, feature_cols) {
  vi_data %>%
    group_by(site, plot_id, Label_12, DOTHI___IN) %>%
    summarise(
      n_pixels = n(),
      across(all_of(feature_cols), ~ mean(.x, na.rm = TRUE),           .names = "{.col}_mean"),
      across(all_of(feature_cols), ~ sd(.x, na.rm = TRUE),             .names = "{.col}_sd"),
      across(all_of(feature_cols), ~ median(.x, na.rm = TRUE),         .names = "{.col}_med"),
      across(all_of(feature_cols), ~ quantile(.x, 0.10, na.rm = TRUE), .names = "{.col}_p10"),
      across(all_of(feature_cols), ~ quantile(.x, 0.90, na.rm = TRUE), .names = "{.col}_p90"),
      across(all_of(feature_cols), ~ safe_cv(.x),                      .names = "{.col}_cv"),
      .groups = "drop")
}


# ==============================================================================
# 7. FEATURE SELECTION
# ==============================================================================

STAT_SUFFIXES <- c(mean = "_mean$", sd = "_sd$", median = "_med$",
                   p10 = "_p10$", p90 = "_p90$", cv = "_cv$")

detect_stat_type <- function(features) {
  out <- rep("other", length(features))
  for (s in names(STAT_SUFFIXES)) out[grepl(STAT_SUFFIXES[[s]], features)] <- s
  out
}

strip_stat_suffix <- function(features)
  sub("_(mean|sd|med|p10|p90|cv)$", "", features)

## Rank features by |Spearman rho| and by univariate permutation importance,
## then average the two ranks. Computed on training crowns only.
rank_features <- function(train_data, feature_cols, rf_trees = RF_TREES_UNIV,
                          verbose = TRUE) {
  na_rate <- vapply(feature_cols,
                    function(f) mean(is.na(train_data[[f]])), numeric(1))
  keep <- feature_cols[na_rate <= NA_THRESHOLD]
  if (verbose)
    cat(sprintf("  Features after NA filter (<=%.0f%%): %d of %d\n",
                NA_THRESHOLD * 100, length(keep), length(feature_cols)))

  y <- train_data$DOTHI___IN

  if (verbose) cat("  Spearman correlations...\n")
  rho <- vapply(keep, function(f) {
    x  <- train_data[[f]]
    ok <- !is.na(x) & !is.na(y) & is.finite(x)
    if (sum(ok) < 10) NA_real_ else cor(x[ok], y[ok], method = "spearman")
  }, numeric(1))

  spearman_df <- data.frame(feature = keep, spearman_r = rho, abs_r = abs(rho)) %>%
    filter(!is.na(spearman_r)) %>%
    mutate(stat_type = detect_stat_type(feature),
           vi_name   = strip_stat_suffix(feature)) %>%
    arrange(desc(abs_r)) %>%
    mutate(spearman_rank = row_number())

  if (verbose) cat(sprintf("  Univariate random forests (%d features)...\n", length(keep)))
  t0 <- Sys.time()
  imp <- vapply(seq_along(keep), function(i) {
    if (verbose && i %% 50 == 0)
      cat(sprintf("    %d / %d\n", i, length(keep)))
    d <- data.frame(y = y, x = train_data[[keep[i]]]) %>%
      filter(!is.na(y), !is.na(x), is.finite(x))
    if (nrow(d) < 20) return(NA_real_)
    tryCatch({
      set.seed(SEED)
      ranger::ranger(y ~ x, data = d, num.trees = rf_trees,
                     importance = "permutation",
                     num.threads = 0)$variable.importance[["x"]]
    }, error = function(e) NA_real_)
  }, numeric(1))
  if (verbose)
    cat(sprintf("  Univariate forests done in %.2f min\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  rf_df <- data.frame(feature = keep, rf_imp = imp) %>%
    filter(!is.na(rf_imp)) %>%
    arrange(desc(rf_imp)) %>%
    mutate(rf_rank = row_number())

  spearman_df %>%
    left_join(select(rf_df, feature, rf_imp, rf_rank), by = "feature") %>%
    filter(!is.na(rf_imp)) %>%
    mutate(mean_rank = (spearman_rank + rf_rank) / 2) %>%
    arrange(mean_rank)
}

## Walk the ranking from best to worst, keeping a feature only if it is not
## collinear with anything already kept.
greedy_prune <- function(ordered_features, data, cutoff = COR_CUTOFF) {
  if (!length(ordered_features)) return(character(0))
  cm <- suppressWarnings(cor(data[, ordered_features, drop = FALSE],
                             use = "pairwise.complete.obs", method = "pearson"))
  selected <- character(0)
  for (f in ordered_features) {
    if (!length(selected)) { selected <- f; next }
    mx <- suppressWarnings(max(abs(cm[f, selected]), na.rm = TRUE))
    if (!is.finite(mx) || mx < cutoff) selected <- c(selected, f)
  }
  selected
}

select_features <- function(ranking, train_data, cutoff = COR_CUTOFF,
                            verbose = TRUE) {
  available <- intersect(ranking$feature, names(train_data))
  ordered   <- ranking %>%
    filter(feature %in% available, !is.na(mean_rank)) %>%
    arrange(mean_rank) %>% pull(feature)
  selected <- greedy_prune(ordered, train_data, cutoff)
  if (verbose) {
    cat(sprintf("  Selected %d of %d candidates (|r| < %.2f)\n",
                length(selected), length(ordered), cutoff))
    print(table(detect_stat_type(selected)))
  }
  selected
}

## Repeat the whole selection procedure on plot-level bootstrap resamples of the
## training plots to quantify how reproducible the selected set is.
bootstrap_feature_stability <- function(crown_data, reference_features,
                                        train_plots, meta_cols,
                                        B = BOOTSTRAP_B, rf_trees = BOOTSTRAP_RF) {
  pool <- crown_data %>%
    filter(plot_id %in% train_plots, !is.na(DOTHI___IN)) %>%
    as.data.frame()
  candidates <- setdiff(names(pool), meta_cols)
  candidates <- candidates[vapply(pool[candidates], is.numeric, logical(1))]
  cat(sprintf("  %d candidate features | %d training crowns | B = %d\n",
              length(candidates), nrow(pool), B))

  set.seed(SEED)
  plots <- sort(unique(pool$plot_id))
  sets  <- vector("list", B)
  t0    <- Sys.time()

  for (b in seq_len(B)) {
    boot_plots <- sample(plots, length(plots), replace = TRUE)
    rows <- do.call(rbind, lapply(boot_plots,
                                  function(p) pool[pool$plot_id == p, , drop = FALSE]))
    usable <- candidates[vapply(candidates, function(f) {
      v <- rows[[f]]
      mean(is.na(v)) <= NA_THRESHOLD && is.finite(sd(v, na.rm = TRUE)) &&
        sd(v, na.rm = TRUE) > 0
    }, logical(1))]
    ranking <- rank_features(rows, usable, rf_trees = rf_trees, verbose = FALSE)
    sets[[b]] <- greedy_prune(ranking$feature, rows)
    if (b %% 5 == 0 || b == B)
      cat(sprintf("    bootstrap %d/%d (%.1f min elapsed, mean set size %.1f)\n",
                  b, B, as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  mean(lengths(sets[seq_len(b)]))))
  }

  all_feats <- sort(unique(c(unlist(sets), reference_features)))
  stability <- data.frame(
    feature = all_feats,
    selection_frequency = round(vapply(all_feats, function(f)
      mean(vapply(sets, function(s) f %in% s, logical(1))), numeric(1)), 4),
    in_reference = all_feats %in% reference_features) %>%
    arrange(desc(selection_frequency))
  save_table(stability, "feature_selection_stability.csv")

  jaccard <- function(a, b) {
    u <- length(union(a, b)); if (u == 0) NA_real_ else length(intersect(a, b)) / u
  }
  pairwise <- unlist(lapply(seq_len(B - 1), function(i)
    vapply((i + 1):B, function(j) jaccard(sets[[i]], sets[[j]]), numeric(1))))
  mean_pw  <- mean(pairwise, na.rm = TRUE)
  mean_ref <- mean(vapply(sets, jaccard, numeric(1), reference_features), na.rm = TRUE)

  cat(sprintf("  Reference set: %d features | bootstrap size %.1f +/- %.1f\n",
              length(reference_features), mean(lengths(sets)), sd(lengths(sets))))
  cat(sprintf("  Mean pairwise Jaccard: %.3f | mean Jaccard vs reference: %.3f\n",
              mean_pw, mean_ref))
  cat(sprintf("  Selected in >=80%% of bootstraps: %d\n",
              sum(stability$selection_frequency >= 0.80)))

  plot_df <- stability %>%
    filter(selection_frequency >= 0.25 | in_reference) %>%
    mutate(label = factor(gsub("_", "-", feature),
                          levels = rev(gsub("_", "-", feature))))
  p <- ggplot(plot_df, aes(selection_frequency, label, fill = in_reference)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed",
               colour = "grey55", linewidth = 0.4) +
    scale_fill_manual(values = c(`TRUE` = "#1565c0", `FALSE` = "#bdbdbd"),
                      name = "In reference\nselection",
                      labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
    scale_x_continuous(name = "Bootstrap selection frequency", limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Feature-selection stability (bootstrap of training plots)",
         subtitle = sprintf("B = %d | mean pairwise Jaccard = %.3f | vs reference = %.3f",
                            B, mean_pw, mean_ref)) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(colour = "grey40", size = 8),
          axis.title.y = element_blank(), axis.text = element_text(size = 8),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
  save_fig(p, "feature_selection_stability.png",
           width = 170, height = max(110, 4 * nrow(plot_df)))

  invisible(list(stability = stability, mean_pairwise_jaccard = mean_pw,
                 mean_jaccard_vs_reference = mean_ref, sets = sets))
}


# ==============================================================================
# 8. MODEL FITTING
# ==============================================================================

## Grid search over mtry and min.node.size using small ordinal forests.
tune_ordfor <- function(train_df, tag = "") {
  n_feat <- ncol(train_df) - 1
  mtry_grid <- sort(unique(pmin(OF_MTRY_GRID, n_feat)))
  mtry_grid <- mtry_grid[mtry_grid >= 1]
  cat(sprintf("  [%s] tuning over mtry {%s} x nodesize {%s}\n", tag,
              paste(mtry_grid, collapse = ", "),
              paste(OF_NODESIZE_GRID, collapse = ", ")))

  # mtry outer, min.node.size inner, so ties resolve as in the original loop
  history <- expand.grid(min.node.size = OF_NODESIZE_GRID, mtry = mtry_grid)[, 2:1]
  t0 <- Sys.time()
  history$mean_perf <- vapply(seq_len(nrow(history)), function(i) {
    set.seed(SEED)
    mod <- ordinalForest::ordfor(
      depvar = "health_class", data = train_df,
      nsets = OF_TUNE$nsets, nbest = OF_TUNE$nbest,
      ntreeperdiv = OF_TUNE$ntreeperdiv, ntreefinal = OF_TUNE$ntreefinal,
      perffunction = "equal", importance = "rps",
      mtry = history$mtry[i], min.node.size = history$min.node.size[i],
      num.threads = 0)
    mean(mod$perffunctionvalues)
  }, numeric(1))
  cat(sprintf("  [%s] tuning complete in %.2f min\n", tag,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  history <- history[order(-history$mean_perf), ]
  best <- list(mtry = history$mtry[1], min.node.size = history$min.node.size[1])
  cat(sprintf("  [%s] best: mtry = %d | min.node.size = %d\n",
              tag, best$mtry, best$min.node.size))
  list(best_params = best, history = history)
}

## Out-of-bag performance of the final forest, for comparison with the test set.
eval_ordfor_oob <- function(model, train_df, tag = "") {
  tryCatch({
    borders <- model$bordersbest
    classes <- model$classes
    idx <- vapply(model$forestfinal$predictions,
                  function(x) max(which(x >= borders[seq_along(classes)])), numeric(1))
    cm <- caret::confusionMatrix(
      factor(classes[idx], levels = classes),
      factor(as.character(train_df$health_class), levels = classes))
    cat(sprintf("  [%s] OOB: Acc = %.4f | BalAcc = %.4f | Kappa = %.4f\n", tag,
                cm$overall["Accuracy"],
                mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
                cm$overall["Kappa"]))
    invisible(list(accuracy = unname(cm$overall["Accuracy"]),
                   balanced_acc = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
                   kappa = unname(cm$overall["Kappa"])))
  }, error = function(e) {
    cat(sprintf("  [%s] OOB evaluation failed: %s\n", tag, conditionMessage(e)))
    invisible(NULL)
  })
}

## Fit one ordinal forest for one granularity and return model, predictions and
## metrics.
fit_ordinal_forest <- function(model_data, train_idx, test_idx, features,
                               midpoints, tag, save_outputs = TRUE) {
  cat(sprintf("\n--- Ordinal Forest, %s ---\n", tag))
  train <- model_data[train_idx, ]
  test  <- model_data[test_idx, ]
  of_train <- as.data.frame(train[, c("health_class", features)])
  of_test  <- as.data.frame(test[, features])
  class_levels <- levels(model_data$health_class)

  sparse <- names(table(of_train$health_class))[table(of_train$health_class) < 5]
  if (length(sparse))
    cat(sprintf("  Sparse training classes (n < 5): %s\n", paste(sparse, collapse = ", ")))

  tuning <- tune_ordfor(of_train, tag)
  params <- tuning$best_params
  if (save_outputs)
    save_table(tuning$history, sprintf("tuning_history_%s.csv", tag))

  set.seed(SEED); t0 <- Sys.time()
  model <- ordinalForest::ordfor(
    depvar = "health_class", data = of_train,
    nsets = OF_FINAL$nsets, nbest = OF_FINAL$nbest,
    ntreeperdiv = OF_FINAL$ntreeperdiv, ntreefinal = OF_FINAL$ntreefinal,
    perffunction = "equal", importance = "rps",
    mtry = params$mtry, min.node.size = params$min.node.size, num.threads = 0)
  cat(sprintf("  Fitted in %.2f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  oob  <- eval_ordfor_oob(model, of_train, tag)
  pred <- predict(model, newdata = of_test)
  pred_num <- midpoints[as.character(pred$ypred)]
  true_num <- test$y_numeric
  if (any(is.na(pred_num)))
    warning(sprintf("[%s] some predictions had no midpoint match", tag))

  cm <- caret::confusionMatrix(
    factor(as.character(pred$ypred), levels = class_levels),
    factor(as.character(test$health_class), levels = class_levels))
  metrics <- compute_metrics(cm, pred_num, true_num,
                             max(table(test$health_class)) / nrow(test))
  metrics$qwk <- compute_qwk(test$health_class, pred$ypred, class_levels,
                             sprintf("OF %s", tag))

  cat(sprintf("  Acc = %.4f | BalAcc = %.4f | QWK = %.4f | MAE = %.2f | W5 = %.3f | W10 = %.3f\n",
              metrics$accuracy, metrics$balanced_acc, metrics$qwk,
              metrics$mae, metrics$within_1, metrics$within_2))

  varimp <- data.frame(feature = names(model$varimp),
                       importance = unname(model$varimp)) %>%
    arrange(desc(importance))

  if (save_outputs) {
    saveRDS(model, out_path(sprintf("ordfor_%s_model.rds", tag)))
    save_table(varimp, sprintf("ordfor_%s_varimp.csv", tag))
    save_table(data.frame(
      Label_12        = test$Label_12,
      DOTHI___IN      = true_num,
      True_class      = as.character(test$health_class),
      Predicted_class = as.character(pred$ypred),
      Predicted_mid   = pred_num,
      Residual        = true_num - pred_num,
      Correct         = as.character(test$health_class) == as.character(pred$ypred)),
      sprintf("ordfor_%s_predictions.csv", tag))

    if (HAS_TREESHAP) {
      shap <- tryCatch({
        unified <- treeshap::ranger.unify(model$forestfinal, of_train[, features])
        as.data.frame(treeshap::treeshap(unified, of_test[, features],
                                         interactions = FALSE, verbose = FALSE)$shaps)
      }, error = function(e) {
        cat(sprintf("  SHAP failed for %s: %s\n", tag, conditionMessage(e))); NULL
      })
      if (!is.null(shap))
        save_table(shap, sprintf("ordfor_%s_shap_values.csv", tag))
    }
  }

  list(model = model, tuning = tuning, params = params, oob = oob,
       pred = pred$ypred, pred_numeric = pred_num, true_numeric = true_num,
       test = test, train = train, cm = cm, metrics = metrics,
       varimp = varimp, class_levels = class_levels)
}

## fabOF expects an unordered factor response; y_numeric is retained for the
## tolerance metrics.
build_fabof_data <- function(model_data, features) {
  df <- model_data %>%
    select(Label_12, DOTHI___IN, y_numeric, health_class, all_of(features)) %>%
    filter(!is.na(DOTHI___IN), complete.cases(.))
  df$health_class_fac <- factor(as.character(df$health_class),
                                levels = levels(df$health_class), ordered = FALSE)
  as.data.frame(df)
}

fit_fabof <- function(model_data, train_idx, test_idx, features, midpoints, tag,
                      mtry = FABOF_MTRY, save_outputs = TRUE) {
  cat(sprintf("\n--- fabOF, %s ---\n", tag))
  train <- build_fabof_data(model_data[train_idx, ], features)
  test  <- build_fabof_data(model_data[test_idx, ],  features)
  class_levels <- levels(train$health_class_fac)

  sparse <- names(table(train$health_class_fac))[table(train$health_class_fac) < 5]
  if (length(sparse))
    cat(sprintf("  Sparse training classes (n < 5): %s\n", paste(sparse, collapse = ", ")))

  set.seed(SEED); t0 <- Sys.time()
  model <- fabOF::fabOF(
    formula = as.formula(paste("health_class_fac ~", paste(features, collapse = " + "))),
    data = train, importance = FALSE,
    ranger.control = list(num.trees = FABOF_NTREES,
                          mtry = min(mtry, length(features)), num.threads = 0))
  cat(sprintf("  Fitted in %.2f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  pred <- predict(model, newdata = test)
  pred_num <- midpoints[as.character(pred)]
  true_num <- test$y_numeric

  cm <- caret::confusionMatrix(
    factor(as.character(pred), levels = class_levels),
    factor(as.character(test$health_class_fac), levels = class_levels))
  metrics <- compute_metrics(cm, pred_num, true_num,
                             max(table(test$health_class_fac)) / nrow(test))
  metrics$qwk <- compute_qwk(test$health_class_fac, pred, class_levels,
                             sprintf("fabOF %s", tag))

  cat(sprintf("  Acc = %.4f | BalAcc = %.4f | QWK = %.4f | MAE = %.2f | W5 = %.3f | W10 = %.3f\n",
              metrics$accuracy, metrics$balanced_acc, metrics$qwk,
              metrics$mae, metrics$within_1, metrics$within_2))

  ## Per-crown predictions, mirroring the Ordinal Forest export above. The
  ## spatial validation script re-scores these on progressively more
  ## separated subsets of the test set, so it needs them per crown.
  if (save_outputs && EXPORT_FABOF_PREDICTIONS) {
    save_table(data.frame(
      Label_12        = test$Label_12,
      DOTHI___IN      = true_num,
      True_class      = as.character(test$health_class_fac),
      Predicted_class = as.character(pred),
      Predicted_mid   = pred_num,
      Residual        = true_num - pred_num,
      Correct         = as.character(test$health_class_fac) == as.character(pred)),
      sprintf("fabof_%s_predictions.csv", tag))
  }

  list(model = model, pred = pred, pred_numeric = pred_num,
       true_numeric = true_num, train = train, test = test,
       cm = cm, metrics = metrics, class_levels = class_levels)
}

## Fit OF and fabOF across all four granularities on one feature set.
run_all_granularities <- function(model_base, train_idx, test_idx, features,
                                  save_outputs = TRUE, label = "combined") {
  grans <- make_granularities(model_base$y_numeric)
  of_fits <- list(); fab_fits <- list()
  for (g in GRAN_ORDER) {
    md <- model_base %>% mutate(health_class = grans[[g]]$fn(DOTHI___IN))
    tag <- if (label == "combined") g else paste(label, g, sep = "_")
    of_fits[[g]] <- fit_ordinal_forest(
      md, train_idx, test_idx, features, grans[[g]]$mids, tag,
      save_outputs = save_outputs)
    fab_fits[[g]] <- fit_fabof(
      md, train_idx, test_idx, features, grans[[g]]$mids, tag,
      save_outputs = save_outputs)
  }
  list(of = of_fits, fabof = fab_fits, granularities = grans)
}

## Flatten a run_all_granularities() result into a tidy results table.
summarise_fits <- function(fits, n_features, extra = list()) {
  rows <- lapply(GRAN_ORDER, function(g) {
    lapply(c("OF", "fabOF"), function(algo) {
      m <- if (algo == "OF") fits$of[[g]]$metrics else fits$fabof[[g]]$metrics
      data.frame(granularity = g, algorithm = algo, n_features = n_features,
                 balanced_acc = round(m$balanced_acc, 4),
                 qwk = round(m$qwk, 4), mae = round(m$mae, 3),
                 within_pm5 = round(m$within_1, 4),
                 within_pm10 = round(m$within_2, 4),
                 accuracy = round(m$accuracy, 4), kappa = round(m$kappa, 4))
    })
  })
  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out %>%
    mutate(granularity = factor(granularity, levels = GRAN_ORDER)) %>%
    arrange(granularity, algorithm)
}

print_results_table <- function(df, title) {
  cat(sprintf("\n%s\n", title))
  hdr <- sprintf("  %-10s %-7s %9s %8s %9s %8s %8s",
                 "Gran", "Algo", "BalAcc", "QWK", "MAE(%)", "W5%", "W10%")
  cat(hdr, "\n"); cat("  ", strrep("-", nchar(hdr) - 2), "\n", sep = "")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cat(sprintf("  %-10s %-7s %9.4f %8.4f %9.3f %8.4f %8.4f\n",
                as.character(r$granularity), r$algorithm, r$balanced_acc,
                r$qwk, r$mae, r$within_pm5, r$within_pm10))
  }
  cat("  ", strrep("=", nchar(hdr) - 2), "\n", sep = "")
}


# ==============================================================================
# 9. FIGURES
# ==============================================================================

STAT_COLOURS <- c(p10 = "#1565c0", p90 = "#42a5f5", cv = "#c62828",
                  sd = "#ef9a9a", mean = "#2e7d32", median = "#81c784",
                  other = "grey50")
STAT_LABELS  <- c("10th percentile (p10)", "90th percentile (p90)",
                  "Coefficient of variation (cv)", "Standard deviation (sd)",
                  "Mean", "Median")

## Crown mean reflectance by severity bin, overall and per trial.
plot_spectral_profiles <- function(pixel_data) {
  band_info <- data.frame(
    band       = RAW_BANDS,
    wavelength = c(444, 475, 531, 560, 650, 668, 705, 717, 740, 842),
    band_label = c("B1\n444", "B2\n475", "B3\n531", "B4\n560", "B6\n650",
                   "B7\n668", "B8\n705", "B9\n717", "B10\n740", "B11\n842"))
  breaks <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 95)
  labels <- c("1-10","11-20","21-30","31-40","41-50",
              "51-60","61-70","71-80","81-95")
  palette <- setNames(colorRampPalette(
    c("#1a9641","#a6d96a","#ffffbf","#fdae61","#d7191c"))(length(labels)), labels)

  build <- function(group_vars) {
    pixel_data %>%
      filter(!is.na(DOTHI___IN)) %>%
      mutate(score_bin = cut(DOTHI___IN, breaks = breaks, labels = labels,
                             right = TRUE, include.lowest = TRUE)) %>%
      filter(!is.na(score_bin)) %>%
      group_by(across(all_of(c(group_vars, "score_bin")))) %>%
      summarise(across(all_of(RAW_BANDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      group_by(across(all_of(c(setdiff(group_vars, "Label_12"), "score_bin")))) %>%
      summarise(across(all_of(RAW_BANDS), ~ mean(.x, na.rm = TRUE)),
                n_crowns = n(), .groups = "drop") %>%
      pivot_longer(all_of(RAW_BANDS), names_to = "band", values_to = "reflectance") %>%
      left_join(band_info, by = "band")
  }

  n_crowns <- length(unique(pixel_data$Label_12[!is.na(pixel_data$DOTHI___IN)]))

  base_layers <- function(p) {
    p + scale_colour_manual(values = palette, name = "DNB severity (%)") +
      scale_x_continuous(name = "Wavelength (nm)", breaks = band_info$wavelength,
                         labels = band_info$band_label,
                         expand = expansion(mult = c(0.03, 0.03))) +
      scale_y_continuous(name = "Mean crown reflectance",
                         labels = scales::label_number(accuracy = 0.01),
                         expand = expansion(mult = c(0.02, 0.05))) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11),
            plot.caption = element_text(size = 7, colour = "grey40", hjust = 0,
                                        lineheight = 1.2),
            strip.text = element_text(face = "bold", size = 11),
            axis.text.x = element_text(size = 8, lineheight = 1.1),
            legend.title = element_text(face = "bold", size = 10),
            legend.key.size = unit(0.45, "cm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))
  }

  p_all <- base_layers(
    ggplot(build("Label_12"),
           aes(wavelength, reflectance, colour = score_bin, group = score_bin)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.2, shape = 16) +
      labs(title = "Crown mean spectral reflectance by DNB severity, both trials",
           caption = sprintf(paste("n = %s crowns. Severity bins are",
                                   "10-percentage-point intervals."),
                             format(n_crowns, big.mark = ","))))
  save_fig(p_all, "spectral_profile_combined.png", width = 190, height = 110)

  p_trial <- base_layers(
    ggplot(build(c("site", "Label_12")) %>%
             mutate(site = recode(site, Site1 = "Trial 1", Site2 = "Trial 2")),
           aes(wavelength, reflectance, colour = score_bin, group = score_bin)) +
      geom_line(linewidth = 0.8) + geom_point(size = 1.8, shape = 16) +
      facet_wrap(~ site, ncol = 2) +
      labs(title = "Crown mean spectral reflectance by DNB severity, per trial",
           caption = paste0(
             "Severity bins are 10-percentage-point intervals, coloured from green ",
             "(1-10%) to red (81-95%). Absolute reflectance differences between ",
             "trials reflect residual between-flight radiometric variation after ",
             "panel-based calibration (Aasen et al., 2018). n = ",
             format(n_crowns, big.mark = ","), " crowns.")))
  save_fig(p_trial, "spectral_profile_per_trial.png", width = 190, height = 110)
}

## Map the plot-level split and export the train/test crowns as shapefiles.
export_split_map <- function(crowns_list, split, eligible_labels) {
  all_crowns <- bind_rows(lapply(crowns_list, function(cr)
    cr %>% st_transform(4326) %>% rename_with(clean_names_fn) %>%
      select(Label_12, plot_id, DOTHI___IN, geometry))) %>%
    mutate(split_role = case_when(plot_id %in% split$train_plots ~ "Train",
                                  plot_id %in% split$test_plots  ~ "Test",
                                  TRUE ~ "Unused"),
           site = substr(plot_id, 1, 2))

  for (role in c("Train", "Test")) {
    sub <- all_crowns %>%
      filter(split_role == role, Label_12 %in% eligible_labels) %>% st_make_valid()
    st_write(sub, file.path(DIR_REF, sprintf("%s_crowns.shp", tolower(role))),
             delete_dsn = TRUE, quiet = TRUE)
  }
  cat("  Split shapefiles written to the Reference folder.\n")

  pal <- c(Train = "#1565c0", Test = "#e65100")
  for (s in c("S1", "S2")) {
    site_crowns <- all_crowns %>% filter(site == s) %>% st_make_valid()
    if (!nrow(site_crowns)) next
    n_tr <- length(unique(site_crowns$plot_id[site_crowns$split_role == "Train"]))
    n_te <- length(unique(site_crowns$plot_id[site_crowns$split_role == "Test"]))
    p <- ggplot() +
      geom_sf(data = site_crowns, aes(fill = split_role, colour = split_role),
              size = 0.1, alpha = 0.75) +
      scale_fill_manual(values = pal, name = "Split assignment",
                        labels = c(Train = sprintf("Train (%d plots)", n_tr),
                                   Test  = sprintf("Test (%d plots)", n_te))) +
      scale_colour_manual(values = pal, guide = "none") +
      labs(title = sprintf("Plot-level train/test split, %s",
                           ifelse(s == "S1", "Trial 1", "Trial 2")),
           subtitle = "Whole plots are assigned to one side of the split, so no plot contributes to both.",
           x = "Longitude", y = "Latitude") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(colour = "grey40", size = 9),
            panel.grid = element_line(colour = "grey90", linewidth = 0.3))
    save_fig(p, sprintf("split_map_%s.png", tolower(s)), width = 200, height = 160)
  }
  invisible(all_crowns)
}

## Full class-by-class confusion grid with tolerance banding.
build_cm_grid <- function(true_cat, pred_cat, midpoints, class_levels) {
  observed <- data.frame(True = as.character(true_cat),
                         Predicted = as.character(pred_cat)) %>%
    count(True, Predicted, name = "Freq")
  expand.grid(True = class_levels, Predicted = class_levels,
              stringsAsFactors = FALSE) %>%
    left_join(observed, by = c("True", "Predicted")) %>%
    mutate(Freq = ifelse(is.na(Freq), 0L, as.integer(Freq)),
           abs_diff = abs(midpoints[True] - midpoints[Predicted])) %>%
    group_by(True) %>%
    mutate(row_n = sum(Freq),
           row_prop = ifelse(row_n > 0, Freq / row_n, NA_real_)) %>%
    ungroup()
}

## Side-by-side OF and fabOF confusion matrices for the 19-class model.
## Fill encodes row-wise proportion (where the data sit); the cell border
## encodes prediction tolerance (how accurate the model is); the y-axis
## carries per-grade support so sparse rows are not over-read.
plot_confusion_matrices <- function(of_fit, fab_fit, midpoints, class_levels,
                                    filename = "confusion_matrices_19class.png") {
  band_levels <- c("Exact match (0%)", "Within \u00b15% (adjacent grade)",
                   "Within \u00b110% (two grades)", "Beyond \u00b110%")
  band_cols <- setNames(c("#b45309", "#f59e0b", "#fcd34d", "grey92"), band_levels)

  prep <- function(true_cat, pred_cat) {
    build_cm_grid(true_cat, pred_cat, midpoints, class_levels) %>%
      mutate(True = factor(True, levels = class_levels),
             Predicted = factor(Predicted, levels = class_levels),
             label = ifelse(Freq > 0, as.character(Freq), ""),
             band = factor(case_when(abs_diff == 0  ~ band_levels[1],
                                     abs_diff <= 5  ~ band_levels[2],
                                     abs_diff <= 10 ~ band_levels[3],
                                     TRUE           ~ band_levels[4]),
                           levels = band_levels),
             dark_bg = !is.na(row_prop) & row_prop > 0.45)
  }

  cm_theme <- theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11, margin = margin(b = 3)),
          plot.subtitle = element_text(colour = "grey40", size = 9, margin = margin(b = 6)),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.title = element_text(face = "bold", size = 9),
          legend.text = element_text(size = 8),
          legend.key.width = unit(0.42, "cm"),
          panel.grid = element_blank(), plot.margin = margin(3, 3, 3, 3, "mm"))

  panel <- function(cm_data, title, metrics, show_y) {
    y_labels <- cm_data %>% distinct(True, row_n) %>% arrange(True)
    p <- ggplot(cm_data, aes(Predicted, True)) +
      geom_tile(aes(fill = row_prop, colour = band), linewidth = 0.45) +
      geom_text(data = subset(cm_data, dark_bg), aes(label = label),
                colour = "white", size = 2.1, fontface = "bold") +
      geom_text(data = subset(cm_data, !dark_bg), aes(label = label),
                colour = "grey20", size = 2.1, fontface = "bold") +
      scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "white",
                          limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                          labels = c("0%", "25%", "50%", "75%", "100%"),
                          name = "Share of crowns\nwithin each true\nseverity grade") +
      scale_colour_manual(values = band_cols, drop = FALSE,
                          name = "Prediction tolerance\n(cell border)") +
      guides(fill = guide_colourbar(order = 1, barheight = unit(3.2, "cm")),
             colour = guide_legend(order = 2,
                                   override.aes = list(fill = "white", linewidth = 1.1))) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE,
                       labels = setNames(sprintf("%s  (n=%d)",
                                                 as.character(y_labels$True), y_labels$row_n),
                                         as.character(y_labels$True))) +
      coord_fixed() +
      labs(title = title, x = "Predicted DNB severity grade (%)",
           subtitle = sprintf("QWK = %.3f | W\u00b15%% = %.3f | W\u00b110%% = %.3f",
                              metrics$qwk, metrics$within_1, metrics$within_2)) +
      cm_theme
    if (CM_SHOW_50PCT_REF && "50" %in% class_levels)
      p <- p + geom_hline(yintercept = match("50", class_levels) + 0.5,
                          colour = "grey35", linetype = "dotted", linewidth = 0.4)
    if (show_y) p + labs(y = "True DNB severity grade (%)")
    else p + labs(y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  cm_of  <- prep(as.character(of_fit$test$health_class), as.character(of_fit$pred))
  cm_fab <- prep(as.character(fab_fit$test$health_class_fac), as.character(fab_fit$pred))

  n_test <- nrow(of_fit$test)
  pct_le50 <- 100 * sum(as.numeric(as.character(of_fit$test$health_class)) <= 50,
                        na.rm = TRUE) / n_test

  combined <-
    (panel(cm_of, "(a) Ordinal Forest (OF)", of_fit$metrics, TRUE) +
       theme(legend.position = "none")) +
    panel(cm_fab, "(b) Frequency-adjusted borders OF (fabOF)", fab_fit$metrics, FALSE) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(caption = paste0(
      "Crown-level Dothistroma needle blight severity, 19-class model. Severity is graded ",
      "in 5-percentage-point steps, so a prediction within \u00b15% falls in the grade ",
      "immediately adjacent to the true grade and one within \u00b110% falls two grades ",
      "away; cell borders encode these tolerance bands.\n",
      "Cell shading gives the proportion of crowns within each true severity grade (row) ",
      "assigned to a given predicted grade, and cell numbers are crown counts. Support per ",
      "true grade appears in the y-axis labels, so proportions based on few crowns can be ",
      "read alongside their counts. ",
      sprintf("Reference data are concentrated at low severity: %.0f%% of test crowns fall ",
              pct_le50),
      "at \u226450% severity, so performance above this threshold rests on sparse ",
      "observations.\n",
      "W\u00b15% and W\u00b110% are the overall proportions of predictions falling within ",
      "those tolerances. Test set: n = ", n_test, " crowns."),
      theme = theme(plot.caption = element_text(size = 7.5, colour = "grey40",
                                                hjust = 0, lineheight = 1.3))) &
    theme(legend.position = "right")

  save_fig(combined, filename, width = 285, height = 175)
}

## Top-20 mean |SHAP|, coloured by the crown-level statistic each feature uses.
plot_shap_importance <- function(shap_csv, filename, title, subtitle) {
  if (!file.exists(out_path(shap_csv))) {
    cat(sprintf("  %s not found - skipping SHAP figure.\n", shap_csv))
    return(invisible(NULL))
  }
  shap <- read.csv(out_path(shap_csv), check.names = FALSE)
  imp <- data.frame(feature = names(shap),
                    mean_shap = colMeans(abs(shap), na.rm = TRUE)) %>%
    arrange(desc(mean_shap)) %>% slice_head(n = 20) %>%
    mutate(stat_type = detect_stat_type(feature),
           label = factor(gsub("_", "-", feature),
                          levels = rev(gsub("_", "-", feature))))
  p <- ggplot(imp, aes(mean_shap, label, fill = stat_type)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = STAT_COLOURS, name = "Statistic type",
                      breaks = names(STAT_COLOURS)[1:6], labels = STAT_LABELS) +
    scale_x_continuous(name = "Mean |SHAP value|",
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, subtitle = subtitle) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(colour = "grey40", size = 8, lineheight = 1.2),
          axis.title.y = element_blank(), axis.text = element_text(size = 9),
          legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
  save_fig(p, filename, width = 160, height = 130)
}

plot_class_distribution <- function(model_data, train_idx, test_idx) {
  dist <- bind_rows(
    model_data[train_idx, ] %>% mutate(split = "Training"),
    model_data[test_idx,  ] %>% mutate(split = "Test")) %>%
    mutate(severity = as.numeric(as.character(health_class))) %>%
    count(severity, split, name = "n") %>%
    mutate(split = factor(split, levels = c("Training", "Test")),
           severity = factor(severity, levels = sort(unique(severity))))
  p <- ggplot(dist, aes(severity, n, fill = split)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.85,
             colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c(Training = "#1565c0", Test = "#e65100"),
                      name = "Dataset") +
    scale_x_discrete(name = "Severity grade (%)") +
    scale_y_continuous(name = "Number of crowns", expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Crown distribution across severity classes",
         subtitle = sprintf("19-class model (n = %d training, n = %d test crowns)",
                            length(train_idx), length(test_idx))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          axis.text = element_text(size = 8),
          legend.title = element_text(face = "bold", size = 9),
          legend.key.size = unit(0.4, "cm"),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  save_fig(p, "class_distribution_train_test.png", width = 190, height = 100)
}


# ==============================================================================
# 10. MAIN PIPELINE
# ==============================================================================

pipeline_start <- Sys.time()
cat("\n======================================================================\n")
cat(" RPBC DOTHISTROMA NEEDLE BLIGHT PIPELINE\n")
cat("======================================================================\n")


## ---- 10.1 Sunlit crown pixels ----------------------------------------------
## Two entry points. "rasters" runs the full preprocessing chain from the
## orthomosaics and crown polygons, caching the extracted pixels on the way.
## "pixels" reads that cache and starts at index calculation, which avoids
## repeating the masking and extraction step on every run.

crowns_s1 <- crowns_s2 <- NULL   # only populated on the raster path

if (START_FROM == "rasters") {

  cat("\n[1] Loading rasters and crown reference data\n")

  raster_s1 <- rast(RASTER_SITE1)
  raster_s2 <- rast(RASTER_SITE2)
  crs(raster_s1) <- "EPSG:4326"   # relabels the CRS, does not reproject (Trial 1 only)

  crowns_s1 <- st_read(CROWNS_SITE1, quiet = TRUE) %>%
    st_transform(4326) %>%
    mutate(plot_id = paste0("S1_", Pset_12)) %>%
    fix_health()

  crowns_s2 <- st_read(CROWNS_SITE2, quiet = TRUE) %>%
    st_make_valid() %>% st_transform(4326) %>%
    rename(DOTHI___IN = Dothi___In, Label_12 = FieldLabel) %>%
    mutate(plot_id = paste0("S2_", Pset)) %>%
    fix_health()

  cat(sprintf("  Trial 1: %d crowns across %d plots (%d dead)\n",
              nrow(crowns_s1), length(unique(crowns_s1$plot_id)),
              sum(crowns_s1$is_dead, na.rm = TRUE)))
  cat(sprintf("  Trial 2: %d crowns across %d plots (%d dead)\n",
              nrow(crowns_s2), length(unique(crowns_s2$plot_id)),
              sum(crowns_s2$is_dead, na.rm = TRUE)))

  ## Reflectance conversion and shadow masking
  cat("\n[2] Converting to reflectance and masking shadow\n")

  raster_s1_ref <- raster_s1 / DN_SCALE
  raster_s2_ref <- raster_s2 / DN_SCALE

  for (r in list(list(raster_s1_ref, "Trial 1"), list(raster_s2_ref, "Trial 2"))) {
    rng <- range(values(r[[1]][[NIR_BAND]]), na.rm = TRUE)
    if (rng[2] > 1.5)
      warning(sprintf("%s NIR reflectance max = %.3f after /%d; input may not be 16-bit DN.",
                      r[[2]], rng[2], DN_SCALE))
  }

  if (RUN_QC_PLOTS) {
    overlap_check <- function(raster_obj, crowns, label, out_file) {
      bb <- st_bbox(crowns)
      pad_x <- (bb[3] - bb[1]) * 0.05; pad_y <- (bb[4] - bb[2]) * 0.05
      r <- crop(raster_obj[[7]], ext(bb[1] - pad_x, bb[3] + pad_x,
                                     bb[2] - pad_y, bb[4] + pad_y))
      v <- values(r, na.rm = TRUE)
      r <- clamp(r, lower = quantile(v, 0.05, na.rm = TRUE),
                 upper = quantile(v, 0.95, na.rm = TRUE))
      png(out_file, width = 3000, height = 2600, res = 300)
      par(mar = c(4, 4, 3, 6))
      plot(r, col = grey.colors(512), axes = TRUE, legend = TRUE)
      plot(vect(crowns), add = TRUE, border = "red",
           col = adjustcolor("red", alpha.f = 0.15), lwd = 0.8)
      title(main = sprintf("Crown overlap check, %s", label))
      dev.off()
    }
    overlap_check(raster_s1, crowns_s1, "Trial 1",
                  file.path(DIR_REF, "overlap_check_trial1.png"))
    overlap_check(raster_s2, crowns_s2, "Trial 2",
                  file.path(DIR_REF, "overlap_check_trial2.png"))
    cat("  Overlap check plots written.\n")
  }

  nir_s1 <- find_nir_threshold(raster_s1_ref, vect(crowns_s1), "Trial 1")
  nir_s2 <- find_nir_threshold(raster_s2_ref, vect(crowns_s2), "Trial 2")
  if (RUN_QC_PLOTS) {
    plot_nir_distribution(nir_s1, "Trial1")
    plot_nir_distribution(nir_s2, "Trial2")
  }

  masks_s1 <- build_shadow_mask(raster_s1_ref, vect(crowns_s1), nir_s1$threshold)
  masks_s2 <- build_shadow_mask(raster_s2_ref, vect(crowns_s2), nir_s2$threshold)

  raster_s1_masked <- mask(crop(raster_s1_ref, masks_s1$clean), masks_s1$clean, inverse = TRUE)
  raster_s2_masked <- mask(crop(raster_s2_ref, masks_s2$clean), masks_s2$clean, inverse = TRUE)

  gdal_opts <- c("COMPRESS=DEFLATE", "PREDICTOR=2")
  writeRaster(masks_s1$clean, file.path(DIR_SHADOW, "Trial1_shadow_mask.tif"),
              overwrite = TRUE, gdal = gdal_opts, datatype = "INT1U")
  writeRaster(masks_s2$clean, file.path(DIR_SHADOW, "Trial2_shadow_mask.tif"),
              overwrite = TRUE, gdal = gdal_opts, datatype = "INT1U")
  writeRaster(raster_s1_masked, out_path("Trial1_reflectance_masked.tif"),
              overwrite = TRUE, gdal = gdal_opts, datatype = "FLT4S")
  writeRaster(raster_s2_masked, out_path("Trial2_reflectance_masked.tif"),
              overwrite = TRUE, gdal = gdal_opts, datatype = "FLT4S")
  cat("  Shadow masks and masked reflectance rasters written.\n")

  ## Extract sunlit pixels within crowns
  cat("\n[3] Extracting pixels within crown boundaries\n")

  pixels <- bind_rows(
    extract_pixels(raster_s1_masked, crowns_s1, "Site1"),
    extract_pixels(raster_s2_masked, crowns_s2, "Site2")) %>%
    filter(!is.na(b1), !is_dead, !is.na(DOTHI___IN)) %>%
    select(all_of(PIXEL_TABLE_COLS))

  ## Cached gzipped so subsequent runs can use START_FROM = "pixels".
  gz <- gzfile(PIXEL_TABLE, "w")
  write.csv(pixels, gz, row.names = FALSE)
  close(gz)
  size_mb <- file.size(PIXEL_TABLE) / 1024^2
  cat(sprintf("  Pixel table written: %s (%s rows, %.1f MB)\n",
              PIXEL_TABLE, format(nrow(pixels), big.mark = ","), size_mb))


} else if (START_FROM == "pixels") {

  cat("\n[1-3] Reading the cached pixel table\n")
  if (!file.exists(PIXEL_TABLE))
    stop("Cached pixel table not found: ", PIXEL_TABLE,
         "\nRun once with START_FROM = 'rasters' to generate it from the imagery.")

  pixels <- read.csv(PIXEL_TABLE, stringsAsFactors = FALSE)
  missing_cols <- setdiff(PIXEL_TABLE_COLS, names(pixels))
  if (length(missing_cols))
    stop("Pixel table is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  pixels <- pixels[, PIXEL_TABLE_COLS]
  cat(sprintf("  Loaded %s pixels from %d crowns across %d plots\n",
              format(nrow(pixels), big.mark = ","),
              length(unique(pixels$Label_12)), length(unique(pixels$plot_id))))

} else {
  stop("START_FROM must be either 'rasters' or 'pixels'.")
}

cat(sprintf("  Valid sunlit pixels: %s\n", format(nrow(pixels), big.mark = ",")))


## ---- 10.2 Vegetation indices and crown-level features ----------------------
cat("\n[4] Calculating vegetation indices and crown-level features\n")

vi_data <- calc_vegetation_indices(pixels)
vi_cols <- names(vi_data)[which(names(vi_data) == "NDVI"):which(names(vi_data) == "VIN_B10")]
cat(sprintf("  Vegetation indices calculated: %d\n", length(vi_cols)))

feature_inputs <- c(vi_cols, RAW_BANDS)
crown_features <- summarise_crowns(vi_data, feature_inputs)
names(crown_features) <- clean_names_fn(names(crown_features))
cat(sprintf("  Crown table: %d crowns x %d features (%d inputs x 6 statistics)\n",
            nrow(crown_features), ncol(crown_features) - 5,
            length(feature_inputs)))

if (EXPORT_GENETICS)
  save_table(crown_features, "crown_features_for_genetic_analysis.csv")

plot_spectral_profiles(vi_data)


## ---- 10.3 Plot-level train/test split (before feature selection) -----------
cat("\n[5] Defining the plot-level train/test split\n")

valid_crowns <- crown_features %>% filter(!is.na(DOTHI___IN))
all_plots <- sort(unique(valid_crowns$plot_id))

set.seed(SEED)
train_plots <- sample(all_plots, round(length(all_plots) * TRAIN_FRACTION))
test_plots  <- setdiff(all_plots, train_plots)

cat(sprintf("  Plots: %d total, %d train, %d test\n",
            length(all_plots), length(train_plots), length(test_plots)))
cat(sprintf("  Crowns: %d train, %d test\n",
            sum(valid_crowns$plot_id %in% train_plots),
            sum(valid_crowns$plot_id %in% test_plots)))


## ---- 10.4 Feature ranking and selection (training crowns only) -------------
cat("\n[6] Ranking and selecting features on training crowns only\n")

META_COLS <- c("site", "plot_id", "Label_12", "DOTHI___IN", "n_pixels")
candidate_features <- setdiff(names(crown_features), META_COLS)
train_crowns <- crown_features %>%
  filter(!is.na(DOTHI___IN), plot_id %in% train_plots)

full_ranking <- rank_features(train_crowns, candidate_features)
save_table(full_ranking, "feature_ranking_full.csv")

selected_features <- select_features(
  full_ranking, train_crowns[, setdiff(names(train_crowns), META_COLS)])
writeLines(selected_features, out_path("selected_features.txt"))
save_table(full_ranking %>% filter(feature %in% selected_features) %>%
             arrange(mean_rank) %>% mutate(selection_rank = row_number()),
           "selected_features_ranking.csv")

cat(sprintf("  Vegetation index features: %d | raw band features: %d\n",
            sum(!grepl("^b[0-9]", selected_features)),
            sum(grepl("^b[0-9]", selected_features))))


## ---- 10.5 Model dataset and split indices ----------------------------------
model_base <- crown_features %>%
  select(Label_12, plot_id, DOTHI___IN, all_of(selected_features)) %>%
  filter(!is.na(DOTHI___IN), complete.cases(.)) %>%
  mutate(y_numeric = as.numeric(DOTHI___IN))

train_idx <- which(model_base$plot_id %in% train_plots)
test_idx  <- which(model_base$plot_id %in% test_plots)
cat(sprintf("\n  Model dataset: %d crowns (%d train, %d test)\n",
            nrow(model_base), length(train_idx), length(test_idx)))

## The split map needs the crown polygons, so it is only produced on the raster
## path. The split itself is fully determined by SEED and is identical either way.
if (!is.null(crowns_s1) && !is.null(crowns_s2)) {
  export_split_map(list(crowns_s1, crowns_s2),
                   list(train_plots = train_plots, test_plots = test_plots),
                   model_base$Label_12)
} else {
  cat("  Crown geometries unavailable - skipping the split map.\n")
}

granularities <- make_granularities(model_base$y_numeric)
severity_levels_19 <- as.character(sort(unique(model_base$y_numeric)))
cat("\n  Class distributions:\n")
for (g in GRAN_ORDER) {
  cat(sprintf("  %s:\n", g))
  print(table(granularities[[g]]$fn(model_base$DOTHI___IN)))
}


## ---- 10.6 Bootstrap feature-selection stability ----------------------------
if (RUN_STABILITY) {
  cat("\n[7] Bootstrap feature-selection stability\n")
  stability <- bootstrap_feature_stability(
    crown_features, selected_features, train_plots, META_COLS)
}


## ---- 10.7 Combined-trial models --------------------------------------------
cat("\n[8] Fitting combined-trial models across all granularities\n")

fits <- run_all_granularities(model_base, train_idx, test_idx, selected_features)

results_combined <- summarise_fits(fits, length(selected_features))
save_table(results_combined, "results_combined_trials.csv")
print_results_table(results_combined, "COMBINED-TRIAL MODEL PERFORMANCE")


## ---- 10.8 Publication figures ----------------------------------------------
cat("\n[9] Generating publication figures\n")

plot_class_distribution(
  model_base %>% mutate(health_class = granularities[["19-class"]]$fn(DOTHI___IN)),
  train_idx, test_idx)

plot_confusion_matrices(
  fits$of[["19-class"]], fits$fabof[["19-class"]],
  granularities[["19-class"]]$mids, severity_levels_19)

plot_shap_importance(
  "ordfor_19-class_shap_values.csv",
  "shap_importance_OF_19class.png",
  "SHAP feature importance, ordinal forest (19-class)",
  sprintf(paste("Mean absolute SHAP value across test crowns (n = %d);",
                "colour denotes crown-level statistic type"),
          length(test_idx)))

## fabOF SHAP is supplementary: values describe the latent continuous score the
## model thresholds internally, not the final ordinal class.
if (HAS_TREESHAP) {
  fab_fg <- fits$fabof[["19-class"]]
  n_shap <- min(300, nrow(fab_fg$test))
  shap_fab <- tryCatch({
    unified <- treeshap::ranger.unify(fab_fg$model$ranger.fit,
                                      fab_fg$train[, selected_features])
    as.data.frame(treeshap::treeshap(unified,
                                     fab_fg$test[seq_len(n_shap), selected_features],
                                     interactions = FALSE, verbose = FALSE)$shaps)
  }, error = function(e) {
    cat(sprintf("  fabOF SHAP failed: %s\n", conditionMessage(e))); NULL
  })
  if (!is.null(shap_fab)) {
    save_table(shap_fab, "fabof_19-class_shap_values.csv")
    plot_shap_importance(
      "fabof_19-class_shap_values.csv",
      "shap_importance_fabOF_19class.png",
      "SHAP feature importance, fabOF (19-class, latent score)",
      sprintf(paste("Values reflect contributions to the latent regression score,",
                    "not the final ordinal class. n = %d test crowns."), n_shap))
  }
}


## ---- 10.9 Site-level performance -------------------------------------------
cat("\n[10] Site-level performance\n")

## MAE and the tolerance metrics use the continuous field score against the
## predicted class midpoint, identically to the overall table.
site_metrics <- function(fit, algorithm, granularity) {
  truth <- if (algorithm == "OF") fit$test$health_class else fit$test$health_class_fac
  df <- data.frame(Label_12 = as.character(fit$test$Label_12),
                   true_class = as.character(truth),
                   pred_class = as.character(fit$pred),
                   pred_mid = fit$pred_numeric,
                   true_cont = fit$true_numeric) %>%
    left_join(crown_features %>%
                transmute(Label_12 = as.character(Label_12), site), by = "Label_12") %>%
    filter(!is.na(site)) %>%
    mutate(abs_diff = abs(pred_mid - true_cont))

  do.call(rbind, lapply(sort(unique(df$site)), function(s) {
    d <- df[df$site == s, ]
    cm <- caret::confusionMatrix(
      factor(d$pred_class, levels = fit$class_levels),
      factor(d$true_class, levels = fit$class_levels))
    data.frame(
      site = s, algorithm = algorithm, granularity = granularity,
      n_crowns = nrow(d),
      balanced_acc = round(mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE), 4),
      qwk = round(compute_qwk(d$true_class, d$pred_class, fit$class_levels), 4),
      mae = round(mean(d$abs_diff, na.rm = TRUE), 3),
      within_pm5 = round(mean(d$abs_diff <= 5, na.rm = TRUE), 4),
      within_pm10 = round(mean(d$abs_diff <= 10, na.rm = TRUE), 4))
  }))
}

site_results <- do.call(rbind, unlist(lapply(GRAN_ORDER, function(g) list(
  site_metrics(fits$of[[g]], "OF", g),
  site_metrics(fits$fabof[[g]], "fabOF", g))), recursive = FALSE)) %>%
  mutate(granularity = factor(granularity, levels = GRAN_ORDER)) %>%
  arrange(granularity, algorithm, site)

overall_rows <- results_combined %>%
  transmute(site = "Overall", algorithm, granularity,
            n_crowns = length(test_idx), balanced_acc, qwk, mae,
            within_pm5, within_pm10)

site_results_full <- bind_rows(overall_rows, site_results) %>%
  mutate(site = factor(site, levels = c("Overall", sort(unique(site_results$site))))) %>%
  arrange(granularity, algorithm, site)
save_table(site_results_full, "results_site_level.csv")

cat("\n  Note: at 3- and 5-class granularity the tolerance metrics approach exact\n")
cat("  class agreement because the midpoints are widely spaced; read them\n")
cat("  alongside QWK rather than on their own.\n")


## ---- 10.10 Per-trial models ------------------------------------------------
## Feature selection is repeated independently within each trial, using the same
## training plots, so the transferability of the selected set can be assessed.
if (RUN_TRIAL_MODELS) {
  cat("\n[11] Per-trial pipelines\n")

  run_trial <- function(site_code, trial_label) {
    cat(sprintf("\n=== %s ===\n", trial_label))
    dir_trial <- out_path(tolower(gsub(" ", "", trial_label)))
    dir.create(dir_trial, showWarnings = FALSE, recursive = TRUE)

    trial_crowns <- crown_features %>%
      filter(site == site_code, !is.na(DOTHI___IN))
    trial_train <- trial_crowns %>% filter(plot_id %in% train_plots)
    cat(sprintf("  Training crowns for selection: %d of %d\n",
                nrow(trial_train), nrow(trial_crowns)))

    ranking <- rank_features(trial_train, candidate_features)
    write.csv(ranking, file.path(dir_trial, "feature_ranking_full.csv"), row.names = FALSE)
    features <- select_features(
      ranking, trial_train[, setdiff(names(trial_train), META_COLS)])
    shared <- intersect(features, selected_features)
    cat(sprintf("  Shared with combined-model selection: %d of %d\n",
                length(shared), length(features)))
    write.csv(data.frame(feature = features,
                         in_combined = features %in% selected_features,
                         is_raw_band = grepl("^b[0-9]", features)),
              file.path(dir_trial, "selected_features.csv"), row.names = FALSE)

    base <- trial_crowns %>%
      select(Label_12, plot_id, DOTHI___IN, all_of(features)) %>%
      filter(complete.cases(.)) %>%
      mutate(y_numeric = as.numeric(DOTHI___IN))

    trial_fits <- run_all_granularities(
      base,
      which(base$plot_id %in% train_plots),
      which(base$plot_id %in% test_plots),
      features, save_outputs = FALSE,
      label = tolower(gsub(" ", "", trial_label)))

    results <- summarise_fits(trial_fits, length(features),
                              extra = list(pipeline = trial_label,
                                           shared_with_combined = length(shared)))
    write.csv(results, file.path(dir_trial, "results.csv"), row.names = FALSE)
    print_results_table(results, sprintf("%s MODEL PERFORMANCE", toupper(trial_label)))

    list(label = trial_label, features = features, results = results)
  }

  trial1 <- run_trial("Site1", "Trial 1")
  trial2 <- run_trial("Site2", "Trial 2")

  ## Feature-set overlap across the three pipelines
  shared_all <- Reduce(intersect, list(selected_features, trial1$features, trial2$features))
  cat("\n  FEATURE-SET OVERLAP\n")
  cat(sprintf("    Combined: %d | Trial 1: %d | Trial 2: %d\n",
              length(selected_features), length(trial1$features), length(trial2$features)))
  cat(sprintf("    Shared across all three: %d | Trial 1 and 2: %d\n",
              length(shared_all), length(intersect(trial1$features, trial2$features))))

  feature_overlap <- data.frame(
    feature = sort(unique(c(selected_features, trial1$features, trial2$features)))) %>%
    mutate(in_combined = feature %in% selected_features,
           in_trial1 = feature %in% trial1$features,
           in_trial2 = feature %in% trial2$features,
           shared_all_three = in_combined & in_trial1 & in_trial2,
           stat_type = detect_stat_type(feature),
           combined_rank = match(feature, selected_features),
           trial1_rank = match(feature, trial1$features),
           trial2_rank = match(feature, trial2$features)) %>%
    arrange(desc(shared_all_three), combined_rank, trial1_rank, trial2_rank)
  save_table(feature_overlap, "feature_overlap_three_pipelines.csv")

  wb <- createWorkbook()
  addWorksheet(wb, "Feature overlap")
  writeData(wb, "Feature overlap", feature_overlap)
  addStyle(wb, "Feature overlap",
           createStyle(fontColour = "#FFFFFF", fgFill = "#1565c0",
                       halign = "CENTER", textDecoration = "bold"),
           rows = 1, cols = seq_len(ncol(feature_overlap)), gridExpand = TRUE)
  shared_rows <- which(feature_overlap$shared_all_three) + 1
  if (length(shared_rows))
    addStyle(wb, "Feature overlap",
             createStyle(textDecoration = "bold", fgFill = "#e3f2fd"),
             rows = shared_rows, cols = seq_len(ncol(feature_overlap)),
             gridExpand = TRUE)
  setColWidths(wb, "Feature overlap", cols = seq_len(ncol(feature_overlap)), widths = "auto")
  saveWorkbook(wb, out_path("feature_overlap_three_pipelines.xlsx"), overwrite = TRUE)
  cat("  Table saved: feature_overlap_three_pipelines.xlsx\n")

  ## Three-pipeline performance comparison
  comparison <- bind_rows(
    results_combined %>% mutate(pipeline = "Combined", shared_with_combined = NA_integer_),
    trial1$results, trial2$results) %>%
    arrange(granularity, algorithm, pipeline)
  save_table(comparison, "results_three_pipeline_comparison.csv")
}


## ---- 10.11 Session information ---------------------------------------------
writeLines(capture.output(sessionInfo()), out_path("sessionInfo.txt"))

cat("\n======================================================================\n")
cat(sprintf(" PIPELINE COMPLETE in %.1f minutes\n",
            as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))))
cat(sprintf(" Outputs written to: %s\n", DIR_ANALYSIS))
cat("======================================================================\n\n")