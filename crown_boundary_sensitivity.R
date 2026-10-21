# ==============================================================================
# Crown boundary sensitivity of the DNB severity models
#
# Companion analysis: the models depend on within-crown distributional
#                     statistics, which depend in turn on which pixels fall
#                     inside the crown boundary. This perturbs the boundaries
#                     and asks whether the features, and the predictions built
#                     on them, survive.
#
# Method    Every crown polygon is eroded and dilated by one and two pixels at
#           each trial's native ground sampling distance. Pixel extraction,
#           index calculation and the within-crown statistics are recomputed
#           from each perturbed boundary set. The model is fitted once on the
#           original training features and is not refitted, so any change in
#           prediction is attributable to the boundaries alone.
# Outputs   Appendix Table 6 and the figures quoted in Sections 2.6 and 4
#
# Reproducibility notes
#   - Edit the CONFIGURATION block below; nothing else needs changing.
#   - Reads the masked reflectance rasters, the crown feature table, the
#     selected-feature list and the saved test predictions written by the main
#     analysis pipeline. Run that first.
#   - Only the indices required by the selected features are evaluated, which
#     is why the full index list is held as expressions rather than computed.
#   - All random operations are seeded via SEED.
#   - Runtime is dominated by pixel extraction, roughly ten minutes per
#     boundary variant on a local disk and considerably more over a network
#     share. Copy the rasters locally first if they are not already.
# ==============================================================================


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

## ---- Paths -----------------------------------------------------------------
DIR_BASE     <- "J:/RPBC_Dothi_Analysis"          # <- edit this
DIR_REF      <- file.path(DIR_BASE, "Reference")
DIR_ANALYSIS <- file.path(DIR_BASE, "Analysis")

## Crown polygons and the shadow-masked reflectance rasters. Not distributed
## with this repository; see the data availability statement in the README.
CROWNS_SITE1 <- file.path(DIR_REF, "Crown_Health_Reference_2026.shp")
CROWNS_SITE2 <- file.path(DIR_REF, "Crown_Health_Reference_2026_site2_clean.shp")

## Artefacts written by the main pipeline, looked for in DIR_ANALYSIS. Output
## names have changed between pipeline versions, so each entry lists the names
## in use and the first that exists is taken. Edit a list if yours differs.
NAMES_MASKED_SITE1   <- c("Trial1_reflectance_masked.tif",
                          "Site1_Reflectance_Masked_11Bands.tif")
NAMES_MASKED_SITE2   <- c("Trial2_reflectance_masked.tif",
                          "Site2_Reflectance_Masked_11Bands.tif")
NAMES_FEATURE_LIST   <- c("selected_features.txt",
                          "combined_selected_features_list.txt")
NAMES_CROWN_FEATURES <- c("crown_features_for_genetic_analysis.csv",
                          "crown_summaries_ext_for_genetics.csv")
NAMES_TEST_PRED_OF   <- c("ordfor_19-class_predictions.csv",
                          "ordfor_finegrained_predictions.csv")

## ---- Run control -----------------------------------------------------------
SEED        <- 42
RUN_PREDICT <- TRUE        # refit fabOF once and score each perturbed test set

## ---- Perturbation ----------------------------------------------------------
## Expressed in pixels so that the perturbation is one sensor element at each
## trial, which is the natural unit of delineation precision. Because the two
## trials differ in ground sampling distance, the same buffer is a larger
## fraction of a Trial 2 crown; the realised change in pixel count is reported.
GSD_M          <- c(Site1 = 0.0275, Site2 = 0.0430)
BUFFERS_PX     <- c(0, -1, 1, -2, 2)
BUFFER_LABELS  <- c("original", "eroded 1 px", "dilated 1 px",
                    "eroded 2 px", "dilated 2 px")

## ---- Model -----------------------------------------------------------------
## As used for the reported 19-class fabOF model. Not re-tuned here.
FABOF_TREES    <- 5000
FABOF_MTRY     <- 10
SEVERITY_SCALE <- seq(5, 95, by = 5)

## ---- Within-crown statistics -----------------------------------------------
CV_MAX         <- 10       # treat a coefficient of variation above this as missing
QUANTILES      <- c(p10 = 0.10, p90 = 0.90)


# ==============================================================================
# 2. LIBRARIES AND HELPERS
# ==============================================================================
library(sf); library(terra); library(dplyr); library(tidyr)
library(caret); library(rlang)
if (RUN_PREDICT) library(fabOF)

clean_names <- function(x) gsub("[^a-zA-Z0-9_]", "_", x)
report      <- function(pass, msg) cat(sprintf("  [%s] %s\n",
                                               ifelse(pass, "PASS", "FAIL"), msg))

## Take the first of the candidate filenames that exists, and say plainly what
## was missing if none do.
resolve_input <- function(candidates, what) {
  found <- file.path(DIR_ANALYSIS, candidates)
  found <- found[file.exists(found)]
  if (!length(found))
    stop(sprintf("Cannot find the %s in %s.\n  Looked for: %s\n  Edit the candidate list in the CONFIGURATION block.",
                 what, DIR_ANALYSIS, paste(candidates, collapse = ", ")), call. = FALSE)
  found[1]
}

qwk <- function(obs, pred, levels_all = SEVERITY_SCALE) {
  o <- factor(as.character(obs),  levels = as.character(levels_all))
  p <- factor(as.character(pred), levels = as.character(levels_all))
  observed <- table(o, p); n <- length(levels_all)
  weights  <- outer(seq_len(n), seq_len(n), function(i, j) (i - j)^2)
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  1 - sum(weights * observed) / sum(weights * expected)
}

score <- function(obs, pred, label) {
  lev <- as.character(SEVERITY_SCALE)
  cm  <- caret::confusionMatrix(factor(as.character(pred), levels = lev),
                                factor(as.character(obs),  levels = lev))
  err <- pred - obs
  data.frame(boundary = label, n_crowns = length(obs),
             bal_acc = round(mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE) * 100, 2),
             qwk     = round(qwk(obs, pred), 3),
             mae     = round(mean(abs(err)), 2),
             bias    = round(mean(err), 2),
             w5      = round(mean(abs(err) <=  5) * 100, 2),
             w10     = round(mean(abs(err) <= 10) * 100, 2))
}

## A coefficient of variation is meaningless for an index whose mean sits near
## zero, so implausible values are returned as missing rather than propagated.
safe_cv <- function(x) {
  v <- x[!is.na(x) & is.finite(x)]
  if (length(v) < 2) return(NA_real_)
  m <- mean(v); s <- sd(v)
  if (is.na(m) || is.na(s) || m == 0) return(NA_real_)
  ratio <- s / abs(m)
  if (ratio > CV_MAX) NA_real_ else ratio
}


# ==============================================================================
# 3. VEGETATION INDICES
#
# The same 82 indices as the main pipeline, held as expressions so that only
# those required by the selected features are evaluated. Band centres for the
# MicaSense RedEdge-P Dual: b1 444, b2 475, b3 531, b4 560, b5 panchromatic
# (unused), b6 650, b7 668, b8 705, b9 717, b10 740, b11 842 nm.
# ==============================================================================
INDEX_EXPR <- rlang::exprs(
  NDVI=      (b11-b7)/(b11+b7),              NDRE=      (b11-b9)/(b11+b9),
  SAVI=      1.5*(b11-b7)/(b11+b7+0.5),      DVI=       b11-b7,
  EVI=       2.5*(b11-b7)/(b11+6*b7-7.5*b2+1),
  EVI2=      2.5*(b11-b7)/(b11+2.4*b7+1),    RDVI=      (b11-b7)/sqrt(b11+b7),
  MSR=       (b11/b7-1)/(sqrt(b11/b7)+1),
  MSAVI=     (2*b11+1-sqrt((2*b11+1)^2-8*(b11-b7)))/2,
  SR=        b11/b7,
  MTVI1=     1.2*(1.2*(b11-b4)-2.5*(b7-b4)),
  MTVI2=     (1.5*(1.2*(b11-b4)-2.5*(b7-b4)))/
             sqrt((2*b11+1)^2-(6*b11-5*sqrt(b7))-0.5),
  MTVI2_den= sqrt((2*b11+1)^2-(6*b11-5*sqrt(b7))-0.5),
  MCARI1=    1.2*(2.5*(b11-b7)-1.3*(b11-b4)),
  SR528_574= b3/b4,                          SR800_635= b11/b6,
  NDVI2=     (b10-b8)/(b10+b8),              NDRE3=     (b11-b10)/(b11+b10),
  NDRE_705=  (b11-b8)/(b11+b8),              NDVI_RE=   (b11-b6)/(b11+b6),
  RE_Slope=  (b10-b8)/35,                    VOG1=      b10/b9,
  GM1=       b10/b4,
  MCARI=     ((b8-b7)-0.2*(b8-b4))*(b8/b7),
  MCARI3=    ((b10-b8)-0.2*(b10-b4))*(b10/b8),
  CARI=      (b8-b7)-0.2*(b8-b4),
  TCARI=     3*((b8-b7)-0.2*(b8-b4)*(b8/b7)),
  OSAVI=     (1+0.16)*(b11-b7)/(b11+b7+0.16),
  TCA_OSA=   TCARI/OSAVI,
  TVI=       0.5*(120*(b10-b4)-200*(b7-b4)),
  PSSR_c=    b11/b2,                         PSNDc=     (b11-b2)/(b11+b2),
  CAR=       b8/b10,                         RARSa=     b7/b8,
  RARS=      b10/b3,                         RGI=       b8/b4,
  RR=        1/b8,                           Carter1=   b8/b1,
  ARI=       (1/b4)-(1/b8),                  ARI2=      (1/b4)-(1/b9),
  CRI550_515=(1/b3)-(1/b4),                  CRI700_515=(1/b3)-(1/b8),
  DCabCxc=   b7/(b4*3*b8),                   DNIRCabCxc=b11/(b4*b8),
  SIPI=      (b11-b1)/(b11+b7),              PSRI=      (b7-b2)/b10,
  NPCI=      (b7-b1)/(b7+b1),                PBI=       b11/b4,
  CVI=       (b11*b8)/(b4^2),                CI=        b10/b8,
  PRI=       (b3-b4)/(b3+b4),                PRIm3=     (b7-b3)/(b7+b3),
  PRIm4=     (b4-b3-b7)/(b4+b3+b7),
  PRIn=      ((b4-b3)/(b4+b3))/(((b11-b7)/sqrt(b11+b7))*(b8/b7)),
  PRI_CI=    ((b4-b3)/(b4+b3))*((b10/b8)-1),
  CCI=       (b3-b6)/(b3+b6),                LIC3=      b1/b10,
  GCC=       (b3+b4)/((b6+b7)+(b3+b4)+(b1+b2)),
  GNDVI=     (b11-b4)/(b11+b4),              GNDVI2=    (b11-b3)/(b11+b3),
  GRVI=      b11/b3,                         NGRDI=     (b4-b7)/(b4+b7),
  VARI=      (b4-b7)/(b4+b7-b2),
  GLI=       (2*b4-b7-b2)/(2*b4+b7+b2),
  DWSI4=     b4/b7,                          R_index=   b8/b7,
  B_index=   b1/b2,                          BGI2=      b1/b4,
  BRI=       b1/b8,                          CN1=       (b2-b7)/(b2+b7),
  CN2=       b2/b7,                          SRPI=      b1/b7,
  LIC4=      b1/b3,                          LIC5=      b8/b3,
  NDWI=      (b3-b11)/(b3+b11),              HI=        (b3-b8)/(b3+b8)-0.5*b8,
  CLS=       (b8-b4)/(b8+b4)-b10,            PMI=       (b3-b4)/(b3+b4)+b9,
  SBRI=      (b4-b3)/(b4+b3)+0.5*b8,
  VIN_B8=    b11/b8,                         VIN_B9=    b11/b9,
  VIN_B10=   b11/b10
)


# ==============================================================================
# 4. INPUTS
# ==============================================================================
cat("\n[4] Loading inputs...\n")

split_severity <- function(crowns) {
  crowns %>%
    mutate(raw        = trimws(as.character(DOTHI___IN)),
           is_dead    = !is.na(raw) & tolower(raw) == "dead",
           DOTHI___IN = case_when(!is.na(raw) & tolower(raw) == "dead" ~ NA_real_,
                                  TRUE ~ suppressWarnings(as.numeric(raw)))) %>%
    select(-raw)
}

crowns_s1 <- st_read(CROWNS_SITE1, quiet = TRUE) %>%
  st_transform(4326) %>% mutate(plot_id = paste0("S1_", Pset_12)) %>% split_severity()
crowns_s2 <- suppressWarnings(st_read(CROWNS_SITE2, quiet = TRUE)) %>%
  st_make_valid() %>% st_transform(4326) %>%
  rename(DOTHI___IN = Dothi___In, Label_12 = FieldLabel) %>%
  mutate(plot_id = paste0("S2_", Pset)) %>% split_severity()

raster_s1 <- rast(resolve_input(NAMES_MASKED_SITE1, "masked Trial 1 raster"))
raster_s2 <- rast(resolve_input(NAMES_MASKED_SITE2, "masked Trial 2 raster"))

FEATURE_LIST   <- resolve_input(NAMES_FEATURE_LIST,   "selected-feature list")
CROWN_FEATURES <- resolve_input(NAMES_CROWN_FEATURES, "crown feature table")
TEST_PRED_OF   <- resolve_input(NAMES_TEST_PRED_OF,   "saved test predictions")

features <- clean_names(readLines(FEATURE_LIST))

crown_features <- read.csv(CROWN_FEATURES, check.names = FALSE)
names(crown_features) <- clean_names(names(crown_features))
crown_features$Label_12 <- as.character(crown_features$Label_12)

pred_of <- read.csv(TEST_PRED_OF)
test_labels <- as.character(pred_of$Label_12)

model_data <- crown_features %>%
  select(Label_12, DOTHI___IN, all_of(features)) %>%
  filter(!is.na(DOTHI___IN), complete.cases(.))
train_labels <- setdiff(model_data$Label_12, test_labels)
cat(sprintf("  %d features | %d modelling crowns = %d train + %d test\n",
            length(features), nrow(model_data), length(train_labels), length(test_labels)))

## Which indices, bands and statistics the selected features actually require.
STAT_SUFFIX  <- "mean|sd|med|p10|p90|cv"
feature_base <- sub(paste0("_(", STAT_SUFFIX, ")$"), "", features)
feature_stat <- sub(paste0("^.*_(", STAT_SUFFIX, ")$"), "\\1", features)
need_stat    <- unique(feature_stat)
need_index   <- intersect(names(INDEX_EXPR), unique(feature_base))
if ("TCA_OSA" %in% need_index) need_index <- union(need_index, c("TCARI", "OSAVI"))
need_index   <- names(INDEX_EXPR)[names(INDEX_EXPR) %in% need_index]  # definition order
need_band    <- intersect(paste0("b", 1:11), unique(feature_base))
cat(sprintf("  evaluating %d indices and %d raw bands, statistics: %s\n",
            length(need_index), length(need_band), paste(need_stat, collapse = ", ")))


# ==============================================================================
# 5. CROWN SUMMARIES UNDER A PERTURBED BOUNDARY
# ==============================================================================
extract_pixels <- function(raster_obj, polygons, site_label) {
  px <- terra::extract(raster_obj, vect(polygons), ID = TRUE)
  colnames(px)[2:12] <- paste0("b", 1:11)
  px$Label_12   <- polygons$Label_12[px$ID]
  px$DOTHI___IN <- polygons$DOTHI___IN[px$ID]
  px$is_dead    <- polygons$is_dead[px$ID]
  px$plot_id    <- polygons$plot_id[px$ID]
  px$site       <- site_label
  px
}

summarise_crowns <- function(px) {
  cols   <- c(need_index, need_band)
  keys   <- c("site", "plot_id", "Label_12", "DOTHI___IN")
  grouped<- px %>% group_by(across(all_of(keys)))
  parts  <- list(grouped %>% summarise(n_pixels = n(), .groups = "drop"))
  add <- function(fn, suffix)
    grouped %>% summarise(across(all_of(cols), fn, .names = paste0("{.col}_", suffix)),
                          .groups = "drop")
  if ("mean" %in% need_stat) parts <- c(parts, list(add(~mean(.x, na.rm = TRUE),   "mean")))
  if ("sd"   %in% need_stat) parts <- c(parts, list(add(~sd(.x, na.rm = TRUE),     "sd")))
  if ("med"  %in% need_stat) parts <- c(parts, list(add(~median(.x, na.rm = TRUE), "med")))
  if ("p10"  %in% need_stat) parts <- c(parts, list(add(~quantile(.x, QUANTILES[["p10"]], na.rm = TRUE), "p10")))
  if ("p90"  %in% need_stat) parts <- c(parts, list(add(~quantile(.x, QUANTILES[["p90"]], na.rm = TRUE), "p90")))
  if ("cv"   %in% need_stat) parts <- c(parts, list(add(~safe_cv(.x),              "cv")))
  out <- Reduce(function(a, b) left_join(a, b, by = keys), parts)
  names(out) <- clean_names(names(out))
  out
}

summaries_for <- function(buffer_px, label) {
  cat(sprintf("\n  --- %s ---\n", label))
  per_site <- list()
  for (site_code in c("Site1", "Site2")) {
    polygons <- if (site_code == "Site1") crowns_s1 else crowns_s2
    raster   <- if (site_code == "Site1") raster_s1 else raster_s2
    distance <- buffer_px * GSD_M[[site_code]]

    if (distance != 0) {
      buffered <- st_transform(polygons, 2193) %>% st_buffer(distance) %>% st_make_valid()
      lost     <- sum(st_is_empty(buffered))
      polygons <- buffered[!st_is_empty(buffered), ] %>% st_transform(4326)
      cat(sprintf("    %s buffer %+.4f m, %d crowns lost\n", site_code, distance, lost))
    } else {
      cat(sprintf("    %s original boundaries\n", site_code))
    }

    px <- extract_pixels(raster, polygons, site_code) %>%
      filter(!is.na(b1), !is_dead, !is.na(DOTHI___IN)) %>%
      mutate(!!!INDEX_EXPR[need_index]) %>%
      mutate(across(all_of(need_index), ~ifelse(is.infinite(.x) | is.nan(.x), NA, .x)))

    per_site[[site_code]] <- summarise_crowns(px)
    rm(px); gc(verbose = FALSE)
  }
  bind_rows(per_site) %>%
    mutate(boundary = label, Label_12 = as.character(Label_12))
}

cat("\n[5] Building boundary variants...\n")
variants <- Map(summaries_for, BUFFERS_PX, BUFFER_LABELS)
names(variants) <- BUFFER_LABELS

## Erosion can extinguish a small crown, so comparisons are restricted to the
## crowns that survive every perturbation.
common <- Reduce(intersect, lapply(variants, function(d) d$Label_12))
variants <- lapply(variants, function(d) d %>% filter(Label_12 %in% common) %>% arrange(Label_12))
cat(sprintf("\n  crowns present under every boundary: %d\n", length(common)))

## The unperturbed variant should reproduce the feature table the models were
## built from. If it does not, the perturbations cannot be interpreted.
recheck <- variants[["original"]] %>% select(Label_12, all_of(features)) %>%
  inner_join(crown_features %>% select(Label_12, all_of(features)),
             by = "Label_12", suffix = c("_new", "_ref"))
recheck_r <- sapply(features, function(f) {
  a <- recheck[[paste0(f, "_ref")]]; b <- recheck[[paste0(f, "_new")]]
  g <- is.finite(a) & is.finite(b)
  suppressWarnings(cor(a[g], b[g]))
})
report(min(recheck_r, na.rm = TRUE) > 0.9999,
       sprintf("unperturbed boundaries reproduce the reference features (min r = %.6f)",
               min(recheck_r, na.rm = TRUE)))


# ==============================================================================
# 6. HOW LARGE WAS THE PERTURBATION?
# Reported as the realised change in pixel count rather than the nominal buffer,
# because the two trials differ in ground sampling distance and in crown size.
# ==============================================================================
cat("\n[6] Realised perturbation...\n")
pixel_counts <- bind_rows(lapply(names(variants), function(v)
  variants[[v]] %>% select(boundary, Label_12, site, n_pixels)))
reference_px <- pixel_counts %>% filter(boundary == "original") %>%
  select(Label_12, n_reference = n_pixels)
perturbation <- pixel_counts %>% filter(boundary != "original") %>%
  left_join(reference_px, by = "Label_12") %>%
  group_by(boundary, site) %>%
  summarise(median_pct_change = round(median(100 * (n_pixels - n_reference) / n_reference), 1),
            .groups = "drop") %>%
  pivot_wider(names_from = site, values_from = median_pct_change)
print(as.data.frame(perturbation), row.names = FALSE)
write.csv(perturbation, file.path(DIR_ANALYSIS, "boundary_perturbation.csv"),
          row.names = FALSE)


# ==============================================================================
# 7. DO THE FEATURE VALUES MOVE?
#
# Feature selection is not repeated: this asks whether the values of the 40
# features already selected are stable, not whether the same features would be
# selected again.
# ==============================================================================
cat("\n[7] Feature stability...\n")
reference <- variants[["original"]]
stability <- lapply(setdiff(names(variants), "original"), function(v) {
  alternative <- variants[[v]]
  do.call(rbind, lapply(features, function(f) {
    a <- reference[[f]]; b <- alternative[[f]]
    g <- is.finite(a) & is.finite(b)
    data.frame(boundary = v, feature = f, n = sum(g),
               pearson  = round(suppressWarnings(cor(a[g], b[g])), 4),
               spearman = round(suppressWarnings(cor(a[g], b[g], method = "spearman")), 4),
               ## Relative change is summarised by its median only. Several
               ## indices cross zero, so the mean and the maximum of a ratio to
               ## the original value are dominated by near-zero denominators
               ## and carry no information about stability.
               median_rel_change_pct = round(
                 100 * median(abs(b[g] - a[g]) / pmax(abs(a[g]), 1e-9)), 2))
  }))
}) %>% bind_rows()

stability_summary <- stability %>% group_by(boundary) %>%
  summarise(features = n(),
            min_pearson    = round(min(pearson, na.rm = TRUE), 3),
            median_pearson = round(median(pearson, na.rm = TRUE), 3),
            n_below_0.95   = sum(pearson < 0.95, na.rm = TRUE),
            median_rel_change_pct = round(median(median_rel_change_pct, na.rm = TRUE), 2),
            .groups = "drop")
print(as.data.frame(stability_summary), row.names = FALSE)

write.csv(stability,         file.path(DIR_ANALYSIS, "boundary_feature_stability.csv"),         row.names = FALSE)
write.csv(stability_summary, file.path(DIR_ANALYSIS, "boundary_feature_stability_summary.csv"), row.names = FALSE)


# ==============================================================================
# 8. DO THE PREDICTIONS MOVE?
#
# The model is fitted once, on the original training features, and is not
# refitted. Only the test input is perturbed.
# ==============================================================================
if (RUN_PREDICT) {
cat("\n[8] Prediction stability...\n")

train_levels <- sort(unique(model_data$DOTHI___IN[model_data$Label_12 %in% train_labels]))
fabof_frame <- function(d) {
  d <- d %>% select(Label_12, DOTHI___IN, all_of(features)) %>% filter(complete.cases(.))
  d$health_class_fac <- factor(as.character(d$DOTHI___IN),
                               levels = as.character(train_levels))
  as.data.frame(d)
}

set.seed(SEED)
model <- fabOF::fabOF(
  formula = as.formula(paste("health_class_fac ~", paste(features, collapse = " + "))),
  data    = fabof_frame(crown_features %>% filter(Label_12 %in% train_labels)),
  importance = FALSE,
  ranger.control = list(num.trees = FABOF_TREES, mtry = FABOF_MTRY, num.threads = 0))

predictions <- lapply(names(variants), function(v) {
  test <- fabof_frame(variants[[v]] %>% filter(Label_12 %in% test_labels))
  data.frame(Label_12 = test$Label_12, observed = test$DOTHI___IN,
             predicted = as.numeric(as.character(predict(model, newdata = test))),
             boundary = v, stringsAsFactors = FALSE)
}) %>% bind_rows()

common_test <- Reduce(intersect, split(predictions$Label_12, predictions$boundary))
predictions <- predictions %>% filter(Label_12 %in% common_test)
cat(sprintf("  test crowns scored under every boundary: %d of %d\n",
            length(common_test), length(test_labels)))

performance <- lapply(names(variants), function(v) {
  d <- predictions %>% filter(boundary == v)
  score(d$observed, d$predicted, v)
}) %>% bind_rows()
print(as.data.frame(performance), row.names = FALSE)

## Agreement with the predictions made from the original boundaries. Exact
## agreement is not reported: the paper argues that it is the wrong lens for a
## 19-level ordinal scale, so the shift is given in severity points instead.
wide <- predictions %>% select(Label_12, boundary, predicted) %>%
  pivot_wider(names_from = boundary, values_from = predicted)
agreement <- lapply(setdiff(names(variants), "original"), function(v)
  data.frame(boundary = v,
             mean_shift_points = round(mean(abs(wide[[v]] - wide$original)), 2),
             within_5_pct      = round(mean(abs(wide[[v]] - wide$original) <= 5) * 100, 2),
             max_shift_points  = max(abs(wide[[v]] - wide$original)))) %>% bind_rows()
print(agreement, row.names = FALSE)

write.csv(performance, file.path(DIR_ANALYSIS, "boundary_performance.csv"), row.names = FALSE)
write.csv(agreement,   file.path(DIR_ANALYSIS, "boundary_agreement.csv"),   row.names = FALSE)
}

cat("\n[done] Crown boundary sensitivity complete.\n")
