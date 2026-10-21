# ==============================================================================
# Spatial validation of the DNB severity models
#
# Companion analysis: how much of the reported predictive performance survives
#                     when training and test crowns are forced apart in space
#
# Checks    1. buffered evaluation of the reported plot-level partition
#           2. spatially blocked cross-validation with exclusion radii
#           3. matched random-removal control separating the spatial effect
#              from the cost of the smaller training set
# Outputs   Appendix Table 4 and the figures quoted in Sections 3.2 and 4
#
# Reproducibility notes
#   - Edit the CONFIGURATION block below; nothing else needs changing.
#   - Reads the crown feature table, the selected-feature list and the saved
#     test predictions written by the main analysis pipeline. Run that first.
#   - No model is refitted for check 1; the saved predictions are re-scored on
#     progressively smaller subsets of the test set.
#   - Checks 2 and 3 refit fabOF only. Hyperparameters are fixed at the values
#     used for the reported model and are not re-tuned, and feature selection
#     is not repeated within folds.
#   - All random operations are seeded via SEED.
#   - Quadratic weighted kappa is computed on the full severity scale, so that
#     a class absent from a subset cannot compress the weighting.
# ==============================================================================


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

## ---- Paths -----------------------------------------------------------------
DIR_BASE     <- "J:/RPBC_Dothi_Analysis"          # <- edit this
DIR_REF      <- file.path(DIR_BASE, "Reference")
DIR_ANALYSIS <- file.path(DIR_BASE, "Analysis")

## Crown polygons. Not distributed with this repository; see the data
## availability statement in the README.
CROWNS_SITE1 <- file.path(DIR_REF, "Crown_Health_Reference_2026.shp")
CROWNS_SITE2 <- file.path(DIR_REF, "Crown_Health_Reference_2026_site2_clean.shp")

## Artefacts written by the main pipeline, looked for in DIR_ANALYSIS. Output
## names have changed between pipeline versions, so each entry lists the names
## in use and the first that exists is taken. Edit a list if yours differs.
NAMES_FEATURE_LIST   <- c("selected_features.txt",
                          "combined_selected_features_list.txt")
NAMES_CROWN_FEATURES <- c("crown_features_for_genetic_analysis.csv",
                          "crown_summaries_ext_for_genetics.csv")
NAMES_TEST_PRED_OF   <- c("ordfor_19-class_predictions.csv",
                          "ordfor_finegrained_predictions.csv")
NAMES_TEST_PRED_FAB  <- c("fabof_19-class_predictions.csv",
                          "fabof_finegrained_predictions.csv")
NAMES_REPORTED       <- c("results_combined_trials.csv",
                          "fabof_vs_ordinalforest_comparison.csv")

## ---- Run control -----------------------------------------------------------
SEED           <- 42
RUN_BLOCKED_CV <- TRUE     # checks 2 and 3; ~25 fabOF fits, roughly 20 minutes
RUN_CONTROL    <- TRUE     # check 3 only; requires RUN_BLOCKED_CV

## ---- Geometry --------------------------------------------------------------
CRS_METRIC     <- 2193         # NZGD2000 / NZTM2000, for distances in metres
BUFFERS_M      <- c(0, 2.5, 5, 7.5, 10, 12.5, 15)   # check 1, separation steps
MIN_TEST_CROWNS<- 40           # skip a buffer retaining fewer crowns than this

## ---- Cross-validation ------------------------------------------------------
BLOCK_SIZE_M   <- 150          # square blocks, several times the autocorrelation range
EXCLUSION_M    <- c(0, 20, 40) # training crowns within this distance of a held-out
                               # block are deleted before refitting
K_FOLDS        <- 5
MIN_TRAIN_KEPT <- 300          # abandon a fold retaining fewer training crowns
CONTROL_REPS   <- 2            # random draws per fold per exclusion radius

## ---- Model -----------------------------------------------------------------
## As used for the reported 19-class fabOF model. Not re-tuned here.
FABOF_TREES    <- 5000
FABOF_MTRY     <- 10
SEVERITY_SCALE <- seq(5, 95, by = 5)   # the field protocol's 19 ordinal grades

## ---- Pseudo-plots ----------------------------------------------------------
## Unscored crowns are grouped under a placeholder plot code that spans the
## whole trial. Such groups are not plots and must not define a CV block.
MAX_PLOT_WIDTH_M <- 40


# ==============================================================================
# 2. LIBRARIES AND HELPERS
# ==============================================================================
library(sf); library(dplyr); library(tidyr); library(caret)
if (RUN_BLOCKED_CV) library(fabOF)

clean_names <- function(x) gsub("[^a-zA-Z0-9_]", "_", x)
report      <- function(pass, msg) cat(sprintf("  [%s] %s\n",
                                               ifelse(pass, "PASS", "FAIL"), msg))

## Take the first of the candidate filenames that exists, and say plainly what
## was missing if none do. Optional inputs return NA rather than stopping.
resolve_input <- function(candidates, what, required = TRUE) {
  found <- file.path(DIR_ANALYSIS, candidates)
  found <- found[file.exists(found)]
  if (length(found)) return(found[1])
  message <- sprintf("Cannot find the %s in %s.\n  Looked for: %s\n  Edit the candidate list in the CONFIGURATION block.",
                     what, DIR_ANALYSIS, paste(candidates, collapse = ", "))
  if (required) stop(message, call. = FALSE)
  cat("  note:", what, "not found, that part of the analysis will be skipped\n")
  NA_character_
}

## Quadratic weighted kappa on a fixed level set. Using the levels present in
## each subset would let a missing intermediate grade shorten the ordinal
## distance between its neighbours and inflate the statistic.
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
  data.frame(
    design    = label,
    n_crowns  = length(obs),
    n_classes = length(unique(obs)),
    bal_acc   = round(mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE) * 100, 2),
    qwk       = round(qwk(obs, pred), 3),
    mae       = round(mean(abs(err)), 2),
    bias      = round(mean(err), 2),
    w5        = round(mean(abs(err) <=  5) * 100, 2),
    w10       = round(mean(abs(err) <= 10) * 100, 2))
}


# ==============================================================================
# 3. INPUTS
# ==============================================================================
cat("\n[3] Loading inputs...\n")

## Severity is recorded as a percentage, with dead crowns flagged in the same
## field. Separate the two so dead crowns do not enter as numeric zeroes.
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

crowns <- bind_rows(crowns_s1 %>% select(Label_12, plot_id),
                    crowns_s2 %>% select(Label_12, plot_id)) %>%
  mutate(Label_12 = as.character(Label_12),
         site     = substr(plot_id, 1, 2)) %>%
  st_transform(CRS_METRIC)

FEATURE_LIST     <- resolve_input(NAMES_FEATURE_LIST,   "selected-feature list")
CROWN_FEATURES   <- resolve_input(NAMES_CROWN_FEATURES, "crown feature table")
TEST_PRED_OF     <- resolve_input(NAMES_TEST_PRED_OF,   "saved OF test predictions")
TEST_PRED_FABOF  <- resolve_input(NAMES_TEST_PRED_FAB,  "saved fabOF test predictions",
                                  required = FALSE)
REPORTED_METRICS <- resolve_input(NAMES_REPORTED,       "reported performance table")

features <- clean_names(readLines(FEATURE_LIST))

crown_features <- read.csv(CROWN_FEATURES, check.names = FALSE)
names(crown_features) <- clean_names(names(crown_features))
crown_features$Label_12 <- as.character(crown_features$Label_12)

pred_of <- read.csv(TEST_PRED_OF)
pred_of$Label_12 <- as.character(pred_of$Label_12)
has_fabof <- !is.na(TEST_PRED_FABOF)
if (has_fabof) {
  pred_fabof <- read.csv(TEST_PRED_FABOF)
  pred_fabof$Label_12 <- as.character(pred_fabof$Label_12)
}

## The modelling data: crowns with a severity grade and a complete feature set.
model_data <- crown_features %>%
  select(Label_12, plot_id, DOTHI___IN, all_of(features)) %>%
  filter(!is.na(DOTHI___IN), complete.cases(.)) %>%
  mutate(y = as.numeric(DOTHI___IN),
         health_class = factor(y, levels = sort(unique(y)), ordered = TRUE))

test_labels  <- pred_of$Label_12
train_labels <- setdiff(model_data$Label_12, test_labels)
cat(sprintf("  %d features | %d modelling crowns = %d train + %d test\n",
            length(features), nrow(model_data),
            length(train_labels), length(test_labels)))
report(length(train_labels) + length(test_labels) == nrow(model_data),
       "train and test partition the modelling data exactly")


# ==============================================================================
# 4. CHECK 1 - BUFFERED EVALUATION OF THE REPORTED PARTITION
#
# Separation is increased by discarding test crowns that lie close to training
# data. The training set never changes. No model is refitted: the saved
# predictions are simply re-scored on each surviving subset.
# ==============================================================================
cat("\n[4] Buffered evaluation...\n")

train_geom <- crowns %>% filter(Label_12 %in% train_labels)
test_geom  <- suppressWarnings(st_centroid(crowns %>% filter(Label_12 %in% test_labels)))
report(nrow(train_geom) == length(train_labels) && nrow(test_geom) == length(test_labels),
       "every modelling crown matched to a polygon")

nearest <- st_nearest_feature(test_geom, train_geom)
test_geom$d_train <- as.numeric(st_distance(test_geom, train_geom[nearest, ],
                                            by_element = TRUE))

cat("  distance from test crown to nearest training crown (m):\n")
print(round(summary(test_geom$d_train), 2))
cat(sprintf("  ceiling: no test crown lies further than %.2f m from training data.\n",
            max(test_geom$d_train)))
cat("  Training and test plots are interleaved, so this check cannot reach the\n")
cat("  autocorrelation range. Check 2 exists for that reason.\n")

buffered <- lapply(BUFFERS_M, function(b) {
  keep <- test_geom$Label_12[test_geom$d_train > b]
  if (length(keep) < MIN_TEST_CROWNS) return(NULL)
  label <- sprintf("buffer %.1f m", b)
  out <- score(pred_of$True_class[pred_of$Label_12 %in% keep],
               pred_of$Predicted_class[pred_of$Label_12 %in% keep],
               label) %>% mutate(model = "OF", min_sep_m = b, .after = design)
  if (has_fabof) {
    k <- pred_fabof$Label_12 %in% keep
    out <- bind_rows(out,
      score(pred_fabof$True_class[k], pred_fabof$Predicted_class[k], label) %>%
        mutate(model = "fabOF", min_sep_m = b, .after = design))
  }
  out
}) %>% bind_rows()
print(as.data.frame(buffered), row.names = FALSE)

## Buffer zero is the reported test set, so it must reproduce the published
## metrics. If it does not, something upstream has changed.
## Column names in the reported table also differ between pipeline versions
## (Balanced_Acc against balanced_acc, Within_pm10pct against within_pm10 and
## so on), so they are normalised before matching.
reported <- read.csv(REPORTED_METRICS)
names(reported) <- sub("pct$", "", gsub("_", "", tolower(names(reported))))
reported <- reported %>% filter(grepl("19", granularity))
verify <- buffered %>% filter(min_sep_m == 0) %>%
  select(model, bal_acc, mae, w10) %>%
  left_join(reported %>% transmute(model   = algorithm,
                                   ref_bal = balancedacc * 100,
                                   ref_mae = mae,
                                   ref_w10 = withinpm10 * 100), by = "model")
report(isTRUE(all(abs(verify$bal_acc - verify$ref_bal) < 0.05,
                  abs(verify$mae     - verify$ref_mae) < 0.05,
                  abs(verify$w10     - verify$ref_w10) < 0.05)),
       "buffer zero reproduces the reported metrics")

write.csv(buffered, file.path(DIR_ANALYSIS, "spatial_buffered_evaluation.csv"),
          row.names = FALSE)


# ==============================================================================
# 5. CHECK 2 AND 3 - BUFFERED SPATIALLY BLOCKED CROSS-VALIDATION
#
# Separation is increased from the other side: whole blocks of plots are held
# out, then every training crown within the exclusion radius of a held-out
# crown is deleted before refitting. The test set is identical across exclusion
# radii, so the only thing that changes is the training data.
#
# That exclusion removes the nearest, and therefore most informative, training
# crowns, so some of the resulting loss is simply a smaller training set. The
# control arm removes the same number of crowns at random to separate the two.
# ==============================================================================
if (RUN_BLOCKED_CV) {
cat("\n[5] Buffered spatially blocked cross-validation...\n")

crown_centroids <- suppressWarnings(st_centroid(
  crowns %>% filter(!is.na(Label_12), Label_12 %in% model_data$Label_12)))

## Drop placeholder plot codes before blocking.
plot_width <- crown_centroids %>% group_by(plot_id) %>% summarise(.groups = "drop") %>%
  mutate(width_m = sapply(st_geometry(.), function(g) {
    bb <- st_bbox(g); max(bb[3] - bb[1], bb[4] - bb[2]) })) %>% st_drop_geometry()
pseudo_plots <- plot_width$plot_id[plot_width$width_m > MAX_PLOT_WIDTH_M]
if (length(pseudo_plots))
  cat("  excluding pseudo-plots:", paste(pseudo_plots, collapse = ", "), "\n")

cv_centroids <- crown_centroids %>% filter(!plot_id %in% pseudo_plots)
cv_data      <- model_data %>% filter(Label_12 %in% cv_centroids$Label_12)
cat(sprintf("  crowns available to cross-validation: %d\n", nrow(cv_data)))

## Blocks are built per trial so that both contribute to every fold, and each
## plot is assigned whole to its modal block so no plot straddles a fold.
make_folds <- function(block_m, k = K_FOLDS, seed = SEED) {
  set.seed(seed)
  per_site <- function(site_code) {
    pts <- cv_centroids %>% filter(site == site_code)
    if (nrow(pts) == 0) return(NULL)
    grid <- st_make_grid(pts, cellsize = block_m, square = TRUE)
    cell <- as.integer(sapply(st_intersects(pts, grid, sparse = TRUE), function(z) z[1]))
    data.frame(plot_id = pts$plot_id, block = paste0(site_code, "_B", cell),
               stringsAsFactors = FALSE)
  }
  assigned <- bind_rows(per_site("S1"), per_site("S2")) %>% filter(!is.na(block))
  per_plot <- assigned %>% count(plot_id, block) %>% group_by(plot_id) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>% select(plot_id, block)
  shuffled <- sample(unique(per_plot$block))
  per_plot$fold <- setNames(rep_len(seq_len(k), length(shuffled)), shuffled)[per_plot$block]
  per_plot
}

fabof_frame <- function(d) {
  d <- d %>% select(Label_12, y, health_class, all_of(features)) %>% filter(complete.cases(.))
  d$health_class_fac <- factor(as.character(d$health_class),
                               levels = levels(d$health_class), ordered = FALSE)
  as.data.frame(d)
}
fabof_formula <- as.formula(paste("health_class_fac ~", paste(features, collapse = " + ")))

fit_and_score <- function(frame, train_keep, test_keep, label, seed = SEED) {
  set.seed(seed)
  model <- fabOF::fabOF(
    formula = fabof_formula,
    data    = fabof_frame(frame %>% filter(Label_12 %in% train_keep)),
    importance = FALSE,
    ranger.control = list(num.trees = FABOF_TREES, mtry = FABOF_MTRY, num.threads = 0))
  test <- fabof_frame(frame %>% filter(Label_12 %in% test_keep))
  score(test$y, as.numeric(as.character(predict(model, newdata = test))), label)
}

folds  <- make_folds(BLOCK_SIZE_M)
framed <- cv_data %>% left_join(folds, by = "plot_id") %>% filter(!is.na(fold))
geom   <- cv_centroids %>% select(Label_12)

cv_results <- list()
for (f in seq_len(K_FOLDS)) {
  held_out  <- framed$Label_12[framed$fold == f]
  available <- framed$Label_12[framed$fold != f]
  g_test    <- geom %>% filter(Label_12 %in% held_out)
  g_train   <- geom %>% filter(Label_12 %in% available)
  if (nrow(g_test) < MIN_TEST_CROWNS) next

  for (radius in EXCLUSION_M) {
    too_close <- if (radius > 0)
      lengths(st_is_within_distance(g_train, g_test, dist = radius)) > 0
      else rep(FALSE, nrow(g_train))
    keep <- g_train$Label_12[!too_close]
    if (length(keep) < MIN_TRAIN_KEPT) next

    kept_geom  <- g_train %>% filter(Label_12 %in% keep)
    separation <- as.numeric(st_distance(
      g_test, kept_geom[st_nearest_feature(g_test, kept_geom), ], by_element = TRUE))

    label <- sprintf("block %dm / exclusion %dm", BLOCK_SIZE_M, radius)
    cat(sprintf("  fold %d  %-28s train %5d (%4.1f%% dropped)  separation min %5.1f median %5.1f m\n",
                f, label, length(keep), 100 * mean(too_close),
                min(separation), median(separation)))

    cv_results[[length(cv_results) + 1]] <-
      fit_and_score(framed, keep, held_out, label) %>%
      mutate(arm = "spatial", fold = f, exclusion_m = radius,
             n_train = length(keep),
             train_dropped_pct = round(100 * mean(too_close), 1),
             min_sep_m    = round(min(separation), 1),
             median_sep_m = round(median(separation), 1), .after = design)

    ## Control: same fold, same test crowns, same number of training crowns
    ## removed, but chosen at random rather than by distance.
    if (RUN_CONTROL && radius > 0) {
      for (rep in seq_len(CONTROL_REPS)) {
        set.seed(SEED + 1000 * f + 10 * radius + rep)
        random_keep <- sample(g_train$Label_12, length(keep))
        cv_results[[length(cv_results) + 1]] <-
          fit_and_score(framed, random_keep, held_out, label) %>%
          mutate(arm = "random", fold = f, exclusion_m = radius,
                 n_train = length(random_keep),
                 train_dropped_pct = round(100 * mean(too_close), 1),
                 min_sep_m = NA_real_, median_sep_m = NA_real_, .after = design)
      }
    }
  }
}
cv_results <- bind_rows(cv_results)

cat("\n  Cross-validation summary\n")
cv_summary <- cv_results %>% group_by(exclusion_m, arm) %>%
  summarise(fits = n(), mean_n_train = round(mean(n_train)),
            min_sep_m    = ifelse(all(is.na(min_sep_m)), NA, round(min(min_sep_m, na.rm = TRUE), 1)),
            median_sep_m = ifelse(all(is.na(median_sep_m)), NA, round(median(median_sep_m, na.rm = TRUE), 1)),
            across(c(bal_acc, qwk, mae, w5, w10),
                   list(mean = ~round(mean(.x), 2), sd = ~round(sd(.x), 2))),
            .groups = "drop") %>% arrange(exclusion_m, arm)
print(as.data.frame(cv_summary), row.names = FALSE)

## Decomposition: how much of the loss under spatial exclusion is simply the
## cost of a smaller training set, and how much is spatial dependence.
if (RUN_CONTROL) {
  cat("\n  Decomposition, relative to no exclusion\n")
  baseline <- cv_results %>% filter(exclusion_m == 0) %>%
    select(fold, mae_base = mae, w10_base = w10)
  decomposition <- cv_results %>% filter(exclusion_m > 0) %>%
    group_by(exclusion_m, fold, arm) %>%
    summarise(mae = mean(mae), w10 = mean(w10), .groups = "drop") %>%
    pivot_wider(names_from = arm, values_from = c(mae, w10)) %>%
    left_join(baseline, by = "fold") %>%
    group_by(exclusion_m) %>%
    summarise(training_size_mae = round(mean(mae_random  - mae_base),   2),
              spatial_mae       = round(mean(mae_spatial - mae_random), 2),
              training_size_w10 = round(mean(w10_random  - w10_base),   2),
              spatial_w10       = round(mean(w10_spatial - w10_random), 2),
              .groups = "drop")
  print(as.data.frame(decomposition), row.names = FALSE)
  cat("  A spatial term near zero would mean the loss was sample size alone.\n")
  write.csv(decomposition,
            file.path(DIR_ANALYSIS, "spatial_decomposition.csv"), row.names = FALSE)
}

write.csv(cv_results, file.path(DIR_ANALYSIS, "spatial_blocked_cv_folds.csv"),   row.names = FALSE)
write.csv(cv_summary, file.path(DIR_ANALYSIS, "spatial_blocked_cv_summary.csv"), row.names = FALSE)
}

cat("\n[done] Spatial validation complete.\n")
