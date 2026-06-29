# This document is a continuation of 02_SCAM_only with a reconstruction to finalize. 
# This code was the first attempt to create a full reconstruction 
# and was further improved in python code, in other folders



# =============================================================================
#   NPi RECONSTRUCTION FORMULA — FINAL MODEL
#   Derived via SCAM approximation with piecewise linear (hinge) terms
#
#   Formula:
#   NPi = 5.771
#       + 0.359 * z_pupil_size  - 0.081 * max(z_pupil_size + 4, 0)
#       - 0.283 * z_pupil_min   - 0.123 * max(z_pupil_min + 4, 0)
#       - 0.235 * z_ch          - 0.046 * z_ch^2
#       + 0.212 * z_max_const_velocity - 0.206 * max(z_max_const_velocity + 3, 0)
#       - 0.049 * z_latency
#
#   Z-scores use healthy control references:
#     pupil_size:         mean = 4.10,  sd = 0.34
#     pupil_min:          mean = 2.70,  sd = 0.21
#     ch (%):             mean = 36.16, sd = 6.04
#     max_const_velocity: mean = 4.05,  sd = 0.90
#     latency:            mean = 0.24,  sd = 0.40
#
#   Performance (10-fold grouped CV):
#     R²   = 0.94
#     RMSE = 0.16 NPi units
# =============================================================================
setwd("L:/Auditdata/CONNECT-ME/DTU/FrederikWeinan_Thesis/R_attempt")

library(tidyverse)
library(scam)

# ── Load and prepare data ─────────────────────────────────────────────────────
dat <- read_csv("NPi_GAM_z_input.csv", show_col_types = FALSE) %>%
  mutate(
    record_id = factor(record_id),
    eye       = factor(eye)
  )

# ── Create patient-level folds ────────────────────────────────────────────────
set.seed(123)
patient_ids <- unique(dat$record_id)
fold_assign <- sample(rep(1:10, length.out = length(patient_ids)))

fold_map <- data.frame(
  record_id = patient_ids,
  fold      = fold_assign
)

dat <- dat %>%
  left_join(fold_map, by = "record_id")

# ── Helper: hinge (hockey stick) function ─────────────────────────────────────
hockey <- function(x, bp) pmax(x - bp, 0)

# ── Helper: apply formula, cap to [0, 5], round to 1 decimal ─────────────────
predict_npi <- function(data) {
  raw <- 5.771 +
    0.359  * data$z_pupil_size         - 0.081 * hockey(data$z_pupil_size, -4) +
    -0.283  * data$z_pupil_min          - 0.123 * hockey(data$z_pupil_min,  -4) +
    -0.235  * data$z_ch                 - 0.046 * data$z_ch^2                   +
    0.212  * data$z_max_const_velocity - 0.206 * hockey(data$z_max_const_velocity, -3) +
    -0.049  * data$z_latency
  
  # Hard cap to [0, 5] then round to 1 decimal
  round(pmin(pmax(raw, 0), 5), 1)
}

# ── 10-fold grouped cross-validation ─────────────────────────────────────────
all_preds <- list()

for (f in 1:10) {
  
  train <- dat %>% filter(fold != f)
  test  <- dat %>% filter(fold == f)
  
  # Fit SCAM on training data (no random effect — fixed effects only)
  scam_fold <- scam(
    npi ~ eye +
      s(z_pupil_size,          bs = "mpi") +
      s(z_pupil_min,           bs = "mpd") +
      s(z_ch,                  bs = "cr")  +
      s(z_max_const_velocity,  bs = "mpi") +
      s(z_dilat_velocity,      bs = "cr")  +
      s(z_latency,             bs = "cr"),
    data = train
  )
  
  # Get SCAM predictions on training data
  train$scam_pred <- as.numeric(predict(scam_fold, newdata = train))
  
  # Fit the interpretable formula to SCAM training predictions
  approx_fold <- lm(
    scam_pred ~
      z_pupil_size + hockey(z_pupil_size, -4) +
      z_pupil_min  + hockey(z_pupil_min,  -4) +
      poly(z_ch, 2, raw = TRUE)               +
      z_max_const_velocity + hockey(z_max_const_velocity, -3) +
      z_latency,
    data = train
  )
  
  # Predict on held-out patients, cap and round
  raw_pred      <- as.numeric(predict(approx_fold, newdata = test))
  test$pred     <- round(pmin(pmax(raw_pred, 0), 5), 1)
  test$fold     <- f
  
  all_preds[[f]] <- test
}

# ── Combine results ───────────────────────────────────────────────────────────
results <- bind_rows(all_preds)

# ── Performance metrics ───────────────────────────────────────────────────────
r2_cv   <- 1 - sum((results$npi - results$pred)^2) /
  sum((results$npi - mean(results$npi))^2)
rmse_cv <- sqrt(mean((results$npi - results$pred)^2))

cat("===========================================\n")
cat("   10-FOLD CV PERFORMANCE (held-out patients)\n")
cat("===========================================\n")
cat("R²:  ", round(r2_cv,   4), "\n")
cat("RMSE:", round(rmse_cv, 4), "NPi units\n\n")

# ── Apply formula directly to full dataset ────────────────────────────────────
dat$npi_predicted <- predict_npi(dat)

# ── Plots ─────────────────────────────────────────────────────────────────────
par(mfrow = c(1, 2))

# 1. CV predicted vs observed
plot(
  results$npi, results$pred,
  pch = 16, col = rgb(0, 0, 0, 0.3),
  xlab = "Observed NPi",
  ylab = "Predicted NPi (10-fold CV)",
  main = "Predicted vs Observed (CV)",
  xlim = c(0, 5), ylim = c(0, 5)
)
abline(0, 1, col = "red", lwd = 2)

# 2. Residuals
plot(
  results$npi, results$pred - results$npi,
  pch = 16, col = rgb(0, 0, 0, 0.3),
  xlab = "Observed NPi",
  ylab = "Residual (Predicted - Observed)",
  main = "Residual Pattern (CV)",
  ylim = c(-1.5, 1.5)
)
abline(h = 0, col = "red", lwd = 2)

par(mfrow = c(1, 1))

# ── Print formula ────────────── ( This is based on the results from above)
cat("===========================================\n")
cat("        NPi RECONSTRUCTION FORMULA\n")
cat("===========================================\n\n")
cat("NPi = 5.771\n")
cat("    + 0.359 * z_pupil_size\n")
cat("    - 0.081 * max(z_pupil_size + 4, 0)\n")
cat("    - 0.283 * z_pupil_min\n")
cat("    - 0.123 * max(z_pupil_min + 4, 0)\n")
cat("    - 0.235 * z_ch\n")
cat("    - 0.046 * z_ch^2\n")
cat("    + 0.212 * z_max_const_velocity\n")
cat("    - 0.206 * max(z_max_const_velocity + 3, 0)\n")
cat("    - 0.049 * z_latency\n\n")
cat("Predictions capped to [0, 5] and rounded to 1 decimal.\n\n")
cat("Z-score reference (healthy controls):\n")
cat("  pupil_size:         mean = 4.10,  sd = 0.34\n")
cat("  pupil_min:          mean = 2.70,  sd = 0.21\n")
cat("  ch (%):             mean = 36.16, sd = 6.04\n")
cat("  max_const_velocity: mean = 4.05,  sd = 0.90\n")
cat("  latency:            mean = 0.24,  sd = 0.40\n")
cat("===========================================\n")
