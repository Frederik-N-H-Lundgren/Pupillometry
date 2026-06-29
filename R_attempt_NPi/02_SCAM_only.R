# This document came from explorative analysis in 01_fit_GAM_zscores and is the initial SCAM attempt.
# In this code the features are fitted while being lilited to a certain degree, to see if there was a 
# biological plausible fit for each feature.


# Clean SCAM only so it is more readable

setwd("L:/Auditdata/CONNECT-ME/DTU/FrederikWeinan_Thesis/R_attempt")

# Data and library loading

library(tidyverse)
library(mgcv)
library(gratia)
library(scam)
library(dplyr)

dat <- read_csv("NPi_GAM_z_input.csv", show_col_types = FALSE) %>%
  mutate(
    record_id = factor(record_id),
    eye = factor(eye)
  )


set.seed(123)

patient_ids <- unique(dat$record_id)
fold_assign <- sample(rep(1:10, length.out = length(patient_ids)))

fold_map <- data.frame(
  record_id = patient_ids,
  fold = fold_assign
)

dat <- dat %>%
  left_join(fold_map, by = "record_id")


# Doing SCAM on all data
scam_z <- scam(
  npi ~ eye +
    s(z_pupil_size, bs = "mpi") +
    s(z_pupil_min, bs = "mpd") +
    s(z_ch, bs = "cr") +
    s(z_max_const_velocity, bs = "mpi") +
    s(z_dilat_velocity, bs = "cr") +
    s(z_latency, bs = "cr") +
    s(record_id, bs = "re"),
  data = dat
)
summary(scam_z)


## Now doing the SCAM with cross val:



all_preds <- list()

for (f in 1:10) {
  
  # split by patient fold
  train <- dat %>% filter(fold != f)
  test  <- dat %>% filter(fold == f)
  
  # ----------------------------
  # fit SCAM on training data only
  # ----------------------------
  scam_fit <- scam(
    npi ~ eye +
      s(z_pupil_size, bs = "mpi") +
      s(z_pupil_min, bs = "mpd") +
      s(z_ch, bs = "cr") +
      s(z_max_const_velocity, bs = "mpi") +
      s(z_dilat_velocity, bs = "cr") +
      s(z_latency, bs = "cr") +
      s(record_id, bs = "re"),
    data = train
  )
  
  # predict on held-out patients
  test$pred <- as.numeric(predict(scam_fit, newdata = test))
  test$fold <- f
  
  all_preds[[f]] <- test
}


# combine predictions

results_scam <- bind_rows(all_preds)

# safety: ensure numeric
results_scam$pred <- as.numeric(results_scam$pred)

# Performance
rmse_scam <- sqrt(mean((results_scam$npi - results_scam$pred)^2))

r2_scam <- 1 - sum((results_scam$npi - results_scam$pred)^2) /
  sum((results_scam$npi - mean(results_scam$npi))^2)

rmse_scam  # 0.08981292 (from first few runs)
r2_scam  # 0.9809098


# plots
plot(
  results_scam$npi,
  results_scam$pred,
  pch = 16,
  col = rgb(0,0,0,0.3),
  xlab = "Observed NPi",
  ylab = "Predicted NPi (SCAM CV)",
  main = "SCAM: Predicted vs Observed"
)

abline(0, 1, col = "red", lwd = 2)

# Another plot:
plot(
  results_scam$npi,
  results_scam$pred - results_scam$npi,
  pch = 16,
  col = rgb(0,0,0,0.3),
  xlab = "Observed NPi",
  ylab = "Residual (Pred - Obs)",
  main = "Residual pattern"
)
abline(h = 0, col = "red", lwd = 2)


plot(scam_fit, pages = 1)

# Now trying to get the functions out
# Extract smooth curves from the full SCAM model 
library(gratia)

sm_pupil_size <- smooth_estimates(scam_z, select = "s(z_pupil_size)") %>% add_confint()
sm_pupil_min  <- smooth_estimates(scam_z, select = "s(z_pupil_min)")  %>% add_confint()
sm_ch         <- smooth_estimates(scam_z, select = "s(z_ch)")         %>% add_confint()
sm_mcv        <- smooth_estimates(scam_z, select = "s(z_max_const_velocity)") %>% add_confint()
sm_dv         <- smooth_estimates(scam_z, select = "s(z_dilat_velocity)")     %>% add_confint()
sm_lat        <- smooth_estimates(scam_z, select = "s(z_latency)")            %>% add_confint()

Plot each smooth so we can see its shape 
par(mfrow = c(2, 3))

plot(sm_pupil_size$z_pupil_size, sm_pupil_size$.estimate, type = "l",
     main = "z_pupil_size", xlab = "z", ylab = "smooth estimate")
abline(h = 0, lty = 2)

plot(sm_pupil_min$z_pupil_min, sm_pupil_min$.estimate, type = "l",
     main = "z_pupil_min", xlab = "z", ylab = "smooth estimate")
abline(h = 0, lty = 2)

plot(sm_ch$z_ch, sm_ch$.estimate, type = "l",
     main = "z_ch", xlab = "z", ylab = "smooth estimate")
abline(h = 0, lty = 2)

plot(sm_mcv$z_max_const_velocity, sm_mcv$.estimate, type = "l",
     main = "z_max_const_velocity", xlab = "z", ylab = "smooth estimate")
abline(h = 0, lty = 2)

plot(sm_dv$z_dilat_velocity, sm_dv$.estimate, type = "l",
     main = "z_dilat_velocity", xlab = "z", ylab = "smooth estimate")
abline(h = 0, lty = 2)

plot(sm_lat$z_latency, sm_lat$.estimate, type = "l",
     main = "z_latency", xlab = "z", ylab = "smooth estimate")
abline(h = 0, lty = 2)

par(mfrow = c(1, 1))

# Fit polynomial approximations to each smooth
# Starting with degree 2 for all — we will adjust based on the plots

approx_pupil_size <- lm(.estimate ~ poly(z_pupil_size, 2, raw = TRUE), data = sm_pupil_size)
approx_pupil_min  <- lm(.estimate ~ poly(z_pupil_min,  2, raw = TRUE), data = sm_pupil_min)
approx_ch         <- lm(.estimate ~ poly(z_ch,         2, raw = TRUE), data = sm_ch)
approx_mcv        <- lm(.estimate ~ poly(z_max_const_velocity, 1, raw = TRUE), data = sm_mcv)
approx_dv         <- lm(.estimate ~ poly(z_dilat_velocity,     1, raw = TRUE), data = sm_dv)
approx_lat        <- lm(.estimate ~ poly(z_latency,            1, raw = TRUE), data = sm_lat)

# Check fit quality for each smooth approximation
cat("--- pupil_size ---\n"); print(summary(approx_pupil_size)$r.squared)
cat("--- pupil_min ---\n");  print(summary(approx_pupil_min)$r.squared)
cat("--- ch ---\n");         print(summary(approx_ch)$r.squared)
cat("--- mcv ---\n");        print(summary(approx_mcv)$r.squared)
cat("--- dv ---\n");         print(summary(approx_dv)$r.squared)
cat("--- lat ---\n");        print(summary(approx_lat)$r.squared)

# ── STEP 4: Predict each smooth's contribution for every row in dat ───────────
contrib_pupil_size <- predict(approx_pupil_size,
                              newdata = data.frame(z_pupil_size = dat$z_pupil_size))
contrib_pupil_min  <- predict(approx_pupil_min,
                              newdata = data.frame(z_pupil_min = dat$z_pupil_min))
contrib_ch         <- predict(approx_ch,
                              newdata = data.frame(z_ch = dat$z_ch))
contrib_mcv        <- predict(approx_mcv,
                              newdata = data.frame(z_max_const_velocity = dat$z_max_const_velocity))
contrib_dv         <- predict(approx_dv,
                              newdata = data.frame(z_dilat_velocity = dat$z_dilat_velocity))
contrib_lat        <- predict(approx_lat,
                              newdata = data.frame(z_latency = dat$z_latency))

# ── STEP 5: Combine into final prediction ─────────────────────────────────────
intercept <- coef(scam_z)[1]

NPi_pred <- intercept +
  contrib_pupil_size +
  contrib_pupil_min  +
  contrib_ch         +
  contrib_mcv        +
  contrib_dv         +
  contrib_lat

# ── STEP 6: Check overall R² of the approximated formula ─────────────────────
r2_approx <- 1 - sum((dat$npi - NPi_pred)^2) / sum((dat$npi - mean(dat$npi))^2)
rmse_approx <- sqrt(mean((dat$npi - NPi_pred)^2))

cat("\nApproximated formula R²:  ", round(r2_approx, 4), "\n")
cat("Approximated formula RMSE:", round(rmse_approx, 4), "\n")

# ── STEP 7: Print the formula coefficients ────────────────────────────────────
cat("\n--- FORMULA COEFFICIENTS ---\n")
cat("Intercept:     ", round(intercept, 4), "\n\n")
cat("z_pupil_size:  ", round(coef(approx_pupil_size), 6), "\n")
cat("z_pupil_min:   ", round(coef(approx_pupil_min),  6), "\n")
cat("z_ch:          ", round(coef(approx_ch),          6), "\n")
cat("z_mcv:         ", round(coef(approx_mcv),         6), "\n")
cat("z_dv:          ", round(coef(approx_dv),           6), "\n")
cat("z_latency:     ", round(coef(approx_lat),          6), "\n")

# ── STEP 8: Plot approximated vs SCAM predictions ────────────────────────────
plot(
  dat$npi, NPi_pred,
  pch = 16, col = rgb(0, 0, 0, 0.3),
  xlab = "Observed NPi",
  ylab = "Approximated Formula Prediction",
  main = "Polynomial Approximation vs Observed"
)
abline(0, 1, col = "red", lwd = 2)






#
#
#
#
#
#
#
#
#


# New attempt with hinge because this might be better:

hockey <- function(x, bp) pmax(x - bp, 0)

# z_pupil_size — hockey stick at -4
approx_pupil_size <- lm(.estimate ~ z_pupil_size + hockey(z_pupil_size, -4),
                        data = sm_pupil_size)
cat("pupil_size R²:", summary(approx_pupil_size)$r.squared, "\n")

# z_pupil_min — hockey stick at -4
approx_pupil_min <- lm(.estimate ~ z_pupil_min + hockey(z_pupil_min, -4),
                       data = sm_pupil_min)
cat("pupil_min R²:", summary(approx_pupil_min)$r.squared, "\n")

# z_ch — quadratic fits well already (R² = 0.987)
approx_ch <- lm(.estimate ~ poly(z_ch, 2, raw = TRUE), data = sm_ch)
cat("ch R²:", summary(approx_ch)$r.squared, "\n")

# z_max_const_velocity — hockey stick at -3
approx_mcv <- lm(.estimate ~ z_max_const_velocity + hockey(z_max_const_velocity, -3),
                 data = sm_mcv)
cat("mcv R²:", summary(approx_mcv)$r.squared, "\n")

# z_dilat_velocity — drop (y-range ~0.04, negligible)
# z_latency — drop (y-range ~0.04, negligible)


contrib_pupil_size <- predict(approx_pupil_size,
                              newdata = data.frame(z_pupil_size = dat$z_pupil_size))

contrib_pupil_min <- predict(approx_pupil_min,
                             newdata = data.frame(z_pupil_min = dat$z_pupil_min))

contrib_ch <- predict(approx_ch,
                      newdata = data.frame(z_ch = dat$z_ch))

contrib_mcv <- predict(approx_mcv,
                       newdata = data.frame(z_max_const_velocity = dat$z_max_const_velocity))


intercept <- coef(scam_z)[1]

NPi_pred <- intercept +
  contrib_pupil_size +
  contrib_pupil_min  +
  contrib_ch         +
  contrib_mcv


r2_approx   <- 1 - sum((dat$npi - NPi_pred)^2) / sum((dat$npi - mean(dat$npi))^2)
rmse_approx <- sqrt(mean((dat$npi - NPi_pred)^2))

cat("\nApproximated formula R²:  ", round(r2_approx, 4), "\n")
cat("Approximated formula RMSE:", round(rmse_approx, 4), "\n")


cat("\n--- FORMULA COEFFICIENTS ---\n")
cat("Intercept:        ", round(intercept, 4), "\n\n")

cat("z_pupil_size:     ", round(coef(approx_pupil_size), 5), "\n")
cat("  (breakpoint at z = -4)\n\n")

cat("z_pupil_min:      ", round(coef(approx_pupil_min), 5), "\n")
cat("  (breakpoint at z = -4)\n\n")

cat("z_ch (quadratic): ", round(coef(approx_ch), 5), "\n\n")

cat("z_mcv:            ", round(coef(approx_mcv), 5), "\n")
cat("  (breakpoint at z = -3)\n\n")


plot(
  dat$npi, NPi_pred,
  pch = 16, col = rgb(0, 0, 0, 0.3),
  xlab = "Observed NPi",
  ylab = "Approximated Formula Prediction",
  main = "Polynomial/Piecewise Approximation vs Observed",
  xlim = c(0, 5), ylim = c(0, 5)
)
abline(0, 1, col = "red", lwd = 2)
















# Trying to fix some intercept thing with the mean of each feature
# Check how much variance the random effect explains
scam_z_no_re <- scam(
  npi ~ eye +
    s(z_pupil_size, bs = "mpi") +
    s(z_pupil_min, bs = "mpd") +
    s(z_ch, bs = "cr") +
    s(z_max_const_velocity, bs = "mpi") +
    s(z_dilat_velocity, bs = "cr") +
    s(z_latency, bs = "cr"),
  data = dat
)

pred_no_re <- predict(scam_z_no_re)
r2_no_re <- 1 - sum((dat$npi - pred_no_re)^2) / sum((dat$npi - mean(dat$npi))^2)
cat("R² without random effect:", round(r2_no_re, 4), "\n")


pred_scam_fixed <- as.numeric(predict(scam_z_no_re, newdata = dat))
hockey <- function(x, bp) pmax(x - bp, 0)

approx_formula <- lm(
  pred_scam_fixed ~
    z_pupil_size + hockey(z_pupil_size, -4) +
    z_pupil_min  + hockey(z_pupil_min, -4)  +
    poly(z_ch, 2, raw = TRUE)               +
    z_max_const_velocity + hockey(z_max_const_velocity, -3) +
    z_dilat_velocity                        +
    z_latency,
  data = dat
)

cat("Approximation R² (vs SCAM predictions):", 
    round(summary(approx_formula)$r.squared, 4), "\n")

# ── Now check against actual NPi ──────────────────────────────────────────────
NPi_pred <- predict(approx_formula, newdata = dat)

r2_vs_npi <- 1 - sum((dat$npi - NPi_pred)^2) / sum((dat$npi - mean(dat$npi))^2)
rmse_vs_npi <- sqrt(mean((dat$npi - NPi_pred)^2))

cat("R² vs observed NPi:  ", round(r2_vs_npi, 4), "\n")
cat("RMSE vs observed NPi:", round(rmse_vs_npi, 4), "\n")

# ── Print the actual formula ───────────────────────────────────────────────────
cat("\n--- FORMULA ---\n")
print(round(coef(approx_formula), 5))

# ── Plot ──────────────────────────────────────────────────────────────────────
plot(
  dat$npi, NPi_pred,
  pch = 16, col = rgb(0,0,0,0.3),
  xlab = "Observed NPi",
  ylab = "Approximated Prediction",
  main = "Direct SCAM Approximation",
  xlim = c(0,5), ylim = c(0,5)
)
abline(0, 1, col = "red", lwd = 2)


all_preds_approx <- list()

for (f in 1:10) {
  train <- dat %>% filter(fold != f)
  test  <- dat %>% filter(fold == f)
  
  # Fit no-RE SCAM on training data
  scam_fold <- scam(
    npi ~ eye +
      s(z_pupil_size, bs = "mpi") +
      s(z_pupil_min, bs = "mpd") +
      s(z_ch, bs = "cr") +
      s(z_max_const_velocity, bs = "mpi") +
      s(z_dilat_velocity, bs = "cr") +
      s(z_latency, bs = "cr"),
    data = train
  )
  
  # Get SCAM predictions on training data
  train$scam_pred <- as.numeric(predict(scam_fold, newdata = train))
  
  # Fit approximation formula to training SCAM predictions
  approx_fold <- lm(
    scam_pred ~
      z_pupil_size + hockey(z_pupil_size, -4) +
      z_pupil_min  + hockey(z_pupil_min, -4)  +
      poly(z_ch, 2, raw = TRUE)               +
      z_max_const_velocity + hockey(z_max_const_velocity, -3) +
      z_dilat_velocity +
      z_latency,
    data = train
  )
  
  # Predict on held-out patients using the simple formula
  test$pred <- as.numeric(predict(approx_fold, newdata = test))
  test$fold <- f
  all_preds_approx[[f]] <- test
}

results_approx <- bind_rows(all_preds_approx)

r2_cv <- 1 - sum((results_approx$npi - results_approx$pred)^2) /
  sum((results_approx$npi - mean(results_approx$npi))^2)
rmse_cv <- sqrt(mean((results_approx$npi - results_approx$pred)^2))

cat("CV R²:  ", round(r2_cv, 4), "\n")
cat("CV RMSE:", round(rmse_cv, 4), "\n")



# TEST 2 MEHTODS TO SEE IF WE CAN GET CLEANER FORMULAR
approx_formula_clean <- lm(
  pred_scam_fixed ~
    z_pupil_size + hockey(z_pupil_size, -4) +
    z_pupil_min  + hockey(z_pupil_min, -4)  +
    poly(z_ch, 2, raw = TRUE)               +
    z_max_const_velocity + hockey(z_max_const_velocity, -3) +
    z_latency,
  data = dat
)

NPi_pred_clean <- predict(approx_formula_clean, newdata = dat)
r2_clean <- 1 - sum((dat$npi - NPi_pred_clean)^2) / 
  sum((dat$npi - mean(dat$npi))^2)
cat("R² without z_dilat_velocity:", round(r2_clean, 4), "\n")

approx_formula_linear <- lm(
  pred_scam_fixed ~
    z_pupil_size +
    z_pupil_min  +
    poly(z_ch, 2, raw = TRUE) +
    z_max_const_velocity +
    z_latency,
  data = dat
)

NPi_pred_linear <- predict(approx_formula_linear, newdata = dat)
r2_linear <- 1 - sum((dat$npi - NPi_pred_linear)^2) / 
  sum((dat$npi - mean(dat$npi))^2)
cat("R² fully linear (no hockey):", round(r2_linear, 4), "\n")



# FINAL RUN!
all_preds_final <- list()

for (f in 1:10) {
  train <- dat %>% filter(fold != f)
  test  <- dat %>% filter(fold == f)
  
  # Fit no-RE SCAM on training data
  scam_fold <- scam(
    npi ~ eye +
      s(z_pupil_size, bs = "mpi") +
      s(z_pupil_min, bs = "mpd") +
      s(z_ch, bs = "cr") +
      s(z_max_const_velocity, bs = "mpi") +
      s(z_dilat_velocity, bs = "cr") +
      s(z_latency, bs = "cr"),
    data = train
  )
  
  # Get SCAM predictions on training data
  train$scam_pred <- as.numeric(predict(scam_fold, newdata = train))
  
  # Fit clean formula to training SCAM predictions
  approx_fold <- lm(
    scam_pred ~
      z_pupil_size + hockey(z_pupil_size, -4) +
      z_pupil_min  + hockey(z_pupil_min, -4)  +
      poly(z_ch, 2, raw = TRUE)               +
      z_max_const_velocity + hockey(z_max_const_velocity, -3) +
      z_latency,
    data = train
  )
  
  # Predict on held-out patients
  test$pred <- as.numeric(predict(approx_fold, newdata = test))
  test$fold <- f
  all_preds_final[[f]] <- test
}

results_final <- bind_rows(all_preds_final)

r2_final <- 1 - sum((results_final$npi - results_final$pred)^2) /
  sum((results_final$npi - mean(results_final$npi))^2)
rmse_final <- sqrt(mean((results_final$npi - results_final$pred)^2))

cat("Final CV R²:  ", round(r2_final, 4), "\n")
cat("Final CV RMSE:", round(rmse_final, 4), "\n")

# ── Refit on ALL data for the final formula coefficients ──────────────────────
pred_scam_fixed <- as.numeric(predict(scam_z_no_re, newdata = dat))

final_formula <- lm(
  pred_scam_fixed ~
    z_pupil_size + hockey(z_pupil_size, -4) +
    z_pupil_min  + hockey(z_pupil_min, -4)  +
    poly(z_ch, 2, raw = TRUE)               +
    z_max_const_velocity + hockey(z_max_const_velocity, -3) +
    z_latency,
  data = dat
)

cat("\n--- FINAL FORMULA COEFFICIENTS ---\n")
print(round(coef(final_formula), 5))

# ── Final plot ────────────────────────────────────────────────────────────────
plot(
  results_final$npi, results_final$pred,
  pch = 16, col = rgb(0,0,0,0.3),
  xlab = "Observed NPi",
  ylab = "Predicted NPi (CV)",
  main = "Final Formula: 10-fold CV",
  xlim = c(0,5), ylim = c(0,5)
)
abline(0, 1, col = "red", lwd = 2)


# SIGMOID ATTEMPT
# Sigmoid-based formula
sigmoid <- function(x) 5 / (1 + exp(-x))

# Build linear predictor first
eta <- with(dat,
            -1.5 +                        # rough starting intercept
              0.36  * z_pupil_size +
              -0.28 * z_pupil_min  +
              -0.24 * z_ch         +
              -0.046* z_ch^2       +
              0.21  * z_max_const_velocity +
              -0.05 * z_latency
)

# Then wrap in sigmoid scaled to 0-5
NPi_sigmoid <- sigmoid(eta)

r2_sigmoid <- 1 - sum((dat$npi - NPi_sigmoid)^2) /
  sum((dat$npi - mean(dat$npi))^2)
cat("Sigmoid formula R²:", round(r2_sigmoid, 4), "\n")



# ── Print the formula cleanly ─────────────────────────────────────────────────
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
cat("Where z = (x - mean) / sd using healthy control references:\n")
cat("  pupil_size:          mean = 4.10,  sd = 0.34\n")
cat("  pupil_min:           mean = 2.70,  sd = 0.21\n")
cat("  ch (%):              mean = 36.16, sd = 6.04\n")
cat("  max_const_velocity:  mean = 4.05,  sd = 0.90\n")
cat("  latency:             mean = 0.24,  sd = 0.40\n")
cat("===========================================\n\n")



## PLOT ALL THE FUNCTIONS:
# ── Plot each term's contribution ────────────────────────────────────────────
hockey <- function(x, bp) pmax(x - bp, 0)

par(mfrow = c(2, 3))

# z_pupil_size
z <- seq(-8, 4, length.out = 200)
contrib <- 0.359 * z - 0.081 * hockey(z, -4)
plot(z, contrib, type = "l", lwd = 2, col = "steelblue",
     main = "Pupil Size contribution",
     xlab = "z_pupil_size", ylab = "NPi contribution")
abline(v = -4, lty = 2, col = "gray")
abline(h = 0,  lty = 2, col = "gray")
text(-4, min(contrib) + 0.3, "hinge at -4", col = "gray40", cex = 0.8)

# z_pupil_min
z <- seq(-8, 4, length.out = 200)
contrib <- -0.283 * z - 0.123 * hockey(z, -4)
plot(z, contrib, type = "l", lwd = 2, col = "darkorange",
     main = "Pupil Min contribution",
     xlab = "z_pupil_min", ylab = "NPi contribution")
abline(v = -4, lty = 2, col = "gray")
abline(h = 0,  lty = 2, col = "gray")
text(-4, max(contrib) - 0.3, "hinge at -4", col = "gray40", cex = 0.8)

# z_ch
z <- seq(-6, 3, length.out = 200)
contrib <- -0.235 * z - 0.046 * z^2
plot(z, contrib, type = "l", lwd = 2, col = "darkgreen",
     main = "CH (%) contribution",
     xlab = "z_ch", ylab = "NPi contribution")
abline(h = 0, lty = 2, col = "gray")

# z_max_const_velocity
z <- seq(-6, 4, length.out = 200)
contrib <- 0.212 * z - 0.206 * hockey(z, -3)
plot(z, contrib, type = "l", lwd = 2, col = "firebrick",
     main = "Max Const. Velocity contribution",
     xlab = "z_max_const_velocity", ylab = "NPi contribution")
abline(v = -3, lty = 2, col = "gray")
abline(h = 0,  lty = 2, col = "gray")
text(-3, min(contrib) + 0.2, "hinge at -3", col = "gray40", cex = 0.8)

# z_latency
z <- seq(-1, 4, length.out = 200)
contrib <- -0.049 * z
plot(z, contrib, type = "l", lwd = 2, col = "purple",
     main = "Latency contribution",
     xlab = "z_latency", ylab = "NPi contribution")
abline(h = 0, lty = 2, col = "gray")

# Overall predicted vs observed
NPi_pred_final <- predict(final_formula, newdata = dat)
NPi_pred_final <- pmin(pmax(NPi_pred_final, 0), 5)

plot(
  dat$npi, NPi_pred_final,
  pch = 16, col = rgb(0,0,0,0.3),
  xlab = "Observed NPi",
  ylab = "Predicted NPi",
  main = "Formula: Predicted vs Observed",
  xlim = c(0,5), ylim = c(0,5)
)
abline(0, 1, col = "red", lwd = 2)

par(mfrow = c(1,1))







