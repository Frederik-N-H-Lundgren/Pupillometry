# This document was the initial GAM exploration. The code in this document is not structured to be ran
# in one continuous run, but is instead comprised of different approached which assisted in later analysis



# ----------------------------
# Create patient-level folds
# ----------------------------

set.seed(123)

patient_ids <- unique(dat$record_id)
fold_assign <- sample(rep(1:10, length.out = length(patient_ids)))

fold_map <- data.frame(
  record_id = patient_ids,
  fold = fold_assign
)

dat <- dat %>%
  left_join(fold_map, by = "record_id")

# ----------------------------
# CV loop
# ----------------------------

all_preds <- list()

for (f in 1:10) {
  
  train <- dat %>% filter(fold != f)
  test  <- dat %>% filter(fold == f)
  
  gam_z <- gam(
    npi ~ eye +
      s(z_pupil_size, k = 6) +
      s(z_pupil_min, k = 6) +
      s(z_ch, k = 6) +
      s(z_const_velocity, k = 6) +
      s(z_max_const_velocity, k = 6) +
      s(z_dilat_velocity, k = 6) +
      s(z_latency, k = 6),
    data = train,
    method = "REML"
  )
  
  # IMPORTANT FIX: force numeric vector (not matrix)
  test$pred <- as.numeric(predict(gam_z, newdata = test))
  
  test$fold <- f
  
  all_preds[[f]] <- test
}

# ----------------------------
# Combine safely
# ----------------------------

results <- bind_rows(all_preds)

# ensure no list-columns remain
results <- results %>%
  mutate(across(where(is.list), ~ as.character(.)))

# final safety check
results$pred <- as.numeric(results$pred)


rmse <- sqrt(mean((results$npi - results$pred)^2))

r2 <- 1 - sum((results$npi - results$pred)^2) /
  sum((results$npi - mean(results$npi))^2)

rmse
r2

# This is we want to see the model in test result
summary(gam_z)

# A plot
plot(
  results$npi,
  results$pred,
  pch = 16,
  col = rgb(0, 0, 0, 0.3),
  xlab = "Observed NPi",
  ylab = "Predicted NPi (10-fold CV)"
)

abline(0, 1, col = "red", lwd = 2)




## Getting smooth curves now:
library(gratia)
library(tidyverse)

sm_pupil_size <- smooth_estimates(gam_z, select = "s(z_pupil_size)")
sm_pupil_min  <- smooth_estimates(gam_z, select = "s(z_pupil_min)")
sm_ch         <- smooth_estimates(gam_z, select = "s(z_ch)")
sm_mcv        <- smooth_estimates(gam_z, select = "s(z_max_const_velocity)")
sm_dv         <- smooth_estimates(gam_z, select = "s(z_dilat_velocity)")


# Dont need to run anymore
# Plot for each of the significant features
plot_smooth <- function(sm, xcol, title) {
  plot(
    sm[[xcol]],
    sm$.estimate,
    type = "l",
    xlab = xcol,
    ylab = "NPi contribution",
    main = title
  )
  abline(h = 0, lty = 2)
}
plot_smooth(sm_pupil_size, "z_pupil_size", "Pupil size effect")
plot_smooth(sm_pupil_min, "z_pupil_min", "Pupil min effect")
plot_smooth(sm_ch, "z_ch", "CH effect")
plot_smooth(sm_mcv, "z_max_const_velocity", "Max velocity effect")
plot_smooth(sm_dv, "z_dilat_velocity", "Dilation velocity")


# Trying to get a closer estimate of each of the function
# So i dont need the full complex functions
deriv_pupil_size <- derivatives(gam_z, select = "s(z_pupil_size)")

deriv_pupil_min <- derivatives(gam_z, select = "s(z_pupil_min)")
range(deriv_pupil_min$.derivative)

deriv_ch <- derivatives(gam_z, select = "s(z_ch)")
range(deriv_ch$.derivative)

deriv_mcv <- derivatives(gam_z, select = "s(z_max_const_velocity)")
range(deriv_mcv$.derivative)

deriv_dv <- derivatives(gam_z, select = "s(z_dilat_velocity)")
range(deriv_dv$.derivative)
sm_dv %>%
  summarise(
    left_mean = mean(.estimate[z_dilat_velocity < -1]),
    right_mean = mean(.estimate[z_dilat_velocity > 1])
  )


# Trying to do a GAM with restrictions to how the functions can look
library(scam)
library(tidyverse)
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
library(scam)
library(dplyr)

set.seed(123)

all_preds <- list()

for (f in 1:10) {
  
  # ----------------------------
  # split by patient fold
  # ----------------------------
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
  
  # ----------------------------
  # predict on held-out patients
  # ----------------------------
  test$pred <- as.numeric(predict(scam_fit, newdata = test))
  test$fold <- f
  
  all_preds[[f]] <- test
}

# ----------------------------
# combine predictions
# ----------------------------
results_scam <- bind_rows(all_preds)

# safety: ensure numeric
results_scam$pred <- as.numeric(results_scam$pred)

# ----------------------------
# metrics (true generalization)
# ----------------------------
rmse_scam <- sqrt(mean((results_scam$npi - results_scam$pred)^2))

r2_scam <- 1 - sum((results_scam$npi - results_scam$pred)^2) /
  sum((results_scam$npi - mean(results_scam$npi))^2)

rmse_scam  # 0.08981292
r2_scam  # 0.9809098

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


plot(scam_z, pages = 1)
#
#
#
#
#
#
#
#
#
# Trying to get a real function based on the smooth curves:
summary(scam_z)$s.table


# FOR MIN

library(mgcv)
library(gratia)
library(segmented)

# -----------------------------
# Extract smooth functions
# -----------------------------

get_smooth <- function(model, term_name) {
  sm <- smooth_estimates(model, smooth = term_name)
  data.frame(
    x = sm[[2]],
    y = sm$.estimate
  )
}

# -----------------------------
# Fit candidate approximations
# -----------------------------

fit_models <- function(df) {
  
  x <- df$x
  y <- df$y
  
  # 1. Linear
  lin <- lm(y ~ x)
  pred_lin <- predict(lin)
  rmse_lin <- sqrt(mean((y - pred_lin)^2))
  
  # 2. Quadratic
  quad <- lm(y ~ poly(x, 2))
  pred_quad <- predict(quad)
  rmse_quad <- sqrt(mean((y - pred_quad)^2))
  
  # 3. Hinge (1 breakpoint)
  base <- lm(y ~ x)
  seg <- try(segmented(base, seg.Z = ~x, npsi = 1), silent = TRUE)
  
  if (!inherits(seg, "try-error")) {
    pred_hinge <- predict(seg)
    rmse_hinge <- sqrt(mean((y - pred_hinge)^2))
  } else {
    rmse_hinge <- NA
  }
  
  list(
    rmse_linear = rmse_lin,
    rmse_quad = rmse_quad,
    rmse_hinge = rmse_hinge
  )
}

# -----------------------------
# Run on all SCAM smooths
# -----------------------------

features <- c(
  "s(z_pupil_size)",
  "s(z_pupil_min)",
  "s(z_ch)",
  "s(z_max_const_velocity)",
  "s(z_dilat_velocity)",
  "s(z_latency)"
)

results <- list()

for (f in features) {
  
  df_smooth <- get_smooth(scam_z, f)
  fits <- fit_models(df_smooth)
  
  results[[f]] <- fits
}

results
































# Rest from here is reconstruction

# Dont need to run anymore
plot(
  sm_pupil_size$z_pupil_size,
  sm_pupil_size$.estimate,
  type = "l",
  lwd = 2,
  xlab = "Pupil size (z-score)",
  ylab = "Contribution to NPi"
)

abline(h = 0, lty = 2)
abline(v = 0, lty = 2, col = "grey")

# Dont need to run anymore
x0_est <- sm_pupil_size$z_pupil_size[
  which.min(abs(sm_pupil_size$.estimate))
]
x0_est


# Need to run from here on

# for the pupil size:
x0 <- -2.18
k_pupil <- 1

Q_pupil_size <- function(x) {
  1 / (1 + exp(-k_pupil * (x - x0)))
}


# pupil min:
x0_min <- 4  # your observed inflection region
k_min <- 1

Q_pupil_min <- function(x) {
  1 / (1 + exp(k_min * (x - x0_min)))
}

# ch:
x0_ch <- -1

Q_ch <- function(x) {
  -(x - x0_ch)^2
}

# Mcv
# Because there is not a clear function that fits this, we keep the smooth curve
Q_mcv <- function(x) {
  approx(
    x = sm_mcv$z_max_const_velocity,
    y = sm_mcv$.estimate,
    xout = x,
    rule = 2
  )$y
}

# Dv
# SAme for this as previous
Q_dv <- function(x) {
  approx(
    x = sm_dv$z_dilat_velocity,
    y = sm_dv$.estimate,
    xout = x,
    rule = 2
  )$y
}



# Now time to put it all together

dat$Q_pupil_size <- Q_pupil_size(dat$z_pupil_size)
dat$Q_pupil_min  <- Q_pupil_min(dat$z_pupil_min)
dat$Q_ch         <- Q_ch(dat$z_ch)
dat$Q_mcv        <- Q_mcv(dat$z_max_const_velocity)
dat$Q_dv         <- Q_dv(dat$z_dilat_velocity)


# adding weights to each Q-function
model <- lm(npi ~ Q_pupil_size + Q_pupil_min + Q_ch + Q_mcv + Q_dv, data = dat)




dat$NPi_hat <- predict(model)


rmse <- sqrt(mean((dat$npi - dat$NPi_hat)^2))
rmse

r2 <- 1 - sum((dat$npi - dat$NPi_hat)^2) /
  sum((dat$npi - mean(dat$npi))^2)
r2

plot(dat$npi, dat$NPi_hat,
     xlab = "Observed NPi",
     ylab = "Predicted NPi",
     pch = 16)
abline(0, 1, col = "red")
