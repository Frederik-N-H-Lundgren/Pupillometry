# This code is an part of the originally provided R document (from article), 
# which has isolated the calculation of T50




# =====================================================
# CLEAN ENVIRONMENT AND LOAD PACKAGES
# =====================================================
rm(list = ls())

packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "ggplot2"
)

install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)


# =====================================================
# 1. READ DATA (your existing function)
# =====================================================
read_all <- function(x) {
  sheet_names <- readxl::excel_sheets(x)
  all_days <- list()
  
  for (i in seq_along(sheet_names)) {
    tmp <- readxl::read_xlsx(x, sheet = sheet_names[i], col_names = FALSE)
    
    if (i == 1) {
      gcs_row <- 3; four_row <- 4; surv_row <- 8; id_row <- 9; start_row <- 10
    } else {
      gcs_row <- 2; four_row <- 3; surv_row <- NA; id_row <- 6; start_row <- 7
    }
    
    pt_ids <- as.character(unlist(tmp[id_row, -1]))
    GCS_vals <- suppressWarnings(as.numeric(unlist(tmp[gcs_row, -1])))
    FOUR_vals <- suppressWarnings(as.numeric(unlist(tmp[four_row, -1])))
    
    meta <- data.frame(pt_id = pt_ids, GCS = GCS_vals, FOUR = FOUR_vals, day = i)
    
    if (i == 1 && !is.na(surv_row)) {
      surv_vals <- as.character(unlist(tmp[surv_row, -1]))
      survival_df <- data.frame(
        pt_id = pt_ids,
        `90-day survival` = ifelse(
          surv_vals == "Y", "survived",
          ifelse(surv_vals == "N", "dead", "Unknown")
        )
      )
    }
    
    pup <- tmp[start_row:nrow(tmp), ]
    colnames(pup)[1] <- "time"
    colnames(pup)[-1] <- pt_ids
    
    pup_long <- pup %>%
      pivot_longer(-time, names_to = "pt_id", values_to = "size") %>%
      mutate(day = i)
    
    df_day <- left_join(pup_long, meta, by = c("pt_id", "day"))
    all_days[[i]] <- df_day
  }
  
  df <- bind_rows(all_days)
  if (exists("survival_df")) df <- left_join(df, survival_df, by = "pt_id")
  else df$`90-day survival` <- "Unknown"
  
  df$time <- as.numeric(df$time)
  df$size <- as.numeric(df$size)
  return(df)
}

# =====================================================
# 2. LOAD DATA
# =====================================================

setwd("L:/Auditdata/CONNECT-ME/DTU/FrederikWeinan_Thesis/Pupilometri") # Folder with the data

dfr <- read_all("Right_manually_cleaned_artefacts.xlsx")
dfr$lateral <- "right"
dfl <- read_all("Left_manually_cleaned_artefacts.xlsx")
dfl$lateral <- "left"
df <- rbind(dfr, dfl)

names(df)[names(df) %in% c("X90.day.survival", "X.90.day.survival", "X_90_day_survival")] <- "90-day survival"
df$day <- as.numeric(df$day)
df <- df[df$day %in% 1:12, ]

# =====================================================
# 3. COMPUTE T50 METRICS (robust version)
# =====================================================

compute_T50 <- function(time, size,
                        t_light_on  = 3,
                        t_light_off = 6,
                        t_late_light_off = 8) {
  
  # Remove NA pairs
  df <- tibble(time = time, size = size) %>%
    filter(!is.na(time), !is.na(size))
  
  if (nrow(df) < 5)
    return(tibble(T50_constr = NA_real_,
                  T50_dilat  = NA_real_))
  
  # =========================
  # 1. T50 CONSTRICTION (PLR)
  # =========================
  df_constr <- df %>%
    filter(time >= t_light_on, time <= t_light_off)
  
  T50_constr <- NA_real_
  
  if (nrow(df_constr) >= 3) {
    max_size <- max(df_constr$size, na.rm = TRUE)
    min_size <- min(df_constr$size, na.rm = TRUE)
    amp <- max_size - min_size
    
    if (amp > 0) {
      constr_thresh <- max_size - 0.5 * amp
      
      # size → time interpolation
      T50_constr <- tryCatch(
        approx(df_constr$size,
               df_constr$time,
               xout = constr_thresh,
               rule = 1)$y,
        error = function(e) NA_real_
      )
    }
  }
  
  # =========================
  # 2. T50 DILATION (LOR)
  # =========================
  df_dilat <- df %>%
    filter(time >= t_late_light_off)
  
  T50_dilat <- NA_real_
  
  if (nrow(df_dilat) >= 3) {
    min_size <- min(df_dilat$size, na.rm = TRUE)
    max_size <- max(df_dilat$size, na.rm = TRUE)
    amp <- max_size - min_size
    
    if (amp > 0) {
      dilat_thresh <- min_size + 0.5 * amp
      
      # reverse to ensure monotonic interpolation
      T50_dilat <- tryCatch(
        approx(rev(df_dilat$size),
               rev(df_dilat$time),
               xout = dilat_thresh,
               rule = 1)$y,
        error = function(e) NA_real_
      )
    }
  }
  
  tibble(
    T50_constr = T50_constr,
    T50_dilat  = T50_dilat
  )
}

# =====================================================
# 4. ##Metadata##
# =====================================================

df_T50 <- df %>%
  group_by(pt_id, day, lateral) %>%
  summarise(
    GCS = first(GCS),
    FOUR = first(FOUR),
    survival = first(`90-day survival`),
    data = list(tibble(time = time, size = size)),
    .groups = "drop"
  ) %>%
  mutate(metrics = map(data, ~ compute_T50(.x$time, .x$size))) %>%
  unnest(metrics)

summary(df_T50$T50_constr)
summary(df_T50$T50_dilat)