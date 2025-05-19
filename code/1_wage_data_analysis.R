###############################################################################
### Author: Alec McClean
### Purpose: Clean and analyze wage panel data
###############################################################################

set.seed(20250513)
options(scipen = 999)
source("0_functions.R")

## --------------------------------------------------
## 0.  Packages & raw data
## --------------------------------------------------

data(wagepan)                              # 545 workers × 8 years (1980-87)
wagepan <- wagepan %>% arrange(nr, year)   # keep rows in tidy panel order
wagepan <- wagepan %>% filter(year <= 1983) # Just use 4 years (1980-'83) to simplify analysis


## Add a numeric time index t = 1 … 3 (1980 --> 1, ..., 1983 --> 4)
wagepan <- wagepan %>%
  mutate(t = year - min(year) + 1)

## --------------------------------------------------
## 1.  Core objects for your pipeline
## --------------------------------------------------

# Subject IDs (length n = 545)
id <- wagepan %>% distinct(nr) %>% arrange(nr) %>% pull(nr)

# Baseline covariates: education, race, ethnicity
baseline_X <- wagepan %>%
  filter(year == 1980) %>%
  arrange(nr) %>%
  select(educ, black, hisp) %>%                        
  as.data.frame()

# Outcome — log wage in 1987
Y <- wagepan %>%
  filter(year == 1983) %>%
  arrange(nr) %>%
  pull(lwage)

# Time-varying covariates: hours worked & experience every year
time_X <- wagepan %>%
  select(nr, t, agric, bus, construc, exper, fin, poorhlth, hours, manuf,
         married, min, nrthcen, nrtheast, occ1:occ9, per, pro, pub, rur, 
         south, tra, trad) %>%
  pivot_wider(
    id_cols    = nr,  
    names_from = t,   
    values_from = -c(nr, t),
    names_glue = "{.value}_t{sprintf('%02d', t)}"
  ) %>%
  arrange(nr) 

# Add wage variables
time_X %<>% left_join(
  wagepan %>% 
  filter(year < 1983) %>%
  select(nr, t, lwage) %>% 
  mutate(t = t+1) %>% # Wages are a covariate for the next year's observation.
  pivot_wider(
    names_from = t,
    values_from = lwage,
    names_glue = "lwage_t{sprintf('%02d', t)}"
  )
  ) %>%
  select(-nr) %>%
  as.data.frame()


# Binary treatment matrix: union status each year
A <- wagepan %>%
  select(nr, t, union) %>%
  pivot_wider(
    id_cols   = nr,
    names_from  = t,
    values_from = union,
    names_glue  = "A_t{sprintf('%02d', t)}"          # A_t01 … A_t08
  ) %>%
  arrange(nr) %>%
  select(-nr) %>%
  as.data.frame()

## --------------------------------------------------
## 3. Run your functions
## --------------------------------------------------

props <- estimate_propensity(
  id,
  baseline_X, 
  time_X, 
  A, 
  nsplits = 5,
  fit             = "sl",
  sl.lib          = list(
    "SL.glmnet.interact",
    "SL.lm",
    "SL.glm",
    "SL.mean",
    "SL.ranger",
    "SL.rpart"
  ))

# Run for always treated
always_tx <- run_flip(
  id              = id,
  baseline_X      = baseline_X,
  time_X          = time_X,
  A               = A,
  Y               = Y,
  target_regime   = rep(1, max(props$time_indices)),
  props           = props,
  flip_func       = function(p) 1 - exp(-20 * p),
  flip_func_deriv = function(p) 20 * exp(-20 * p),
  fit             = "sl",
  sl.lib          = list(
    "SL.glmnet.interact",
    "SL.lm",
    "SL.glm",
    "SL.mean",
    "SL.ranger",
    "SL.rpart"
  ),
  outcome_family  = gaussian(),
  type            = "seq_dr",
  verbose         = TRUE
)

# 3) Run for "never treated"
never_tx <- run_flip(
  id              = id,
  baseline_X      = baseline_X,
  time_X          = time_X,
  A               = A,
  Y               = Y,
  target_regime   = rep(0, max(props$time_indices)),  
  props           = props,
  flip_func       = function(p) 1 - exp(-20 * p),
  flip_func_deriv = function(p) 20 * exp(-20 * p),
  fit             = "sl",
  sl.lib          = list(
    "SL.glmnet.interact",
    "SL.lm",
    "SL.glm",
    "SL.mean",
    "SL.ranger",
    "SL.rpart"
  ),
  outcome_family  = gaussian(),
  type            = "seq_dr",
  verbose         = TRUE
)

# 4) Compute point estimate and 95% CI
num_ifs <- always_tx$ifdat$ifvalues - never_tx$ifdat$ifvalues
num_est <- mean(num_ifs)
num_sd  <- sd(num_ifs) / sqrt(length(num_ifs))


save(props, always_tx, never_tx, 
     file = "../output/wage-analysis-main-output.RData")

cat("Propensity score estimates by timepoint:\n",
    paste(capture.output(apply(props$prop_scores[, , 2], 2, summary)), collapse = "\n"),
    "\n")



########################################
### Estimate the change in treatment

always_tx_dist <- list()
for (t in 1:max(props$time_indices)) {
  cat("\nTimepoint: ", t)
  always_tx_dist[[t]] <- estimate_avg_treatment(
    id              = id,
    baseline_X      = baseline_X,
    time_X          = time_X[, props$time_indices < t],
    A               = A[, props$A_indices < t, drop = F],
    final_A         = A[, props$A_indices == t],
    target_regime   = rep(1, t),  
    props           = props,
    flip_func       = function(p) 1 - exp(-20 * p),
    flip_func_deriv = function(p) 20 * exp(-20 * p),
    fit             = "sl",
    sl.lib          = list(
      "SL.glmnet.interact",
      "SL.lm",
      "SL.glm",
      "SL.mean",
      "SL.ranger",
      "SL.rpart"
    ),
    outcome_family  = gaussian(),
    type            = "seq_dr",
    verbose         = TRUE
  )
}

never_tx_dist <- list()
for (t in 1:max(props$time_indices)) {
  cat("\nTimepoint: ", t)
  never_tx_dist[[t]] <- estimate_avg_treatment(
    id              = id,
    baseline_X      = baseline_X,
    time_X          = time_X[, props$time_indices < t],
    A               = A[, props$A_indices < t, drop = F],
    final_A         = A[, props$A_indices == t],
    target_regime   = rep(0, t),  
    props           = props,
    flip_func       = function(p) 1 - exp(-20 * p),
    flip_func_deriv = function(p) 20 * exp(-20 * p),
    fit             = "sl",
    sl.lib          = list(
      "SL.glmnet.interact",
      "SL.lm",
      "SL.glm",
      "SL.mean",
      "SL.ranger",
      "SL.rpart"
    ),
    outcome_family  = gaussian(),
    type            = "seq_dr",
    verbose         = TRUE
  )
}


save(always_tx_dist, never_tx_dist, 
     file = "../output/wage-analysis-avg-trtment.RData")

### Combine to calculate average per-timepoint absolute change in treatment
num_timepoints <- max(props$time_indices)
trt_diff_info <- data.frame()
denom_est <- 0
denom_ifs <- rep(0, length(id))
for (t in seq_len(num_timepoints)) {
  # Pull the ifdat data frames
  always_ifvalues <- always_tx_dist[[t]]$ifdat$ifvalues
  never_ifvalues  <- never_tx_dist[[t]]$ifdat$ifvalues
  
  # Extract the ifvalues column
  this_time_est <- abs(mean(always_ifvalues - never_ifvalues)) 
  this_time_ifs <- sign(mean(always_ifvalues - never_ifvalues)) * 
    (always_ifvalues - never_ifvalues)
  
  denom_est <- denom_est + this_time_est / num_timepoints
  denom_ifs <- denom_ifs + this_time_ifs / num_timepoints
  
  # Add time, point estimate, variance estimate
  trt_diff_info %<>% bind_rows(
    data.frame(
      t = t,
      ptest = this_time_est,
      sdest = sd(this_time_ifs) / sqrt(length(this_time_ifs))
    )
  )
}

denom_sd <- sd(denom_ifs) / sqrt(length(denom_ifs))



ratio_ifs <- num_ifs / denom_est - (num_est / denom_est^2) * denom_ifs 
ratio_sd <- sd(ratio_ifs) / sqrt(length(ratio_ifs))



##############################
### Clean and save

# 1) Plot of propensity scores at each timepoint
prop_plot_dat <- as.data.frame(props$prop_scores[,,2])
colnames(prop_plot_dat) <- c(1:num_timepoints)
prop_score_plot <- prop_plot_dat %>%
  pivot_longer(cols = all_of(colnames(prop_plot_dat))) %>%
  ggplot(aes(x = name, y = value)) + 
  geom_boxplot() + 
  theme_clean() +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1)) + 
  labs(x = "Timepoint",
       y = "Propensity score distribution")

ggsave("../figures/wage-prop-score-box.png",
       plot = prop_score_plot,
       width = 8, height = 4)

# 2) Plot with difference in treatments and their 95% CI
trt_diff_plot <- trt_diff_info %>%
  mutate(ci_ll = ptest - qnorm(0.975) * sdest,
         ci_ul = ptest + qnorm(0.975) * sdest) %>%
  ggplot(aes(x = t, y = ptest)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_ll, ymax = ci_ul)) +
  theme_clean() +
  labs(x = "Timepoint",
       y = "Absolute difference in number of treatments")

ggsave("../figures/wage-trt-diff-info.png",
       plot = trt_diff_plot,
       width = 6, height = 4)

# 3) Output the key point estimate information
cat(
  " Mean difference in potential outcomes: ", round(num_est, 3), "\n",
  "95% CI: [", 
  round(num_est - qnorm(0.975) * num_sd, 3), 
  ", ", 
  round(num_est + qnorm(0.975) * num_sd, 3), 
  "]\n",
  "Per-timepoint change in average absolute number of treatments: ", round(denom_est, 3), "\n",
  "95% CI: [", 
  round(denom_est - qnorm(0.975) * denom_sd, 3), 
  ", ", 
  round(denom_est + qnorm(0.975) * denom_sd, 3), 
  "]\n",
  "Ratio: ", round(num_est / denom_est, 3), "\n",
  "95% CI: [", 
  round(num_est / denom_est - qnorm(0.975) * ratio_sd, 3), 
  ", ", 
  round(num_est / denom_est + qnorm(0.975) * ratio_sd, 3), 
  "]\n"
)





