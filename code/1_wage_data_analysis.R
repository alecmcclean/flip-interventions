###############################################################################
### Author: Alec McClean
### Purpose: Clean and analyze wage panel data; last run on 2025/06/09
###############################################################################

set.seed(20250513)
options(scipen = 999)

## --------------------------------------------------
## 0. Raw data and load package
## --------------------------------------------------

data(wagepan)                              # 545 workers × 8 years (1980-87)
wagepan <- wagepan %>% arrange(nr, year)   # keep rows in tidy panel order
wagepan <- wagepan %>% filter(year <= 1983) # Just use 4 years (1980-'83) to simplify analysis

## Add a numeric time index t = 1 … 3 (1980 --> 1, ..., 1983 --> 4)
wagepan <- wagepan %>%
  mutate(t = year - min(year) + 1)

# Load development branch from lmtp
devtools::install_github("alecmcclean/lmtp@flip", force = T)
library(lmtp)

## --------------------------------------------------
## 1. Clean data for analysis
## --------------------------------------------------

# 1a) baseline covariates (educ, black, hisp) + id
baseline_df <- wagepan %>%
  filter(year == 1980) %>%
  select(nr, educ, black, hisp)

# 1b) outcome Y in final year (1983)
Y_df <- wagepan %>%
  filter(year == 1983) %>%
  select(nr, Y = lwage)

# 1c) treatment matrix A_t01, ..., A_t04 (here 4 years)
A_df <- wagepan %>% filter(year <= 1983) %>%
  select(nr, t, union) %>%
  pivot_wider(
    names_from   = t,
    values_from  = union,
    names_prefix = "A_t"
  )

# 1d) time‐varying covariates L*_t01,..., L*_t04
L_df <- wagepan %>%
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
L_df %<>% left_join(
  wagepan %>% 
    filter(year < 1983) %>%
    select(nr, t, lwage) %>% 
    mutate(t = t+1) %>% # Wages are a covariate for the next year's observation.
    pivot_wider(
      names_from = t,
      values_from = lwage,
      names_glue = "lwage_t{sprintf('%02d', t)}"
    ),
  by = "nr"
  ) 

# 1e) join all pieces
dat <- baseline_df %>%
  inner_join(A_df, by = "nr") %>%
  inner_join(L_df, by = "nr") %>%
  inner_join(Y_df, by = "nr")


## -----------------------------------------------------------
## 2. Call longitudinal interventional effect estimator
## -----------------------------------------------------------

time_vary_list <- lapply(1:4, function(j) {
  grep(paste0("_t", sprintf("%02d", j), "$"),
       colnames(L_df), value = TRUE)
})

res <- life_sdr(
  data               = dat,
  id                 = "nr",           # your subject ID column
  trt                = setdiff(colnames(A_df), "nr"),      
  outcome            = "Y",            
  baseline           = setdiff(colnames(baseline_df), "nr"), 
  time_vary          = time_vary_list,
  cens               = NULL,           # no censoring here
  compete            = NULL,           # no competing events
  overlap            = FALSE,          # or TRUE for overlap‐weighting
  trimming_threshold = 0,            
  smoothing_constant = 20,              
  k                  = Inf,            
  outcome_type       = "continuous",   
  learners_outcome = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.lm", "SL.glm", "SL.rpart"),
  learners_trt     = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.lm", "SL.glm", "SL.rpart"),
  folds              = 5,
  control            = lmtp_control(), 
  num_timepoints     = 4               
)



##############################
### Output results

pt_est <- mean(res$auxiliary_info$overall_num_ifs)
sd_est <- sd(res$auxiliary_info$overall_num_ifs) / sqrt(nrow(dat))
cat("Numerator (difference in POs) estimate: ", pt_est, 
    "\n95% CI: ", pt_est - qnorm(0.975) * sd_est, pt_est + qnorm(0.975) * sd_est)
    
pt_est <- mean(res$auxiliary_info$overall_denom_ifs)
sd_est <- sd(res$auxiliary_info$overall_denom_ifs) / sqrt(nrow(dat))
cat("Denominator (per-timepoint average change in treatments) estimate: ", pt_est, 
    "\n95% CI: ", pt_est - qnorm(0.975) * sd_est, pt_est + qnorm(0.975) * sd_est)

pt_est <- res$estimate@x
sd_est <- sd(res$ifvalues) / sqrt(nrow(dat))
cat("Overall point estimate:ß", pt_est, 
    "\n95% CI: ", pt_est - qnorm(0.975) * sd_est, pt_est + qnorm(0.975) * sd_est)

# 1) Plot of propensity scores at each timepoint
prop_plot_dat <- as.data.frame(res$auxiliary_info$always_treated_avg_PO$prop_scores)
colnames(prop_plot_dat) <- c(1:4)
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
trt_diff_plot <- res$auxiliary_info$trt_diff_info %>%
  mutate(ci_ll = ptest - qnorm(0.975) * sdest,
         ci_ul = ptest + qnorm(0.975) * sdest) %>%
  ggplot(aes(x = t, y = ptest)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_ll, ymax = ci_ul)) +
  theme_clean() +
  labs(x = "Timepoint",
       y = "Mean difference in number of treatments")

ggsave("../figures/wage-trt-diff-info.png",
       plot = trt_diff_plot,
       width = 6, height = 4)


