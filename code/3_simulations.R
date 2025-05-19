###############################################################################
### Author: Alec McClean
### Purpose: Simulations to illustrate flip estimator
###############################################################################

source("2_simulation_functions.R")   
set.seed(20250507)

############################
### Calculate ground truth

# Generate large dataset
large_data <- generate_data(n = 10^8) 

# Compute estimands using LLN
with(large_data, {
  ## Overlap weights: Always treated
  q1 <- overlap_always(pi1)
  q2 <- overlap_always(pi2)
  summand <- Y * (A2 * safe_divide(q2, pi2) + (1 - A2) * safe_divide(1 - q2, 1 - pi2)) *
    (A1 * safe_divide(q1, pi1) + (1 - A1) * safe_divide(1 - q1, 1 - pi1))
  overlap_always_tx_true <<- mean(summand)
  
  ## Overlap weights: Never treated
  q1 <- overlap_never(pi1)
  q2 <- overlap_never(pi2)
  summand <- Y * (A2 * safe_divide(q2, pi2) + (1 - A2) * safe_divide(1 - q2, 1 - pi2)) *
    (A1 * safe_divide(q1, pi1) + (1 - A1) * safe_divide(1 - q1, 1 - pi1))
  overlap_never_tx_true <<- mean(summand)
  
  ## Smooth trim: Always treated
  q1 <- smooth_trim_always(pi1)
  q2 <- smooth_trim_always(pi2)
  summand <- Y * (A2 * safe_divide(q2, pi2) + (1 - A2) * safe_divide(1 - q2, 1 - pi2)) *
    (A1 * safe_divide(q1, pi1) + (1 - A1) * safe_divide(1 - q1, 1 - pi1))
  smooth_trim_always_tx_true <<- mean(summand)
  
  ## Smooth trim: Never treated
  q1 <- smooth_trim_never(pi1)
  q2 <- smooth_trim_never(pi2)
  summand <- Y * (A2 * safe_divide(q2, pi2) + (1 - A2) * safe_divide(1 - q2, 1 - pi2)) *
    (A1 * safe_divide(q1, pi1) + (1 - A1) * safe_divide(1 - q1, 1 - pi1))
  smooth_trim_never_tx_true <<- mean(summand)
})

overlap_true <- overlap_always_tx_true - overlap_never_tx_true
smooth_trim_true <- smooth_trim_always_tx_true - smooth_trim_never_tx_true

# truth lookup table
truth_tbl <- tibble(
  group = c("overlap", "smooth_trim"),
  truth = c(overlap_true, smooth_trim_true)
)

rm(large_data)
gc()


####################################
### Estimate flip effects 

sample_sizes <- c(5, 10, 20, 40, 80, 150, 300, 600)
alpha_ps     <- c(0.1, 0.3, 0.5)
alpha_ms     <- c(0.1, 0.25, 0.4, 0.5)
NUM_ITERS <- 250

transform_groups <- list(
  overlap = list(
    always = list(q  = overlap_always,
                  IF = overlap_always_IF),
    never  = list(q  = overlap_never,
                  IF = overlap_never_IF)
  ),
  smooth_trim = list(
    always = list(q  = smooth_trim_always,
                  IF = smooth_trim_always_IF),
    never  = list(q  = smooth_trim_never,
                  IF = smooth_trim_never_IF)
  )
)

results_df <- data.frame(
  group        = character(),
  n            = numeric(),
  alpha_ps     = numeric(),
  alpha_m      = numeric(),
  iter = numeric(),
  plugin_diff  = numeric(),
  EIF_se       = numeric(),
  EIF_ci_lower = numeric(),
  EIF_ci_upper = numeric(),
  stringsAsFactors = FALSE
)

for (grp_name in names(transform_groups)) {
  grp <- transform_groups[[grp_name]]
  cat("\nEstimand: ", grp_name)
  for (n in sample_sizes) {
    
    cat("\nSample size: ", n)
    for (alpha_ps_curr in alpha_ps) {
      for (alpha_m in alpha_ms) {
        ests <- list(); eifs <- list()
        
        for (iter in 1:NUM_ITERS) {
          dat_raw <- generate_data(n)
          
          for (regime in c("always", "never")) {
            q_fun  <- grp[[regime]]$q
            IF_fun <- grp[[regime]]$IF
            
            # 1. Calculate the true m1_0 & m1_1
            dat <- dat_raw
            dat$m1_0 <- generate_m1(dat$pi2_0, dat$m2_00, dat$m2_01, q_fun)
            dat$m1_1 <- generate_m1(dat$pi2_1, dat$m2_10, dat$m2_11, q_fun)
            
            # 2. "Estimate" nuisances by adding noise
            dat <- estimate_nuisances(dat, alpha_ps_curr, alpha_m)
            
            # 3. plugin & EIF
            ests[[regime]] <- compute_plugin(dat, q_fun)
            eifs[[regime]] <- compute_eif(dat, q_fun, IF_fun)
          }
          
          # 4. contrast + CI
          plugin_diff <- ests$always - ests$never
          diff_if     <- eifs$always - eifs$never
          one_step    <- plugin_diff + mean(diff_if)
          se          <- sd(diff_if) / sqrt(n)
          ci_lower    <- one_step - 1.96*se
          ci_upper    <- one_step + 1.96*se
          
          # 5. record
          results_df <- rbind(
            results_df,
            data.frame(
              group        = grp_name,
              n            = n,
              alpha_ps     = alpha_ps_curr,
              alpha_m      = alpha_m,
              iter = iter,
              one_step     = one_step,
              EIF_se       = se,
              EIF_ci_lower = ci_lower,
              EIF_ci_upper = ci_upper,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }
}

coverage_df <- results_df %>%
  left_join(truth_tbl, by = "group") %>%
  group_by(group, n, alpha_ps, alpha_m) %>%
  summarise(
    coverage = mean(truth >= EIF_ci_lower & truth <= EIF_ci_upper),
    .groups   = "drop"
  )

coverage_plot_df <- coverage_df %>%
  mutate(
    se_cov  = sqrt(coverage * (1 - coverage) / NUM_ITERS),   
    ci_low  = pmax(0, coverage - 1.96 * se_cov),
    ci_high = pmin(1, coverage + 1.96 * se_cov),
    ## rename the estimands
    group   = recode(group,
                     overlap      = "Overlap weights",
                     smooth_trim  = "Smooth trim weights")
  )

## custom facet labellers
ps_lab  <- function(x) sprintf("PS rate:  n^{-%s}", x)   # rows
m_lab   <- function(x) sprintf("m rate:  n^{-%s}",  x)   # columns

custom_labeller <- labeller(
  alpha_ps = ps_lab,
  alpha_m  = m_lab
)

### Create plot
p <- ggplot(coverage_plot_df,
       aes(x = n, y = coverage, colour = group, group = group)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.05) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous(breaks = sample_sizes,
                     trans = "log10",
                     minor_breaks = NULL) +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_colour_manual(values = c("Overlap weights" = "#1b9e77",
                                 "Smooth trim weights" = "#d95f02")) +
  labs(
    x = "Sample size (log scale)",
    y = "Coverage",
    colour = NULL
  ) +
  facet_grid(alpha_ps ~ alpha_m, labeller = custom_labeller) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.spacing    = unit(1, "lines"),
    strip.text.x     = element_text(size = 10),
    strip.text.y     = element_text(size = 10)
  )

ggsave("../figures/simulations-results.png", width = 8, height = 6)
rm(list = ls(all = T))
gc()
