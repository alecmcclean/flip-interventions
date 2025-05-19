###############################################################################
### 3a_simulation_functions.R
###  – All transformations, nuisances, EIFs, data‐generation, noise, helpers
###  - for the simulations.
###############################################################################

logit <- function(x) log(x / (1-x))
expit <- function(x) 1 / (1 + exp(-x))

#-- intervention weights and derivatives
overlap_wt      <- function(x) x * (1 - x)
overlap_wt_diff <- function(x) 1 - 2*x

#-- overlap interventions + centered EIF
overlap_always    <- function(pi) pi + overlap_wt(pi)*(1 - pi)
overlap_always_IF <- function(pi, b, A) {
  (2*(b==1)-1)*((A==1)-pi)*(1 - overlap_wt(pi) + overlap_wt_diff(pi)*(1 - pi))
}

overlap_never    <- function(pi) pi - overlap_wt(1 - pi)*pi
overlap_never_IF <- function(pi, b, A) {
  (2*(b==0)-1)*((A==0)-(1-pi))*(1 - overlap_wt(1-pi) + overlap_wt_diff(1-pi)*pi)
}

#-- smooth-trim interventions + centered EIF
smooth_trim_wt      <- function(x)    (1 - exp(-10*x))*(1 - exp(-10*(1-x)))
smooth_trim_wt_diff <- function(x) 10*(exp(-10*x) - exp(-10*(1-x)))

smooth_trim_always    <- function(pi) pi + smooth_trim_wt(pi)*(1 - pi)
smooth_trim_always_IF <- function(pi, b, A) {
  (2*(b==1)-1)*((A==1)-pi)*(1 - smooth_trim_wt(pi) + smooth_trim_wt_diff(pi)*(1 - pi))
}

smooth_trim_never    <- function(pi) pi - smooth_trim_wt(1 - pi)*pi
smooth_trim_never_IF <- function(pi, b, A) {
  (2*(b==0)-1)*((A==0)-(1-pi))*(1 - smooth_trim_wt(1-pi) + smooth_trim_wt_diff(1-pi)*pi)
}

#-- Data generation (with all potentials)
generate_data <- function(n) {
  X1  <- runif(n)
  pi1 <- ifelse(X1 < 0.1, 0, ifelse(X1 > 0.9, 1, (X1-0.1)/0.8))
  A1  <- rbinom(n,1,pi1)
  X2_0 <- 0.5*X1;    X2_1 <- 0.5*(X1+1)
  pi2_0 <- ifelse(X2_0 < 0.1, 0, ifelse(X2_0 > 0.9, 1, (X2_0-0.1)/0.8))
  pi2_1 <- ifelse(X2_1 < 0.1, 0, ifelse(X2_1 > 0.9, 1, (X2_1-0.1)/0.8))
  A2_0 <- rbinom(n,1,pi2_0); A2_1 <- rbinom(n,1,pi2_1)
  m2_00 <- X1 + 0 + X2_0 + 0
  m2_01 <- X1 + 0 + X2_0 + 1
  m2_10 <- X1 + 1 + X2_1 + 0
  m2_11 <- X1 + 1 + X2_1 + 1
  X2  <- ifelse(A1==1, X2_1, X2_0)
  pi2 <- ifelse(A1==1, pi2_1, pi2_0)
  A2  <- ifelse(A1==1, A2_1, A2_0)
  m2  <- ifelse(A1==1,
                ifelse(A2==1, m2_11, m2_10),
                ifelse(A2==1, m2_01, m2_00))
  Y   <- m2 + rnorm(n,0,1)
  data.frame(X1,pi1,A1,
             X2_0,X2_1,pi2_0,pi2_1,A2_0,A2_1,
             m2_00,m2_01,m2_10,m2_11,
             X2,pi2,A2,m2,Y)
}

#-- Sequential regression at t=1
generate_m1 <- function(pi2, m2_0, m2_1, q_fun) {
  q2 <- q_fun(pi2)
  q2*m2_1 + (1-q2)*m2_0
}

#-- Estimate the outcome models
generate_outcome_noise <- function(n, alpha=1/3) {
  if (alpha == 0) return(rep(0,n))
  3 * rnorm(n, mean = n^(-1 * alpha), sd = n^(-0.5 * alpha))
}

#-- Estimate the propensity score models
generate_treatment_noise <- function(pi, alpha=1/3) {
  if (alpha == 0) return(pi)
  n <- length(pi)
  varepsilon <- 1e-06
  pi <- pmax(1e-06, pmin(pi, 1-1e-06))
  return(
    expit(
      logit(pi) + 2 * rnorm(n, mean = n^(-1 * alpha), sd = n^(-0.5 * alpha))
    )
  )
}

#-- Set 0/0 = 0. 
safe_divide <- function(a,b) { x <- a/b; x[is.nan(x)] <- 0; x }

#-- Estimate nuisance functions
estimate_nuisances <- function(dat, alpha_ps, alpha_m) {
  n <- nrow(dat)
  dat$pi1hat   <- generate_treatment_noise(dat$pi1 ,alpha_ps)
  dat$pi2hat_0 <- generate_treatment_noise(dat$pi2_0, alpha_ps)
  dat$pi2hat_1 <- generate_treatment_noise(dat$pi2_1, alpha_ps)
  dat$pi2hat   <- ifelse(dat$A1==1,dat$pi2hat_1,dat$pi2hat_0)
  dat$m1hat_0  <- dat$m1_0 + generate_outcome_noise(n,alpha_m)
  dat$m1hat_1  <- dat$m1_1 + generate_outcome_noise(n,alpha_m)
  dat$m2hat_00 <- dat$m2_00 + generate_outcome_noise(n,alpha_m)
  dat$m2hat_01 <- dat$m2_01 + generate_outcome_noise(n,alpha_m)
  dat$m2hat_10 <- dat$m2_10 + generate_outcome_noise(n,alpha_m)
  dat$m2hat_11 <- dat$m2_11 + generate_outcome_noise(n,alpha_m)
  return(dat)
}

#-- Centered EIF over two timepoints
compute_eif <- function(dat, q_fun, IF_fun) {
  
  # Interventional prop scores and ratios
  q1 <- q_fun(dat$pi1hat);  q2 <- q_fun(dat$pi2hat)
  r1 <- ifelse(dat$A1==1, safe_divide(q1,dat$pi1hat),  safe_divide(1-q1,1-dat$pi1hat))
  r2 <- ifelse(dat$A2==1, safe_divide(q2,dat$pi2hat),  safe_divide(1-q2,1-dat$pi2hat))
  
  # Pseudo-outcome estimates
  pseudo0 <- q1*dat$m1hat_1 + (1-q1)*dat$m1hat_0
  pseudo1 <- ifelse(dat$A1==1,
                    q2*dat$m2hat_11+(1-q2)*dat$m2hat_10,
                    q2*dat$m2hat_01+(1-q2)*dat$m2hat_00)
  pseudo2 <- dat$Y
  
  # observed nuisances
  seq_reg0 <- mean(q1*dat$m1hat_1 + (1-q1)*dat$m1hat_0)
  seq_reg1 <- ifelse(dat$A1==1, dat$m1hat_1, dat$m1hat_0)
  seq_reg2 <- ifelse(dat$A1==1, ifelse(dat$A2==1,dat$m2hat_11,dat$m2hat_10),
                     ifelse(dat$A2==1,dat$m2hat_01,dat$m2hat_00))
  
  # sequential regression residuals
  resid2 <- r1*r2*(pseudo2 - seq_reg2)
  resid1 <- r1*(pseudo1 - seq_reg1)
  resid0 <- pseudo0 - seq_reg0
  
  
  # Estimating Q residuals
  phi1_b1  <- IF_fun(dat$pi1hat, b=1, A=dat$A1)
  phi1_b0  <- IF_fun(dat$pi1hat, b=0, A=dat$A1)
  term_phi1 <- dat$m1hat_1*phi1_b1 + dat$m1hat_0*phi1_b0
  phi2_b1   <- IF_fun(dat$pi2hat, b=1, A=dat$A2)
  phi2_b0   <- IF_fun(dat$pi2hat, b=0, A=dat$A2)
  m2_b1     <- ifelse(dat$A1==1,dat$m2hat_11,dat$m2hat_01)
  m2_b0     <- ifelse(dat$A1==1,dat$m2hat_10,dat$m2hat_00)
  term_phi2 <- r1*(m2_b1*phi2_b1 + m2_b0*phi2_b0)
  
  return(resid0 + resid1 + resid2 + term_phi1 + term_phi2)
}

#-- Point‐estimate (plugin) for contrast under one transform
compute_plugin <- function(dat, q_fun) {
  q1 <- q_fun(dat$pi1hat)
  mean(q1 * dat$m1hat_1 + (1 - q1) * dat$m1hat_0)
}
