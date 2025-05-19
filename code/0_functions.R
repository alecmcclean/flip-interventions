###############################################################################
### Author: Alec McClean
### Purpose: Functions for longitudinal smooth trimming with censoring
###############################################################################

### Helper functions
#-- Set 0/0 = 0. 
safe_divide <- function(a,b) { x <- a/b; x[is.nan(x)] <- 0; x }

# Helper function to assert that all the time-varying variables follow
# the correct form: they should end in "_XX"
check_colnames_pattern <- function(mat) {
  if (!all(grepl(".*_\\d{2}$", colnames(mat)))) {
    stop("Error: Not all column names end in '_XX' where XX are exactly two digits.")
  }
  invisible(TRUE)  # Return TRUE invisibly if the check passes
}

# --- Helper function for model fitting ---
fit_model <- function(train_data, eval_data, outcome_var, predictors, fit, sl.lib, family) {
  if (fit == "rf") {
    mod <- ranger(as.formula(paste0(outcome_var, " ~ .")), data = train_data)
    preds <- predict(mod, data = eval_data)$predictions
  } else if (fit == "sl") {
    mod <- suppressWarnings(
      SuperLearner(
        Y = train_data[[outcome_var]], 
        X = train_data[predictors],
        newX = eval_data[predictors],
        SL.library = sl.lib,
        family = family,
        verbose = F)
    )
    
    preds <- mod$SL.predict
  }
  
  return(preds)
}

# Validate that time-varying variables have proper column names and extract time indices.
validate_time_columns <- function(data, var_name) {
  cols <- colnames(data)
  valid <- grepl(".*_t\\d{2}$", cols)
  if (!all(valid)) {
    stop(sprintf("Error: Not all column names in %s end with '_tXX' (two digits). Problem columns: %s", 
                 var_name, paste(cols[!valid], collapse = ", ")))
  }
  # Extract the XX part after '_t' for use as time index
  as.numeric(sub(".*_t(\\d{2})$", "\\1", cols))
}

### Helper SL learners
# 1) Interaction screener: pick top mains by correlation
screen.interaction <- function(Y, X, family, top = 5, ...) {
  # 1) pick the top 'top' columns by Pearson screening
  keep_idx <- screen.corP(Y, X, family = family, P = top)
  # 2) build a logical vector of length ncol(X)
  sel <- rep(FALSE, ncol(X))
  sel[keep_idx] <- TRUE
  sel
}

# 2) glmnet interaction learner: expands all pairwise terms internally
SL.glmnet.interact <- function(Y, X, newX, family, obsWeights, ...) {
  form   <- as.formula("~ (.)^2")
  Xmat   <- model.matrix(form, data = X)[, -1]
  newMat <- model.matrix(form, data = newX)[, -1]
  fit    <- glmnet::cv.glmnet(
    x      = Xmat,
    y      = Y,
    weights= obsWeights,
    family = family$family
  )
  pred   <- predict(fit, newx = newMat, s = "lambda.min", type = "response")
  out    <- list(pred = as.vector(pred), object = fit, formula = form)
  class(out) <- "SL.glmnet.interact"
  return(out)
}
  
predict.SL.glmnet.interact <- function(object, newdata, ...) {
  newMat <- model.matrix(object$formula, data = newdata)[, -1]
  predict(object$object, newx = newMat, s = "lambda.min", type = "response")
}


#’ @title Estimate Time‐Varying Propensity Scores via Cross‐Fitting
#’ @description
#’ Performs V‐fold cross‐validation (cross‐fitting) to estimate P(Aₜ=0) and P(Aₜ=1)
#’ for each subject and each timepoint t=1…T.  Uses either ranger (“rf”) or
#’ a SuperLearner ensemble (“sl”).
#’
#’ @param id           Vector of subject IDs (length n).
#’ @param baseline_X   Data.frame of baseline covariates (n × p₀), or NULL.
#’ @param time_X       Data.frame of time‐varying covariates; names end in `_tXX`.
#’ @param A            Data.frame of binary treatments; names end in `_tXX`.
#’ @param nsplits      Integer ≥ 2; number of CV splits.
#’ @param fit          “rf” or “sl”.
#’ @param sl.lib       Character vector of SuperLearner learners (if fit = “sl”).
#’ @param tx_family    glm family for the treatment model (e.g. binomial()).
#’ @param verbose      Logical; show progress bar if TRUE.
#’
#’ @return A list with:
#’   \item{prop_scores  }{n × T × 2 array of P(Aₜ=0) and P(Aₜ=1).}
#’   \item{folds        }{length‐n vector of CV fold assignments.}
#’   \item{time_indices }{extracted numeric timepoints from `time_X` colnames.}
#’   \item{A_indices    }{extracted numeric timepoints from `A` colnames.}
#’
#’ @export
estimate_propensity <- function(
    id,
    baseline_X,
    time_X,
    A,
    nsplits   = 2,
    fit       = "rf",
    sl.lib    = c("SL.gam","SL.glm","SL.mean","SL.ranger"),
    tx_family = binomial(),
    verbose   = TRUE
) {
  if (fit == "rf") require(ranger)
  n <- length(id)
  if (nrow(baseline_X) != n & is.null(baseline_X) != 0) stop("Row count of baseline_X does not match length of id.")
  if (nrow(time_X) != n) stop("time_X must have n rows")
  if (nrow(A)      != n) stop("A must have n rows")
  if (!fit %in% c("rf","sl")) stop("fit must be 'rf' or 'sl'")
  
  time_indices   <- validate_time_columns(time_X, "time_X")
  A_indices      <- validate_time_columns(A,      "A")
  num_timepoints <- max(time_indices)
  if (num_timepoints > 99) stop("num_timepoints > 99 not supported")
  
  prop_scores <- array(NA_real_, dim = c(n, num_timepoints, 2))
  folds       <- sample(rep(1:nsplits, length.out = n))
  
  if (verbose) {
    message("Estimating propensity scores…")
    pb <- utils::txtProgressBar(min = 0, max = num_timepoints, style = 3)
  }
  
  for (t in seq_len(num_timepoints)) {
    estdat <- cbind(
      time_X[, time_indices <= t, drop = FALSE],
      A[,      A_indices    <= t, drop = FALSE]
    )
    if (nrow(baseline_X) != 0) {
      estdat <- cbind(estdat, baseline_X)
    } 
    
    trt_var    <- sprintf("A_t%02d", t)
    predictors <- setdiff(colnames(estdat), trt_var)
    
    for (split in seq_len(nsplits)) {
      train_idx <- which(folds != split)
      eval_idx  <- which(folds == split)
      
      train_data <- as.data.frame(estdat[train_idx, , drop = FALSE])
      eval_data  <- as.data.frame(estdat[eval_idx,  , drop = FALSE])
      
      prop_preds <- fit_model(
        train_data  = train_data,
        eval_data   = eval_data,
        outcome_var = trt_var,
        predictors  = predictors,
        fit         = fit,
        sl.lib      = sl.lib,
        family      = tx_family
      )
      
      prop_scores[eval_idx, t, 2] <- prop_preds
    }
    
    prop_scores[, t, 1] <- 1 - prop_scores[, t, 2]
    
    if (verbose) utils::setTxtProgressBar(pb, t)
  }
  
  if (verbose) close(pb)
  
  list(
    prop_scores  = prop_scores,
    folds        = folds,
    time_indices = time_indices,
    A_indices    = A_indices
  )
}


#’ @title Run sequential regressions & influence function estimation for  
#' flip effects
#’ @description
#’ Given pre‐computed propensity scores (`props`), fits cross‐fitted pseudo‐outcome
#’ regressions at each timepoint under both A=0 and A=1, constructs plug‐in or
#’ sequential DR‐augmented estimates, and returns influence function values.
#’
#’ @param id             Vector of subject IDs (length n).
#’ @param baseline_X     Data.frame of baseline covariates (n × p₀), or NULL.
#’ @param time_X         Data.frame of time‐varying covariates; names end in `_tXX`.
#’ @param A              Data.frame of binary treatments; names end in `_tXX`.
#’ @param Y              Outcome vector of length n.
#’ @param target_regime  Integer vector (0/1) of length T.
#’ @param props          List output from `estimate_propensity()`, containing:
#’                        - `prop_scores` (n × T × 2 array),
#’                        - `folds` (length‐n CV folds),
#’                        - `time_indices`, `A_indices`.
#’ @param flip_func      Function mapping a propensity score to flip probability.
#’ @param flip_func_deriv Function giving derivative of `flip_func`.
#’ @param fit            `"rf"` or `"sl"` (same as in `estimate_propensity`).
#’ @param sl.lib         Character vector of SuperLearner learners (if `fit = "sl"`).
#’ @param outcome_family glm family for pseudo‐outcome regressions (e.g. `gaussian()`).
#’ @param type           `"mult"` for plug‐in only, or `"seq_dr"` for sequential DR.
#’ @param verbose        Logical; if `TRUE`, shows a progress bar.
#’
#’ @return A list with components:
#’   \describe{
#’     \item{ifdat}{data.frame with `id` and the final influence‐function values.}
#’     \item{prop_scores}{matrix (n × T) of P(Aₜ = 1).}
#’     \item{seq_regs}{array (n × T × 2) of cross‐fitted predictions under A=0,1.}
#’     \item{q_scores}{array (n × T × 2) of interventional propensity scores.}
#’     \item{phi_scores}{array (n × T × 2) of EIF adjustments for DR.}
#’   }
#’ @export
run_flip <- function(
    id,
    baseline_X,
    time_X,
    A,
    Y,
    target_regime,
    props,
    flip_func,
    flip_func_deriv,
    fit,
    sl.lib,
    outcome_family,
    type    = "mult",
    verbose = TRUE
) {
  n              <- length(id)
  prop_scores    <- props$prop_scores
  folds          <- props$folds
  splits <- unique(folds)
  time_indices   <- props$time_indices
  A_indices      <- props$A_indices
  num_timepoints <- length(target_regime)
  
  # Pre‐allocate
  q_scores   <- array(NA_real_, c(n, num_timepoints, 2))
  phi_scores <- array(NA_real_, c(n, num_timepoints, 2))
  ratios     <- matrix(NA_real_, n, num_timepoints)
  seq_regs   <- array(NA_real_, c(n, num_timepoints, 2))
  pseudo     <- NULL
  
  # compute total steps for progress bar
  extra_if <- if (type == "mult") num_timepoints else 0
  total    <- num_timepoints + num_timepoints + extra_if
  
  if (verbose) {
    message("Running flip: building sequential regressions and influence functions.")
    pb   <- utils::txtProgressBar(min = 0, max = total, style = 3)
    step <- 0
  }
  
  # 1) Build interventional propensity scores & ratios
  for (t in seq_len(num_timepoints)) {
    target_tx    <- target_regime[t]
    flip_p       <- flip_func(prop_scores[, t, target_tx + 1])
    flip_p_deriv <- flip_func_deriv(prop_scores[, t, target_tx + 1])
    
    q_scores[, t, 2] <- prop_scores[, t, 2] * (1 - flip_p) + (target_tx == 1) * flip_p
    q_scores[, t, 1] <- 1 - q_scores[, t, 2]
    
    trt_var <- sprintf("A_t%02d", t)
    A_t     <- A[[trt_var]]

    ratios[, t] <- 
      safe_divide(q_scores[, t, ][cbind(seq_len(n), A_t + 1)], 
                  prop_scores[, t, ][cbind(seq_len(n), A_t + 1)])

    # phi_scores
    for (lvl in 0:1) {
      phi_scores[, t, lvl + 1] <- (2 * (lvl == target_tx) - 1) *
        ((A_t == target_tx) - prop_scores[, t, target_tx + 1]) *
        (1 - flip_p + flip_p_deriv * (1 - prop_scores[, t, target_tx + 1]))
    }

    if (verbose) {
      step <- step + 1
      utils::setTxtProgressBar(pb, step)
    }
  }
  
  # 2) Sequential regressions
  for (t in num_timepoints:1) {
    estdat <- cbind(
      time_X[, time_indices <= t, drop = FALSE],
      A[,      A_indices    <= t, drop = FALSE]
    )
    if (nrow(baseline_X) != 0) {
      estdat <- cbind(estdat, baseline_X)
    }     
    if (t == num_timepoints) pseudo <- as.vector(Y)
    outdat      <- cbind(estdat, pseudo)
    outcome_var <- "pseudo"
    predictors  <- setdiff(colnames(estdat), outcome_var)
    trt_var     <- sprintf("A_t%02d", t)
    
    for (split in splits) {
      train_idx <- which(folds != split)
      eval_idx  <- which(folds == split)
      
      train_data <- outdat[train_idx, , drop = FALSE] %>%
        as.data.frame() %>%
        dplyr::select(all_of(c(predictors, outcome_var)))
      eval_base  <- outdat[eval_idx, , drop = FALSE] %>%
        as.data.frame() %>%
        dplyr::select(all_of(predictors))
      
      tmp0 <- eval_base; tmp0[, trt_var] <- 0
      tmp1 <- eval_base; tmp1[, trt_var] <- 1
      
      preds0 <- fit_model(
        train_data   = train_data,
        eval_data    = tmp0,
        outcome_var  = outcome_var,
        predictors   = predictors,
        fit          = fit,
        sl.lib       = sl.lib,
        family       = outcome_family
      )
      preds1 <- fit_model(
        train_data   = train_data,
        eval_data    = tmp1,
        outcome_var  = outcome_var,
        predictors   = predictors,
        fit          = fit,
        sl.lib       = sl.lib,
        family       = outcome_family
      )
      
      seq_regs[eval_idx, t, 1] <- preds0
      seq_regs[eval_idx, t, 2] <- preds1
    }
    
    # Plug-in estimate for pseudo-outcome
    pseudo <- rowSums(seq_regs[, t, ] * q_scores[, t, ])
    
    if (type == "seq_dr") {
      
      pseudo <- pseudo + rowSums(seq_regs[, t, ] * phi_scores[, t, ])
      
      for (s in t:num_timepoints) {
        
        for (k in t:s) {
          if (k == t) {
            ratio_product <- ratios[, k]
            
          } else {
            ratio_product <- ratio_product * ratios[, k]
          }
        }
        
        if (s < num_timepoints) {
          resid <- 
            rowSums(seq_regs[, s+1, ] * (q_scores[, s+1, ] + phi_scores[, s+1, ])) -
            seq_regs[, s, ][cbind(1:n, A[, sprintf("A_t%02d", s)] + 1)]
        } else {
          resid <- as.vector(Y) - seq_regs[, s, ][cbind(1:n, A[, sprintf("A_t%02d", s)]+1)]
        }
        
        pseudo <- pseudo + ratio_product * resid
      }
    }
    
    if (verbose) {
      step <- step + 1
      utils::setTxtProgressBar(pb, step)
    }
  }
  
  # 3) Only if type = "mult" and plug-in pseudo-outcome was used
  if (type == "mult") {
    plugin      <- mean(pseudo)
    centered_if <- rowSums(seq_regs[, 1, ] * q_scores[, 1, ]) - plugin
    
    for (t in 2:(num_timepoints + 1)) {
      rp <- rep(1, n)
      for (s in 1:(t - 1)) rp <- rp * ratios[, s]
      
      resid <- if (t <= num_timepoints) {
        rowSums(seq_regs[, t, ] * (q_scores[, t, ] + phi_scores[, t, ])) -
          seq_regs[, t - 1, ][cbind(seq_len(n), A[, sprintf("A_t%02d", t - 1)] + 1)]
      } else {
        as.vector(Y) -
          seq_regs[, num_timepoints, ][cbind(seq_len(n), A[, sprintf("A_t%02d", num_timepoints)] + 1)]
      }
      
      centered_if <- centered_if + rp * resid
      
      if (verbose) {
        step <- step + 1
        utils::setTxtProgressBar(pb, step)
      }
    }
    
    ifvalues <- plugin + centered_if
  } else {
    # for seq_dr, pseudo is already the IF
    ifvalues <- pseudo
  }
  
  if (verbose) close(pb)
  
  list(
    ifdat       = data.frame(id = id, ifvalues = ifvalues),
    prop_scores = prop_scores[, , 2],
    seq_regs    = seq_regs,
    q_scores    = q_scores,
    phi_scores  = phi_scores
  )
}



#’ @title Run sequential regressions & influence function estimation for  
#' average treatments
#’ @description
#’ Given pre‐computed propensity scores (`props`), fits cross‐fitted pseudo‐outcome
#’ regressions at each timepoint under both A=0 and A=1, constructs plug‐in or
#’ sequential DR‐augmented estimates, and returns influence function values.
#’
#’ @param id             Vector of subject IDs (length n).
#’ @param baseline_X     Data.frame of baseline covariates (n × p₀), or NULL.
#’ @param time_X         Data.frame of time‐varying covariates; names end in `_tXX`.
#’ @param A              Data.frame of binary treatments; names end in `_tXX`.
#’ @param Y              Outcome vector of length n.
#’ @param target_regime  Integer vector (0/1) of length T.
#’ @param props          List output from `estimate_propensity()`, containing:
#’                        - `prop_scores` (n × T × 2 array),
#’                        - `folds` (length‐n CV folds),
#’                        - `time_indices`, `A_indices`.
#’ @param flip_func      Function mapping a propensity score to flip probability.
#’ @param flip_func_deriv Function giving derivative of `flip_func`.
#’ @param fit            `"rf"` or `"sl"` (same as in `estimate_propensity`).
#’ @param sl.lib         Character vector of SuperLearner learners (if `fit = "sl"`).
#’ @param outcome_family glm family for pseudo‐outcome regressions (e.g. `gaussian()`).
#’ @param type           `"mult"` for plug‐in only, or `"seq_dr"` for sequential DR.
#’ @param verbose        Logical; if `TRUE`, shows a progress bar.
#’
#’ @return A list with components:
#’   \describe{
#’     \item{ifdat}{data.frame with `id` and the final influence‐function values.}
#’     \item{prop_scores}{matrix (n × T) of P(Aₜ = 1).}
#’     \item{seq_regs}{array (n × T × 2) of cross‐fitted predictions under A=0,1.}
#’     \item{q_scores}{array (n × T × 2) of interventional propensity scores.}
#’     \item{phi_scores}{array (n × T × 2) of EIF adjustments for DR.}
#’   }
#’ @export
estimate_avg_treatment <- function(
    id,
    baseline_X,
    time_X,
    A,
    final_A, # final treatment value
    target_regime,
    props,
    flip_func,
    flip_func_deriv,
    fit,
    sl.lib,
    outcome_family,
    type    = "mult",
    verbose = TRUE
) {

  num_timepoints <- length(target_regime) - 1
  final_target <- target_regime[num_timepoints + 1]
  
  ### Get final timepoint q_scores and phi_scores
  final_prop_scores <- props$prop_scores[, num_timepoints+1,]
  flip_p       <- flip_func(final_prop_scores[, final_target + 1])
  flip_p_deriv <- flip_func_deriv(final_prop_scores[, final_target + 1])
  
  # Update here: these need to take into account the target treatment
  final_q_scores <- final_prop_scores[,2] * (1 - flip_p) +
    (final_target == 1) * flip_p * (1 - final_prop_scores[,2])

  # phi_scores
  final_phi_scores <- ((Y == 1) - final_prop_scores[,2]) *
      (1 - flip_p + flip_p_deriv * (1 - final_prop_scores[,2]))
    
  # If there's only one timepoint, just return
  if (num_timepoints == 0) {
    return(
      list(
        ifdat       = data.frame(id = id, ifvalues = final_q_scores + final_phi_scores),
        prop_scores = final_prop_scores,
        q_scores    = final_q_scores,
        phi_scores  = final_phi_scores
      )
    )
  }
  
  time_indices   <- props$time_indices[props$time_indices <= num_timepoints]
  A_indices      <- props$A_indices[props$A_indices <= num_timepoints]
  n              <- length(id)
  prop_scores    <- props$prop_scores
  folds          <- props$folds
  splits <- unique(folds)
  
  # Pre‐allocate
  q_scores   <- array(NA_real_, c(n, num_timepoints, 2))
  phi_scores <- array(NA_real_, c(n, num_timepoints, 2))
  ratios     <- matrix(NA_real_, n, num_timepoints)
  seq_regs   <- array(NA_real_, c(n, num_timepoints, 2))
  pseudo     <- NULL
  
  # compute total steps for progress bar
  extra_if <- if (type == "mult") num_timepoints else 0
  total    <- num_timepoints + num_timepoints + extra_if
  
  if (verbose) {
    message("Running flip: building sequential regressions and influence functions.")
    pb   <- utils::txtProgressBar(min = 0, max = total, style = 3)
    step <- 0
  }
  
  # 1) Build interventional propensity scores & ratios
  for (t in seq_len(num_timepoints)) {
    target_tx    <- target_regime[t]
    flip_p       <- flip_func(prop_scores[, t, target_tx + 1])
    flip_p_deriv <- flip_func_deriv(prop_scores[, t, target_tx + 1])
    
    q_scores[, t, 2] <- prop_scores[, t, 2] * (1 - flip_p) + (target_tx == 1) * flip_p
    q_scores[, t, 1] <- 1 - q_scores[, t, 2]
    
    trt_var <- sprintf("A_t%02d", t)
    A_t     <- A[[trt_var]]

    ratios[, t] <- 
      safe_divide(q_scores[, t, ][cbind(seq_len(n), A_t + 1)], 
                  prop_scores[, t, ][cbind(seq_len(n), A_t + 1)])

    # phi_scores
    for (lvl in 0:1) {
      phi_scores[, t, lvl + 1] <- (2 * (lvl == target_tx) - 1) *
        ((A_t == target_tx) - prop_scores[, t, target_tx + 1]) *
        (1 - flip_p + flip_p_deriv * (1 - prop_scores[, t, target_tx + 1]))
    }

    if (verbose) {
      step <- step + 1
      utils::setTxtProgressBar(pb, step)
    }
  }
  
  # 2) Sequential regressions
  for (t in num_timepoints:1) {
    estdat <- cbind(
      time_X[, time_indices <= t, drop = FALSE],
      A[,      A_indices    <= t, drop = FALSE]
    )
    if (nrow(baseline_X) != 0) {
      estdat <- cbind(estdat, baseline_X)
    }     

    if (t == num_timepoints) pseudo <- final_q_scores
    if (type == "seq_dr") pseudo <- pseudo + final_phi_scores
    
    outdat      <- cbind(estdat, pseudo)
    outcome_var <- "pseudo"
    predictors  <- setdiff(colnames(estdat), outcome_var)
    trt_var     <- sprintf("A_t%02d", t)
    
    for (split in splits) {
      train_idx <- which(folds != split)
      eval_idx  <- which(folds == split)
      
      train_data <- outdat[train_idx, , drop = FALSE] %>%
        as.data.frame() %>%
        dplyr::select(all_of(c(predictors, outcome_var)))
      eval_base  <- outdat[eval_idx, , drop = FALSE] %>%
        as.data.frame() %>%
        dplyr::select(all_of(predictors))
      
      tmp0 <- eval_base; tmp0[, trt_var] <- 0
      tmp1 <- eval_base; tmp1[, trt_var] <- 1
      
      preds0 <- fit_model(
        train_data   = train_data,
        eval_data    = tmp0,
        outcome_var  = outcome_var,
        predictors   = predictors,
        fit          = fit,
        sl.lib       = sl.lib,
        family       = outcome_family
      )
      preds1 <- fit_model(
        train_data   = train_data,
        eval_data    = tmp1,
        outcome_var  = outcome_var,
        predictors   = predictors,
        fit          = fit,
        sl.lib       = sl.lib,
        family       = outcome_family
      )
      
      seq_regs[eval_idx, t, 1] <- preds0
      seq_regs[eval_idx, t, 2] <- preds1
    }
    
    # Plug-in estimate for pseudo-outcome
    pseudo <- rowSums(seq_regs[, t, ] * q_scores[, t, ])
    
    if (type == "seq_dr") {
      
      pseudo <- pseudo + rowSums(seq_regs[, t, ] * phi_scores[, t, ])
      
      for (s in t:num_timepoints) {
        
        for (k in t:s) {
          if (k == t) {
            ratio_product <- ratios[, k]
            
          } else {
            ratio_product <- ratio_product * ratios[, k]
          }
        }
        
        if (s < num_timepoints) {
          resid <- 
            rowSums(seq_regs[, s+1, ] * (q_scores[, s+1, ] + phi_scores[, s+1, ])) -
            seq_regs[, s, ][cbind(1:n, A[, sprintf("A_t%02d", s)] + 1)]
        } else {
          resid <- final_q_scores + final_phi_scores -
            seq_regs[, s, ][cbind(1:n, A[, sprintf("A_t%02d", s)]+1)]
        }
        
        pseudo <- pseudo + ratio_product * resid
      }
    }
    
    if (verbose) {
      step <- step + 1
      utils::setTxtProgressBar(pb, step)
    }
  }
  
  # 3) Only if type = "mult" and plug-in pseudo-outcome was used
  if (type == "mult") {
    plugin      <- mean(pseudo)
    centered_if <- rowSums(seq_regs[, 1, ] * q_scores[, 1, ]) - plugin
    
    for (t in 2:(num_timepoints + 1)) {
      rp <- rep(1, n)
      for (s in 1:(t - 1)) rp <- rp * ratios[, s]
      
      resid <- if (t <= num_timepoints) {
        rowSums(seq_regs[, t, ] * (q_scores[, t, ] + phi_scores[, t, ])) -
          seq_regs[, t - 1, ][cbind(seq_len(n), A[, sprintf("A_t%02d", t - 1)] + 1)]
      } else {
        final_q_scores + final_phi_scores -
          seq_regs[, num_timepoints, ][cbind(seq_len(n), A[, sprintf("A_t%02d", num_timepoints)] + 1)]
      }
      
      centered_if <- centered_if + rp * resid
      
      if (verbose) {
        step <- step + 1
        utils::setTxtProgressBar(pb, step)
      }
    }
    
    ifvalues <- plugin + centered_if
  } else {
    # for seq_dr, pseudo is already the IF
    ifvalues <- pseudo
  }
  
  if (verbose) close(pb)
  
  list(
    ifdat       = data.frame(id = id, ifvalues = ifvalues),
    prop_scores = prop_scores[, , 2],
    seq_regs    = seq_regs,
    q_scores    = q_scores,
    phi_scores  = phi_scores
  )
}



