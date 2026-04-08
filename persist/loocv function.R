loocv_glm_subsets <- function(response, predictors, data, family = "binomial") {
  # ---- Interpret response & predictors ----
  if (is.character(response)) {
    resp_name <- response[1]
  } else {
    resp_name <- deparse(substitute(response))
  }
  
  if (!resp_name %in% names(data)) {
    stop("Response '", resp_name, "' is not a column in `data`.")
  }
  
  preds <- as.character(predictors)
  missing_cols <- setdiff(c(resp_name, preds), names(data))
  if (length(missing_cols) > 0) {
    stop("These variables are missing from `data`: ",
         paste(missing_cols, collapse = ", "))
  }
  
  p <- length(preds)
  if (p < 1L) stop("Need at least one predictor.")
  
  # ---- Handle family argument ----
  fam <- if (is.character(family)) {
    fam_fun <- get(family, mode = "function", envir = parent.frame())
    fam_fun()
  } else {
    family
  }
  
  # ---- Helper: coerce response to numeric for scoring ----
  coerce_y <- function(y) {
    if (is.matrix(y) && ncol(y) == 2) return(y[, 1] / rowSums(y))
    if (is.factor(y)) return(as.numeric(y) - 1L)
    as.numeric(y)
  }
  
  # ---- Helper: LOOCV for a single model formula ----
  loocv_single <- function(formula, dat, family) {
    mf_all <- model.frame(formula, data = dat)
    y_raw  <- model.response(mf_all)
    y      <- coerce_y(y_raw)
    n_obs  <- nrow(mf_all)
    
    preds_vec <- rep(NA_real_, n_obs)
    
    for (i in seq_len(n_obs)) {
      fit_i <- tryCatch(
        glm(formula,
            data   = mf_all[-i, , drop = FALSE],
            family = family),
        error = function(e) NULL
      )
      
      if (!is.null(fit_i)) {
        preds_vec[i] <- predict(
          fit_i,
          newdata = mf_all[i, , drop = FALSE],
          type = "response"
        )
      }
    }
    
    if (all(is.na(preds_vec))) return(NA_real_)
    mean((y - preds_vec)^2, na.rm = TRUE)   # Brier/MSE
  }
  
  # ---- Build all non-empty subsets ----
  subset_list <- unlist(
    lapply(1:p, function(m) combn(preds, m, simplify = FALSE)),
    recursive = FALSE
  )
  
  results <- list()
  k <- 1L
  
  ## ---- 1. Intercept-only model ----
  int_form_str <- paste(resp_name, "~ 1")
  int_form <- as.formula(int_form_str)
  
  int_fit <- tryCatch(
    glm(int_form, data = data, family = fam),
    error = function(e) NULL
  )
  int_aic <- if (is.null(int_fit)) NA_real_ else AIC(int_fit)
  
  cv_int <- loocv_single(int_form, data, fam)
  
  results[[k]] <- data.frame(
    formula      = int_form_str,
    predictors   = "(intercept only)",
    n_predictors = 0,
    model_type   = "intercept",
    loocv_score  = cv_int,
    AIC          = int_aic,
    stringsAsFactors = FALSE
  )
  k <- k + 1L
  
  ## ---- 2. All models: additive AND full interactions for each subset ----
  for (j in seq_along(subset_list)) {
    vars <- subset_list[[j]]
    pred_str <- paste(vars, collapse = ",")
    n_pred   <- length(vars)
    
    ## 2a. Additive model: y ~ x1 + x2 + ...
    rhs_add   <- paste(vars, collapse = " + ")
    form_add_str <- paste(resp_name, "~", rhs_add)
    form_add     <- as.formula(form_add_str)
    
    fit_add <- tryCatch(
      glm(form_add, data = data, family = fam),
      error = function(e) NULL
    )
    aic_add <- if (is.null(fit_add)) NA_real_ else AIC(fit_add)
    
    cv_add <- loocv_single(form_add, data, fam)
    
    results[[k]] <- data.frame(
      formula      = form_add_str,
      predictors   = pred_str,
      n_predictors = n_pred,
      model_type   = "additive",
      loocv_score  = cv_add,
      AIC          = aic_add,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    
    ## 2b. Full interaction model: y ~ x1 * x2 * ...
    rhs_int   <- paste(vars, collapse = " * ")
    form_int_str <- paste(resp_name, "~", rhs_int)
    form_int     <- as.formula(form_int_str)
    
    fit_int <- tryCatch(
      glm(form_int, data = data, family = fam),
      error = function(e) NULL
    )
    aic_int <- if (is.null(fit_int)) NA_real_ else AIC(fit_int)
    
    cv_int <- loocv_single(form_int, data, fam)
    
    results[[k]] <- data.frame(
      formula      = form_int_str,
      predictors   = pred_str,
      n_predictors = n_pred,
      model_type   = "full_interactions",
      loocv_score  = cv_int,
      AIC          = aic_int,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  
  out <- do.call(rbind, results)
  
  # sort by loocv_score
  out[order(out$loocv_score), ]
}
