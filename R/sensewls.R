
#' Sensitivity analysis to omitted confounders for weighted least squares estimators of treatment effects
#'
#' @description This function performs the weighted sensitivity analysis for a weighted least squares
#'  estimator as described in [Wainstein and Hazlett (2025)](https://www.arxiv.org/abs/2508.02954). It returns an object of class `sensewls`
#'  containing a variety of sensitivity statistics. An output object can use the `plot_sensewls()`
#'  function for visualization.
#'
#' @param df A dataframe containing `treatment`, `outcome`, and `covars`.
#' @param treatment A string indicating the column name of the treatment variable.
#' @param outcome A string indicating the column name of the outcome variable.
#' @param covars A character vector of column names of the relevant covariates in the dataframe.
#' @param bounding_covars A list of character vectors containing covariates to benchmark by. May be a character vector if only one set of bounding covariates is desired -- note in this case that `sensewls()` will benchmark on all of the covariates in the provided vector *combined*, and *not* separately.
#' @param kd A numeric vector of considered covariate bound strengths for association with the treatment variable. Default value is `1`.
#' @param ky A numeric vector of considered covariate bound strengths for association with the outcome variable. Default value is `kd`.
#' @param w A weights object, containing either
#' (1) a numeric vector of weights
#' (2) an appropriate weighting function, or
#' (3) a string for one of the preset weighting functions.
#' Three presets are currently available, `prop_score` (inverse propensity score weights, with a logistic regression for the treatment variable), `matchit` (through the `MatchIt` package: https://cran.r-project.org/web/packages/MatchIt/index.html)
#' and `ebal` (through the `ebal` package: https://cran.r-project.org/web/packages/ebal/index.html).
#' @param normalize A single logical, indicating whether or not to normalize the weights (and `semiweights`) when either (1) a numeric vector of weights, or (2) a custom weight function is provided for `w`.
#' Default is `FALSE`. But weights are normalized for all the preset weighting functions.
#' @param semiweights List of semiweights, with the ith entry corresponding
#' to the semiweights formed when constructing weights while omitting the
#' ith set of covariates in the `bounding_covars` list. May be one numeric vector of semiweights if `bounding_covars` is only one character vector of bounding covariates.
#' Defaults to `NULL`. `semiweights` should be `NULL` unless `w` is a numeric vector of specified weights.
#' @param estimand A string indicating the desired estimand. Only required if using preset weights. Can be `"ATT"`, `"ATC"`, or `"ATE"`. Default is `NULL`.
#' @param alpha Number between 0 and 1. Level of significance used for RV_alpha, and (1 - alpha)*100% confidence intervals.
#' Defaulted to alpha = 0.05 (and 95% confidence intervals).
#' @param inference Whether to obtain uncertainty estimates. If set to the logical
#' `"TRUE`, then a non-parametric bootstrap is performed that fixes the weights, and samples them
#' along with the data in each bootstrap sample. If set to the character `"reestimate"`,
#' a non-parametric bootstrap is performed that re-estimates the weights in each bootstrap sample. This is only allowed if
#' a preset weighting function or an appropriate, user-supplied weighting function is chosen for `w`.
#' The percentile method is used for all bootstrap confidence intervals.
#' If set to the logical `FALSE`, then no uncertainty estimates are provided. Default is `TRUE`.
#' @param cluster An optional string indicating a column from the data frame `df` for which to perform a cluster bootstrap if desired.
#' `inference` must be set to `TRUE` or `"reestimate"`. Defaulted to `NULL`.
#' @param B Number of bootstrap samples if used. Defaulted to 500.
#' @param par Boolean for whether bootstrapping should be parallelized through
#' the `parallel` library. Defaulted to `FALSE`.
#' @param ncpus Number of cpus used if `par` is set to `TRUE`. Defaulted to `NULL` and
#' must be supplied if `par` is set to true.
#' @param ... Additional arguments passed to preset weighting functions.
#' Currently supports providing `caliper` and `ratio` to the [MatchIt::matchit()] weighting function
#' used by setting `w = 'matchit'`.
#'
#'
#' @returns
#' A sensewls object containing a variety of sensitivity statistics. The output object can be plotted using the
#'  `plot_sensewls()` function. Significance levels and `q` threshold for robustness values and confidnece intervals
#'  can be changed by applying the `alter_sensewls()` function.
#'
#'
#'
#'
#' @example
#' # Darfur Dataset (from the `sensemakr` package) with the `matchit` weighting function.
#' /dontrun{
#' library(sensemakr)
#' data(darfur)
#'
#' sensewls(
#'    df = darfur,
#'    treatment = "directlyharmed",
#'    outcome = "peacefactor",
#'    covars = c("village",
#'               "female",
#'               "age",
#'               "farmer_dar",
#'               "herder_dar",
#'               "pastvoted",
#'               "hhsize_darfur"),
#'    bounding_covars = "female",
#'    kd = 1:3,
#'    w = "matchit",
#'    estimand = "ATE",
#'    inference = TRUE,
#'    cluster = "village",
#'    par = T,
#'    ncpus = 7,
#'    caliper = 0.1,
#'    ratio = 3)
#' }
#'
#'
#'
#' @import dplyr
#'
#' @export
#'
#'
#'
sensewls <- function(df, treatment, outcome, covars,
                bounding_covars, kd = 1, ky = kd,
                w, semiweights = NULL, normalize = FALSE, estimand = NULL,
                inference=TRUE, alpha = 0.05, cluster = NULL,
                B = 500, par = F, ncpus = NULL,
                ...)
  {

  ##### ERROR/WARNING MESSAGES
  ### check basic argument quality
  # df
  if(
    !is.data.frame(df)
  ){stop("`df` must be a dataframe.\n")}

  # Treatment
  if(
    !is.character(treatment)
  ){stop("`treatment` must be a single character string.\n")}
  if(
    is.character(treatment) & length(treatment)>1
  ){stop("`treatment` must be a single character string.\n")}
  if(
    !(treatment %in% names(df))
  ){stop("`treatment` must be variable name in `df`.\n")}
  if(
    !is.numeric(df[, treatment])
  ){stop("`treatment` variable must be numeric.\n")}

  # outcome
  if(
    !is.character(outcome)
  ){stop("`outcome` must be a single character string.\n")}
  if(
    is.character(outcome) & length(outcome)>1
  ){stop("`outcome` must be a single character string.\n")}
  if(
    !(outcome %in% names(df))
  ){stop("`outcome` must be variable name in `df`.\n")}
  if(
    !is.numeric(df[, outcome])
  ){stop("`outcome` variable must be numeric.\n")}

  # covars
  if(
    !is.character(covars)
  ){stop("`covars` must be a character vector.\n")}
  if(
    sum(!(covars %in% names(df)))>0
  ){stop("`covars` must be variable names in `df`.\n")}

  # bounding covars
  if(
    !(is.character(bounding_covars) | is.list(bounding_covars))
  ){stop("`bounding_covars` must be a character vector, or a list of character vectors.\n")}
  if(is.list(bounding_covars)){
    user_list_bounding_covars <- T
    for(listnum in 1:length(bounding_covars)){
      if(!is.character(bounding_covars[[listnum]])) stop("`bounding_covars` must be a character vector, or a list of character vectors.\n")
      if(
        sum(!(bounding_covars[[listnum]] %in% covars))>0
      ){stop("`bounding_covars` must be contained in `covars`.\n")}
    }
  }
  if(is.character(bounding_covars)){
    user_list_bounding_covars <- F
    if(
      sum(!(bounding_covars %in% covars))>0
    ){stop("`bounding_covars` must be contained in `covars`.\n")}
    bounding_covars <-list(bounding_covars)
  }

  # kd
  if(!is.numeric(kd)){stop("`kd` must be a numeric vector.\n")}
  if(sum(kd<=0)>0){stop("`kd` must only contain positive numbers.\n")}

  # ky
  if(!is.numeric(ky)){stop("`ky` must be a numeric vector.\n")}
  if(sum(ky<=0)>0){stop("`ky` must only contain positive numbers.\n")}

  # alpha
  if(
    !is.numeric(alpha)
  ){stop("`alpha` must be a single number.\n")}
  if(
    is.numeric(alpha) & length(alpha)>1
  ){stop("`alpha` must be a single number.\n")}
  if(alpha<=0 | alpha >=1){stop("`alpha` must be a number between 0 and 1.\n")}

  # weights
  if(!is.numeric(w) & !is.function(w) & !is.character(w)){
    stop("`w` must be a vector of numeric weights, a weighting function, or a preset weight string.\n")
  }
  if(is.numeric(w)){
    if(sum(w<0)>0) stop("`w` must be a vector of positive weights.\n")
    if(length(w)!=nrow(df)) stop("`w` must be the same length as the data.\n")
  }
  if(is.character(w)){
    if(length(w) > 1 | !(w %in% c("prop_score", "matchit", "ebal"))) stop("If preset weights are desired, `w` must be `'prop_score'`, `'matchit'` or `'ebal'`.\n")
  }
  if(is.character(w)){
    unique_values <- sort(unique(df[, treatment]))
    if(length(unique_values)!=2 | min(unique_values)!=0 | max(unique_values)!=1) stop("If preset weights are desired, `treatment` must be binary, with 1=treated and 0=control.\n")
  }

  # normalize
  if(
    !is.logical(normalize)
  ){stop("`normalize` must be a single logical.\n")}
  if(
    is.logical(normalize) & length(normalize)>1
  ){stop("`normalize` must be a single logical.\n")}
  if(
    (is.numeric(w) | is.function(w)) & normalize==TRUE
  ){
    unique_values <- sort(unique(df[, treatment]))
    if(length(unique_values)!=2 | min(unique_values)!=0 | max(unique_values)!=1) stop("In order to normalize weights, `treatment` must be binary, with 1=treated and 0=control.\n")
  }

  # estimand
  if(is.character(w)){
    if(
      !is.character(estimand)
    ){stop("`estimand` must be a single character string.\n")}
    if(
      is.character(estimand) & length(estimand)>1
    ){stop("`estimand` must be a single character string.\n")}
    if(
      !(estimand %in% c("ATT", "ATC", "ATE"))
    ){stop("`estimand` must be one of ATT, ATC, or ATE.\n")}
  }

  # inference
  if(
    !(is.logical(inference) | is.character(inference))
  ){stop("`inference` must be either `TRUE`, `FALSE`, or `'reestimate'`.\n")}
  if(
    is.logical(inference) & length(inference)>1
  ){stop("`inference` must be either `TRUE`, `FALSE`, or `'reestimate'`.\n")}
  if(
    is.character(inference) & length(inference)>1
  ){stop("`inference` must be either `TRUE`, `FALSE`, or `'reestimate'`.\n")}
  if(is.numeric(w) & inference=="reestimate"){
    stop("Bootstrap cannot re-estimate weights if they are specified with `w`. `inference` must be set to `TRUE` or `FALSE`.\n")
  }

  # cluster
  if( !is.null(cluster) & (inference==TRUE | inference=="reestimate") ){ # only check cluster if it's not NULL
    if(
      !is.character(cluster)
    ){stop("`cluster` must be a single character string.\n")}
    if(
      is.character(cluster) & length(cluster)>1
    ){stop("`cluster` must be a single character string.\n")}
    if(
      !(cluster %in% names(df))
    ){stop("`cluster` must be variable name in `df`.\n")}
  }

  # B, par, and cpus
  if(inference==TRUE | inference=="reestimate"){
    if(
      !is.numeric(B)
    ){stop("`B` must be a single number.\n")}
    if(
      is.numeric(B) & length(B)>1
    ){stop("`B` must be a single number.\n")}
    if(
      B < 1 | (B - floor(B) > 1e-10)
    ){stop("`B` must be a positive integer.\n")}
  }

  # par
  if(inference==TRUE | inference=="reestimate"){
    if(
      !is.logical(par)
    ){stop("`par` must be a single logical.\n")}
    if(
      is.logical(par) & length(par)>1
    ){stop("`par` must be a single logical.\n")}
  }

  # ncpus
  if( (inference==TRUE | inference=="reestimate") & par==TRUE ){
    if(
      !is.numeric(ncpus)
    ){stop("`ncpus` must be a single number.\n")}
    if(
      is.numeric(ncpus) & length(ncpus)>1
    ){stop("`ncpus` must be a single number.\n")}
    if(
      ncpus < 1 | (ncpus - floor(ncpus) > 1e-10)
    ){stop("`ncpus` must be a positive integer.\n")}
  }

  # semiweights
  if(is.numeric(w)){
    if(user_list_bounding_covars==F){
      if(is.list(semiweights)){
        if(length(bounding_covars[[1]])==length(semiweights))  stop("`semiweights` must be a numeric vector if `bounding_covars` is input as a character vector.\nNumber of covariates in `bounding_covars` vector is the same as the number of weight vectors in the `semiweights` list. If you meant to benchmark on each variable in `bounding_covars` individually, then `bounding_covars` should be input as a list, with each variable being its own element in the list.\n")
      }
      if(!is.numeric(semiweights)) stop("`semiweights` must be a numeric vector if `bounding_covars` is input as a character vector.\n")
      if(length(w)!=length(semiweights)) stop("`semiweights` must be the same length as the data.\n")
      if(sum(semiweights<0)>0) stop("`semiweights` must be a vector of positive weights.\n")
      semiweights <- list(semiweights)
    }
    if(user_list_bounding_covars==T){
      if(!is.list(semiweights)) stop("`semiweights` must be a list.\n")
      if(length(bounding_covars)!=length(semiweights)) stop("`semiweights` must be a list of the same length as `bounding_covars`.\n")
      for(listnum in 1:length(semiweights)){
        if(!is.numeric(semiweights[[listnum]])) stop("Weight vectors in `semiweights` must be numeric vectors.\n")
        if(length(w)!=length(semiweights[[listnum]])) stop("Weight vectors in `semiweights` must be the same length as the data.\n")
        if(sum(semiweights[[listnum]]<0)>0) stop("Weight vectors in `semiweights` must be non-negative.\n")
      }
    }
  }
  if(user_list_bounding_covars==F){
    warning("`bounding_covars` input as a charactor vector. Note that benchmarking will be done using the *combined* strength of the variables in the `bounding_covars` vector, and *not* each variable individually. To benchmark on each variable individually, `bounding_covars` should be input as a list, with each variable being its own element in the list.\n")
  }
  if(!is.numeric(w)){
    if(!is.null(semiweights)) stop("`semiweights` should be `NULL` if `w` is not a vector of prespecified weights.\n")
  }





  ##### SET-UP
  ### Define a string for the CI_method
  if(inference==TRUE){
    CI_method <- "fixed_weights"
  }
  if(inference=="reestimate"){
    CI_method <- "reestimate_weights"
  }
  if(inference==FALSE){
    CI_method <- "none"
  }

  ### Set up data
  # Only keep the variables we need
  df <- df[, colnames(df) %in% c(treatment, outcome, covars)]

  ### Get weights and semi-weights normalized/in the proper format if necessary
  # If manual or preset function is given
  if(!is.numeric(w)){
    ### Weights
    # "sanitize" weight function, which (i) grabs the correct preset weight function, if necessary, or (ii) normalizes the weights from the manual weight function, if desired
    weight_function <- wf_sanitize(w, estimand, normalize, ...)

    # apply the function to get a vector of weights
    weights <- weight_function(treatment,
                               covars,
                               df)
    ### Semiweights
    # If a list of bounding covariates is provided
    if(is.list(bounding_covars)){
      # Pull out the bounding covariates from the covariates list -- the rest will inform the semiweights
      semiweight_covars <-
        lapply(
          bounding_covars,
          function(x){setdiff(covars, x)}
        )

      # Get the semiweights -- here the semiweights will be a list of weights
      semiweights <-
        lapply(
          semiweight_covars,
          weight_function,
          treatment = treatment,
          df = df
        )
    }
    # If single vector of bounding covariates is provided
    if(is.character(bounding_covars)){
      # Pull out the bounding covariates from the covariates list -- the rest will inform the semiweights
      semiweight_covars <- setdiff(covars, bounding_covars)

      # Get the semiweights -- here the semiweights will be a vector of weights
      if(length(semiweight_covars)!=0){ # if there are still covariates left
        semiweights <- weight_function(treatment, semiweight_covars, df)
      }
      if(length(semiweight_covars)==0){ # if there are no covariates left
        semiweights <- rep(1, nrow(df))
      }
    }
  }

  # If a numeric vector of weights is given
  if(is.numeric(w)){
    # Unlike above, we don't need a weight function. Define it as NA for consistency
    weight_function <- NA

    # Just redefine the w as weights, but normalized to sum to n
    weights <- (w / sum(w)) * nrow(df)

    # Normalize the weights and semiweights if needed
    if(normalize==TRUE){
      # normalize weights
      weights <- w_normalize(w, treatment, df)

      # normalize semiweights
      if(is.list(bounding_covars)){ # if bounding_covars is a list
        semiweights <-
          lapply(
            semiweights,
            w_normalize,
            treatment = treatment,
            df = df
          )
      }
      if(is.character(bounding_covars)){ # if bounding_covars is just one vector of benchmarking variables
        semiweights <- w_normalize(semiweights, treatment, df)
      }
    }

    # If don't want normalization, but semiweights still need to sum to n
    if(normalize==FALSE){
      if(is.list(bounding_covars)){ # if bounding_covars is a list
        for(index in 1:length(bounding_covars)){
          semiweights[[index]] <- (semiweights[[index]] / sum(semiweights[[index]])) * nrow(df)
        }
      }
      if(is.character(bounding_covars)){ # if bounding_covars is just one vector of benchmarking variables
        semiweights <- (semiweights / sum(semiweights)) * nrow(df)
      }
    }
  }








  ##### (1)
  #####
  ##### compute partials for each strength and covariate bound
  ##### taking in semiweights list and weights vector
  ##### df is filtered to relevant columns only

  ### Get bounds data frame
  cov_bound_partials <- mapply(
    function(bounding_cov, w2){ # defining a function that will compute the bounds based on (i) vector of bounding covariate names, and (ii) semiweights vector
      covariate_bound(
        bounding_cov, w2,
        treatment, outcome, weights,
        kd, ky, df)
    },
    bounding_covars,
    semiweights,
    SIMPLIFY = F) %>% do.call("rbind", .) # combining everything into a dataset that has all necessary information on covariate bounds

  ### build WLS model to check for significance of weighted treatment effect
  model <- do.call("lm",
                   list(as.formula(paste(outcome, "~.", sep = "")),
                        data = df,
                        weights = weights))

  ### pass model into `wls_bias_stats` -- calculates relavant statistics needed for bias, CI, and stat sig calculation.
  bias_stats <- wls_bias_stats(treatment, outcome, df, model = model, weights)
  wls_sign <- sign(bias_stats[1]) # will be needed later

  ### compute adjusted estimates based on covariate bounds
  adjusted_estimates <- wls_sign * (abs(bias_stats['wls_estimate']) - wls_bias_standardized(cov_bound_partials[,'bound.y'],
                cov_bound_partials[,'bound.d']
                ) * bias_stats['sd_ratio'])

  ###  Inference SE, CI, and Stat sig calculations
  # If you actually wanted inference
  if(CI_method %in% c("reestimate_weights", "fixed_weights")){
    # Run bootstrap -- calculates bootstrapped wls_estimates and sd ratios
    if(CI_method == "reestimate_weights"){
      boot_stats <- wls_bias_stats_boot(weight_function, B, strata = cluster, par = par, ncpus = ncpus,
                                        treatment = treatment, outcome = outcome, df= df,
                                        covars = covars)
    }
    if(CI_method == "fixed_weights"){
      boot_stats <- wls_bias_stats_boot_fixed(weights, B, strata = cluster, par = par, ncpus = ncpus,
                                        treatment = treatment, outcome = outcome, df= df,
                                        covars = covars)
    }

    # calculate SE from the bootstrap
    wls_se_est <- sd(boot_stats[,1])

    # Make a function that will calculate a CI. Creating a function so we can "alter" the output later based on `alpha` chosen
    calc_adjusted_CIs <- function(alpha){
      # If there is more than one benchmarking bound we need
      if(nrow(cov_bound_partials)>1){
        temp <- boot_stats %>%
            apply(1,
              function(r){
                wls_sign * (abs(r[1]) - wls_bias_standardized(cov_bound_partials[,'bound.y'],
                                     cov_bound_partials[,'bound.d']) * r[2])
              }) %>%
            apply(1, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm=TRUE) %>%
            t()
        return(temp)
      }
      # If there is only one benchmarking bound we need
      if(nrow(cov_bound_partials)==1){
        temp <- boot_stats %>%
          apply(1,
                function(r){
                  wls_sign * (abs(r[1]) - wls_bias_standardized(cov_bound_partials[,'bound.y'],
                                                        cov_bound_partials[,'bound.d']) * r[2])
                }) %>%
          as.numeric()
        temp <- quantile(temp, probs = c(alpha / 2, 1 - alpha / 2), na.rm=TRUE)
        temp <- t(as.matrix(temp))
        return(temp)
      }
    }

    # Apply the function
    adjusted_CIs <- calc_adjusted_CIs(alpha)
  }

  # if you don't want inference. But will still need to make sure named output is in the right format/dimension to compile
  if(CI_method == "none"){
    # SE should be NA
    wls_se_est <- NA

    # CI function should just output NA upper and lower bounds (cols=2), the same number of times as we have bounds (rows=same # as adjusted_estimates)
    calc_adjusted_CIs <- function(alpha){
      return(matrix(NA, nrow=length(adjusted_estimates), ncol=2))
    }

    # Actual CIs should just print out the same thing as the calc_adjusted_CIs function would give you.
    adjusted_CIs <- matrix(NA, nrow=length(adjusted_estimates), ncol=2)
  }

  ### Compile the Benchmarking results information. We will compile straight WLS results later.
  bound_covar_stats <- cbind(adjusted_CIs,
                             adjusted_estimates,
                             cov_bound_partials)

  colnames(bound_covar_stats) <- c("Lower",
                                   "Upper",
                                   "Adjusted_Est",
                                   "wR2.DZ.X",
                                   "wR2.YZ.DX",
                                   "Bounding_Var",
                                   "Strength")

  rownames(bound_covar_stats) <- row.names(cov_bound_partials)







  ##### (2) Compute Robustness Values
  #####
  ### RV_{q=1}
  # Get necessary partialed out vectors
  wYperpX <- lm(as.formula(
                  paste(outcome, "~ . -", treatment, sep = "") # note dataframe has been thinned out
                  ),
                data = df,
                weights=weights)$residuals

  wYperpXD <- lm(as.formula(
                   paste(outcome, "~ .",  sep = "") # note that the dataframe has been thinned out
                  ),
                data = df,
                weights=weights)$residuals

  # Calculate necessary weighted variances -- NOTE THAT THIS REQUIRES THE WEIGHTS TO SUM TO N
  wvar.wYperpX <- mean(weights * wYperpX^2)
  wvar.wYperpXD <- mean(weights * wYperpXD^2)

  # Calculate required R2 for extreme scenario
  wR2.YD.X <- (wvar.wYperpX - wvar.wYperpXD) / wvar.wYperpX

  # Make a robustness value function that calculates the required RV (needed to alter q if desired), and use it
  calc_RV_q <- function(q){
    omega_q <- q * sqrt(abs(wR2.YD.X)) / sqrt(1 - wR2.YD.X)
    (sqrt(omega_q^4 + 4*omega_q^2) - omega_q^2)/2
  }
  RV_q <- calc_RV_q(q=1)


  ### RV_alpha
  # If we actually want inference
  if(CI_method %in% c("reestimate_weights", "fixed_weights")){
    # Make a function -- we'll need it in order to alter alpha later.
    calc_RV_alpha <- function(alpha, grid_resolution = 1000){

        # First check significance
        CI <- quantile(boot_stats[,1], probs = c(alpha/2, 1 - alpha/2))
        sig_bool <- CI[1] > 0 | CI[2]<0 #true if significant, false otherwise

        # If significant, grab RV_alpha
        if(sig_bool==TRUE){
          grid_strengths <- seq(0, 1, length.out = grid_resolution)
          CI_grid <- boot_stats %>%
            apply(1,
                  function(r){
                    wls_sign * (abs(r[1]) - wls_bias_standardized(grid_strengths,
                                         grid_strengths) * r[2])
                  }) %>%
            apply(1, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm=TRUE) %>%
            t()

          contains_zero <- apply(CI_grid, 1, function(r){r[1] < 0 & r[2] > 0})
          out <- grid_strengths[which.max(contains_zero)]
          return(out)
        }

        # If not significant, return 0
        if(sig_bool==FALSE){
          return(0)
        }
      }

    # Apply it
    RV_alpha <- calc_RV_alpha(alpha)
  }

  # If we don't want inference. But still need to output correctly named objects so code works below.
  if(CI_method=="none"){
    calc_RV_alpha <- function(alpha){
      return(NA)
    }
    RV_alpha <- NA
  }

  ### Compile robustness value results
  robustness_values <- c(wR2.YD.X, RV_q, RV_alpha)
  names(robustness_values) <- c("wR2.YD.X", "RV_{q=1}", paste("RV_{alpha=", alpha, "}", sep=""))








  ###### (3) Compile Straight WLS results
  ### Make CI function
  # If we actually want inference. Making a function because the alter function needs it
  if(CI_method %in% c("reestimate_weights","fixed_weights")){
    calc_wls_CI <- function(alpha){
      return(quantile(boot_stats[,1], probs = c(alpha/2, 1 - alpha/2)))
    }
  }
  # If we don't want inference, but still require a function that returns something of the correct dimension
  if(CI_method=="none"){
    calc_wls_CI <- function(alpha){
      return(c(NA, NA))
    }
  }

  ### Calculate CI and Compile
  wls_data <- c(bias_stats[1], wls_se_est,  calc_wls_CI(alpha))
  names(wls_data) <- c("WLS_estimate", "Est_SE", "Lower", "Upper")










  ##### (4) Output
  ### Compile in a list
  sensewls_output <- list("wls_data" = wls_data,
       "RVs" = robustness_values,
       "Covariate_Bound_Estimates" = bound_covar_stats,
       "Inference_Method" = CI_method,
       "alpha" = alpha,
       "weights" = weights,
       "semiweights" = semiweights,
       "calc_wls_CI" = calc_wls_CI,
       "calc_RV_q" = calc_RV_q,
       "calc_RV_alpha" = calc_RV_alpha,
       "calc_adjusted_CIs" = calc_adjusted_CIs,
       'plot_data' = bias_stats)

  ### Assign class
  class(sensewls_output) <- "sensewls"

  ### Return
  return(sensewls_output)
}




