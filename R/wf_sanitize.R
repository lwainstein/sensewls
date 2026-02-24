### Helper functions
# Effective sample size
ESS <- function(w){
  sum(w)^2 / sum(w^2)
}


#' Sanitize preset or user-supplied weight functions
#' @description Grab appropriate preset weight function, or edit weight functions to (i) normalize weights if desired, or (ii) output out uniform weights if no covariates are provided.
#'
#'
#' @param w A weights object, containing either
#' (1) an appropriate weighting function, or
#' (2) a string for one of the preset weighting functions.
#' @param estimand A string indicating the desired estimand. Only required if using preset weights. Can be `"ATT"`, `"ATC"`, or `"ATE"`. Default is `NULL`.
#' @param normalize A single logical, whether or not to normalize the weights when a custom weight function is provided for `w`.
#' Weights are normalized for all the preset weighting functions.
#' @param ... Additional arguments passed to preset weighting functions.
#' Currently supports providing `caliper` and `ratio` to the [MatchIt::matchit()] weighting function
#' used by setting `w = 'matchit'`.
#'

#'
#'
#' @returns
#' Appropriate weight function that (if desired) normalizes weights.
#'
#' @import dplyr
#'
#'
#' @export


# Normalize preset or custom weight functions
wf_sanitize <- function(w, estimand = NULL, normalize, ...){
  weight_function <- if(typeof(w) == "character"){
    wf_string(w, estimand, ...)
  } else {
    w
  }

  # constructed sanitized weighting function, if needed
  if( (normalize==TRUE & is.function(w)) | is.character(w) ){
    sanitized_weight_function <- function(treatment, covars, df){
      if(length(covars)!=0){
        weights <- weight_function(treatment, covars, df)
        d <- df[, treatment]
        ess0 <- ESS(weights[d==0])
        ess1 <- ESS(weights[d==1])
        w.new <- rep(0, length(weights))
        w.new[d==0] <-(weights[d==0] / sum(weights[d==0])) * (ess0 / (ess0 + ess1))
        w.new[d==1] <-(weights[d==1] / sum(weights[d==1])) * (ess1 / (ess0 + ess1))
        w.new <- length(d) * w.new
        return(w.new)
      }
      if(length(covars)==0){
        w.new <- rep(1, nrow(df))
        return(w.new)
      }
    }
  }
  if( normalize==FALSE & is.function(w) ){
    sanitized_weight_function <- function(treatment, covars, df){
      if(length(covars)!=0){
        weights <- weight_function(treatment, covars, df)
        weights <- (weights / sum(weights)) * nrow(df)
        return(weights)
      }
      if(length(covars)==0){
        w.new <- rep(1, nrow(df))
        return(w.new)
      }
    }
  }
  return(sanitized_weight_function)
}
