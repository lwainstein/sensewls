



#' Calculate WLS estimate and weighted SD ratio needed for bias
#'
#' @description Calculate WLS estimate and weighted SD ratio needed for bias
#'
#' @param model Weighted `lm` object. If `NULL` (the default), function will presume
#'  that a dataframe `df` is supplied instead
#' @param df Input dataframe which contains  `outcome` and  `treatment` columns along with additional covariates
#' @param treatment String name of the of treatment variable in supplied model or dataframe
#' @param outcome String indicating name of treatment variable in supplied model or dataframe
#' @param weights weights vector, weights summing to n
#'
#' @returns Vector of wls estimate and weighted SD ratio
#'
#' @import dplyr
#'
#'
#' @export

wls_bias_stats <- function(treatment, outcome, df, model = NULL, weights){

  if(is.null(model)){
   model <- do.call("lm",
                    list(as.formula(paste0("`", outcome, "`~.")),
                         data = df,
                         weights = weights)) #fit weighted OLS
  }

  wls_estimate <- coef(model)[treatment]

  sd_y_perp_xd <- sqrt(mean(weights * model$residuals^2))

  d_perp_x <- do.call("lm",
                      list(as.formula(paste0("`", treatment, "` ~ . - `", outcome, "`" )),
                           data =  df,
                           weights = weights))

  #calculate weighted standard error of residuals of D ~ X
  #using same weights from original weighted model
  sd_d_perp_x <-  sqrt(mean(weights * d_perp_x$residuals^2))

  #return vector of biases
  sd_ratio <- sd_y_perp_xd / sd_d_perp_x

  out <- c(wls_estimate, sd_ratio)
  names(out) <- c("wls_estimate", "sd_ratio")
  out
}







