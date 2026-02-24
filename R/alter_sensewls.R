#' Altering sensewls objects
#'
#' @description This function allows one to change the
#' significance levels (`alpha`) and `q` threshold for robustness values and confidence intervals of
#' a `sensewls` object.
#'
#' @param sensewls_obj A `sensewls` object, the output from the `sensewls()` function.
#' @param alpha Updated level of significance used for RV_alpha, and (1 - alpha)*100% confidence intervals.
#' @param q Updated `q` threshold for RV_q.
#'
#'
#' @returns
#' A sensewls object containing a variety of sensitivity statistics, with updated significance levels, confidence intervals, and
#' RV_q.
#'
#' @import dplyr
#'
#' @export

alter_sensewls <- function(sensewls_obj, alpha, q){

  ### WLS data
  wls_data <- sensewls_obj$wls_data
  wls_CI <- sensewls_obj$calc_wls_CI(alpha)
  wls_data <- c(wls_data[1:2], wls_CI)
  names(wls_data) <- c("WLS_estimate", "Est_SE", "Lower", "Upper")

  ### Robustness Values
  rv_data <- sensewls_obj$RVs
  rv_data[2] <- sensewls_obj$calc_RV_q(q)
  rv_data[3] <- sensewls_obj$calc_RV_alpha(alpha)
  names(rv_data) <- c("wR2.YD.X",
                      paste("RV_{q=", q, "}", sep = ""),
                      paste("RV_{alpha=", alpha, "}", sep = ""))

  ### Adjusted CIs
  CIs_matrix <- sensewls_obj$calc_adjusted_CIs(alpha)
  cbes <- sensewls_obj$Covariate_Bound_Estimates
  cbes[,'Lower'] <- CIs_matrix[,1]
  cbes[,'Upper'] <- CIs_matrix[,2]

  sensewls_obj$wls_data <- wls_data
  sensewls_obj$RVs <- rv_data
  sensewls_obj$Covariate_Bound_Estimates <- cbes
  sensewls_obj$alpha <- alpha
  sensewls_obj
}





