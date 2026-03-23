#' Altering sensewls objects
#'
#' @description This function allows one to change the
#' significance levels (`alpha`) and `q` threshold for robustness values and confidence intervals of
#' a `sensewls` object.
#'
#' @param sensewls_obj A `sensewls` object, the output from the `sensewls()` function.
#' @param alpha Single number from 0 to 1. Updated level of significance used for RV_alpha, and (1 - alpha)*100% confidence intervals. If `NULL`, then the significance level from `sensewls_obj` is left unchanged. Defaults to `NULL`.
#' @param q Single number greater than 0. Updated `q` threshold for RV_q. If `NULL`, then the `q` threshold from `sensewls_obj` is left unchanged. Defaults to `NULL`.
#'
#'
#' @returns
#' A sensewls object containing a variety of sensitivity statistics, with updated significance levels, confidence intervals, and
#' RV_q.
#'
#' @import dplyr
#'
#' @export

alter_sensewls <- function(sensewls_obj, alpha = NULL, q = NULL){

  ##### ERROR/WARNING MESSAGES
  ### check argument quality for sensewls_obj
  # sensewls_obj
  if(
    class(sensewls_obj)!="sensewls"
  ){stop("`sensewls_obj` must be a `sensewls` object.\n")}

  # alpha
  if(!is.null(alpha)){
    if(
      !is.numeric(alpha)
    ){stop("`alpha` must be a single number, or `NULL`.\n")}
    if(
      is.numeric(alpha) & length(alpha)>1
    ){stop("`alpha` must be a single number, or `NULL`.\n")}
    if(alpha<=0 | alpha >=1){stop("`alpha` must be a number between 0 and 1, or `NULL`.\n")}
  }

  # q
  if(!is.null(q)){
    if(
      !is.numeric(q)
    ){stop("`q` must be a single number, or `NULL.\n")}
    if(
      is.numeric(q) & length(q)>1
    ){stop("`q` must be a single number, or `NULL.\n")}
    if(q<=0){stop("`q` must be a positive number, or `NULL.\n")}
  }


  ### Not altering object
  if(is.null(alpha) & is.null(q)){
    warning("Neither `alpha` or `q` has been altered. Returning original `sensewls_obj`.")
  }

  ### Altering object, if desired
  if(!is.null(alpha) | !is.null(q)){

    ### WLS data
    wls_data <- sensewls_obj$wls_data
    if(!is.null(alpha)){
      wls_CI <- sensewls_obj$calc_wls_CI(alpha)
      wls_data <- c(wls_data[1:2], wls_CI)
      names(wls_data) <- c("WLS_estimate", "Est_SE", "Lower", "Upper")
    }

    ### Robustness Values
    rv_data <- sensewls_obj$RVs
    if(!is.null(q)){
      rv_data[2] <- sensewls_obj$calc_RV_q(q)
      names(rv_data)[2] <- paste("RV_{q=", q, "}", sep = "")
    }
    if(!is.null(alpha)){
      rv_data[3] <- sensewls_obj$calc_RV_alpha(alpha)
      names(rv_data)[3] <-  paste("RV_{alpha=", alpha, "}", sep = "")
    }

    ### Adjusted CIs
    cbes <- sensewls_obj$Covariate_Bound_Estimates
    if(!is.null(alpha)){
      CIs_matrix <- sensewls_obj$calc_adjusted_CIs(alpha)
      cbes[,'Lower'] <- CIs_matrix[,1]
      cbes[,'Upper'] <- CIs_matrix[,2]
    }

    ### Put everything together
    sensewls_obj$wls_data <- wls_data
    sensewls_obj$RVs <- rv_data
    sensewls_obj$Covariate_Bound_Estimates <- cbes
    sensewls_obj$alpha <- alpha
  }


  ### Return object
  return(sensewls_obj)
}





