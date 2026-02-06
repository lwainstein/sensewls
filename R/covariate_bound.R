
### Helper functions
# R_{wv}^2(D ~ Xsub | Xrem)
partial_DXsub.Xrem <- function(d,y,X,Xrem, wv){
  wvDperpX <- lm(d ~ . - y, data = X, weights = wv)$residuals
  wvDperpXrem <- lm(d ~ . - y, data = Xrem, weights = wv)$residuals

  wvvar.DperpX <- mean(wv * wvDperpX^2)
  wvvar.DperpXrem <- mean(wv * wvDperpXrem^2)

  wvR2.DXsub.Xrem <- (wvvar.DperpXrem - wvvar.DperpX) / wvvar.DperpXrem
  wvR2.DXsub.Xrem
}


#' Calculating covariate bounds
#'
#' @description Calculate covariate bound for a single vector of bounding covariates
#'
#' @param bounding_cov Character vector of bounding covariates
#' @param w2 semiweights vector for associated bounding covariate. Should sum to n.
#' @param treatment Character vector indicating treatment column in `df`
#' @param `outcome` A string indicating the column name of the outcome variable in `df`
#' @param w weights vector (already calculated). Should sum to n.
#' @param kd treatment strength ratio
#' @param ky outcome strength ratio
#' @param df Source Dataframe -- prefiltered to contain only relevant covariates
#'
#'
#'
#' @returns Data frame of covariate bounds
#' @import dplyr
#'
#' @export
#'
#'
covariate_bound <- function(bounding_cov, w2, treatment, outcome, w, kd, ky, df){
  X <- df
  Xrem <- df[, !(colnames(df) %in% c(bounding_cov))]
  d <- df[, treatment]
  y <- df[, outcome]

  ### Treatment Bounds
  wDperpX <- lm(as.formula(paste(treatment, '~ . -', outcome)),
                data = X, weights = w)$residuals
  wDperpXrem <- lm(as.formula(paste(treatment, '~ . -', outcome)),
                   data = Xrem, weights = w)$residuals

  wvar.DperpX <- mean(w * wDperpX^2)
  wvar.DperpXrem <- mean(w* wDperpXrem^2)

  wR2.DXsub.Xrem <- (wvar.DperpXrem - wvar.DperpX) / wvar.DperpXrem

  w2DperpX <- lm(as.formula(paste(treatment, '~ . -', outcome)),
                 data = X, weights = w2)$residuals
  w2DperpXrem <- lm(as.formula(paste(treatment, '~ . -', outcome)),
                    data = Xrem, weights = w2)$residuals

  w2var.DperpX <- mean(w2 * w2DperpX^2)
  w2var.DperpXrem <- mean(w2 * w2DperpXrem^2)

  w2R2.DXsub.Xrem <- (w2var.DperpXrem - w2var.DperpX) / w2var.DperpXrem


  ### Outcome Bounds
  ## R_{wv}^2(Y~Xsub|D,X)
  wYperpXD <- lm(as.formula(paste(outcome, '~ .')),
                 data = X, weights=w)$residuals
  wYperpXremD <- lm(as.formula(paste(outcome, '~ .')),
                    data = Xrem, weights=w)$residuals

  wvar.wYperpXD <- mean(w * wYperpXD^2)
  wvar.wYperpXremD <- mean(w * wYperpXremD^2)
  wR2.YXsub.XremD <- (wvar.wYperpXremD - wvar.wYperpXD) / wvar.wYperpXremD

  ## Loop through the kd, ky combinations
  bound.d <- NULL
  bound.y <- NULL
  strengths <- NULL
  for(k1 in kd){
  for(k2 in ky){
    wR2.ZXsub.XremD <-  k1 * w2R2.DXsub.Xrem / (1 - k1 * w2R2.DXsub.Xrem) *
                        wR2.DXsub.Xrem / (1 - wR2.DXsub.Xrem)

    eta <- ( sqrt(k2) + sqrt(wR2.ZXsub.XremD) ) / sqrt(1-wR2.ZXsub.XremD)

    temp.d <- k1 * w2R2.DXsub.Xrem / (1 - wR2.DXsub.Xrem)
    temp.d <- ifelse(temp.d > 1, 1, temp.d)
    bound.d <- c(bound.d, temp.d)

    temp.y <- eta^2 * wR2.YXsub.XremD / (1 - wR2.YXsub.XremD)
    temp.y <- ifelse(temp.y > 1, 1, temp.y)
    bound.y <- c(bound.y, temp.y)

    strengths <- c(strengths, paste("kd=", k1, ", ", "ky=", k2, sep=""))
  }
  }

  data.frame(bound.d, bound.y, bounding_cov = paste(bounding_cov, collapse = " "), strengths)
}







