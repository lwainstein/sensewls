

#' Calculate standardized bias based on sensitivity parameters
#'
#' @description Calculate standardized bias based on sensitivity parameters

#' @param r2wyz Vector of weighted squared partial correlations between Y and Z (Conditional on D and X)
#' @param r2wdz Vector of weighted squared partial correlations between D and Z (Conditional on X)
#' @returns Vector of Z-related standardized bias factors for bias calculation
#'
#'
#' @import dplyr
#'
#' @export
#'


wls_bias_standardized <- function(r2wyz, r2wdz){
  temp <- sqrt(r2wyz) * sqrt(r2wdz) / sqrt(1 - r2wdz)

  temp[is.nan(temp)] <- NA
  temp[is.na(temp)] <- NA
  temp[is.infinite(temp)] <- NA

  return(temp)
}






