### Helper functions
# Effective sample size
ESS <- function(w){
  sum(w)^2 / sum(w^2)
}


#' Normalize vector of weights per Wainstein and Hazlett (2025).
#'
#' @description Normalize vector of weights per Wainstein and Hazlett (2025).
#'
#'
#' @param w A weights vector.
#' @param `df` A dataframe containing `treatment`.
#' @param `treatment` A string indicating the columnname of the treatment variable.
#'

#'
#'
#' @returns
#' Normalized weight vector
#' @import dplyr
#'
#' @export



# Normalize weight vectors
w_normalize <- function(w, treatment, df){
    d <- df[, treatment]
    ess0 <- ESS(w[d==0])
    ess1 <- ESS(w[d==1])
    w.new <- rep(0, length(w))
    w.new[d==0] <-(w[d==0] / sum(w[d==0])) * (ess0 / (ess0 + ess1))
    w.new[d==1] <-(w[d==1] / sum(w[d==1])) * (ess1 / (ess0 + ess1))
    w.new <- length(d) * w.new
    return(w.new)
}

