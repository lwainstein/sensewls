


#' Entropy Balancing Weighting Function
#'
#'
#' @description Create Entropy Balancing weighting function
#'
#'
#'
#' @param estimand Character for either ATT, ATC, or ATE
#'
#'
#'
#'
#' @returns Appropriate Entropy Balancing weight function.
#'
#' @import ebal
#' @import dplyr
#'
#' @export



wf_ebal <- function(estimand){
  function(treatment, covars, df){
    n <- nrow(df)
    weights <- rep(1, n)
    d <- df[,treatment]
    switch(estimand,
           ATT = {
             weights[d==0] <- ebal::ebalance(df[,treatment],
                                       as.matrix(df[,covars]))$w
           },
           ATC = {
             weights[d==1] <- ebal::ebalance(1-df[,treatment],
                                       as.matrix(df[,covars]))$w
           },
           ATE = {
             n1 <- sum(d)
             n0 <- n - n1
             X <- as.matrix(df[, covars])
             d_weights0 <- c(rep(1, n), rep(0, n0))
             X_weights0 <- rbind(X, X[d==0, ])
             output0 <- ebal::ebalance(
               Treatment=d_weights0,
               X=X_weights0
             )
             d_weights1 <- c(rep(1, n), rep(0, n1))
             X_weights1 <- rbind(X, X[d==1, ])
             output1 <- ebal::ebalance(
               Treatment=d_weights1,
               X=X_weights1
             )
             weights[d==0] <- output0$w
             weights[d==1] <- output1$w
           })
    return(weights)
  }
}
