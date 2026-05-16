
#' Inverse Propensity Score Weighting Function
#'
#'
#' @description Create Inverse Propensity Score weighting function
#'
#'
#'
#' @param estimand Character for either ATT, ATC, or ATE
#'
#'
#'
#'
#' @returns Appropriate Inverse Propensity Score weight function.
#'
#' @import dplyr
#'
#'
#' @export

wf_prop_score <- function(estimand){
  function(treatment, covars, df){
    psmodel <- glm(as.formula(paste0("`", treatment,
                                    "` ~ ",
                                    paste(paste0("`", covars, "`"), collapse = "+"))),
                   family = binomial(link="logit"),
                   data = df)

    d <- df[, treatment]
    weights <- rep(1, nrow(df))
    switch(estimand,
           ATE = {
             weights[d==0] <- 1 / (1 - psmodel$fitted.values[d==0])
             weights[d==1] <- 1 / psmodel$fitted.values[d==1]
           },
           ATC = {
             weights[d==1] <- (1 - psmodel$fitted.values[d==1]) /
               psmodel$fitted.values[d==1]
           },
           ATT = {
             weights[d==0] <- (psmodel$fitted.values[d==0]) /
               (1 - psmodel$fitted.values[d==0])
           })
    return(weights)
  }
}

