
#' Propensity Score Matching Weighting Function
#'
#'
#' @description Create MatchIt weighting function
#'
#'
#'
#' @param estimand Character for either ATT or ATC
#' @param caliper Caliper for matching
#' @param ratio Number of units to match to
#'
#'
#'
#'
#' @returns Appropriate MatchIt weight function.
#'
#' @import MatchIt
#' @import dplyr
#'
#' @export



wf_matchit <- function(estimand, caliper = 0.1, ratio = 5){
  function(treatment, covars, df){
    MatchIt::matchit(paste0("`", treatment,
                  "`~",
                  paste(paste0("`", covars, "`"), collapse = " + ")) %>%
              as.formula,
            data = df,
            caliper = caliper,
            ratio = ratio,
            estimand = estimand,
            distance = "glm",
            link = "logit",
    )$weights
  }
}


