#' Grab preset weight function
#'
#' @description Grab appropriate preset weight function
#'
#'
#' @param w String for one of the preset weighting functions.
#' @param `estimand` A string indicating the desired estimand. Only required if using preset weights. Can be `"ATT"`, `"ATC"`, or `"ATE"`. Default is `NULL`.
#' @param ... Additional arguments passed to preset weighting functions.
#' Currently supports providing `caliper` and `ratio` to the [MatchIt::matchit()] weighting function
#' used by setting `w = 'matchit'`.
#'

#'
#'
#' @returns
#' Appropriate preset weight function.
#'
#'
#' @import dplyr
#'



#' @export



wf_string <- function(w, estimand, ...){
  switch(w,
         matchit = {wf_matchit(estimand, ...)},
         ebal = {wf_ebal(estimand, ...)},
         prop_score = {wf_prop_score(estimand)},
         stop("Unsupported"))
}
