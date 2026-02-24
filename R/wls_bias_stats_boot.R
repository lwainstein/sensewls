

#' Bootstrap WLS bias and weighted SD ratio
#'
#' @description Bootstrap WLS bias and weighted SD ratio
#'
#' @param weight_function A weighting function, already sanitized with `wf_sanitize`
#' @param B Number of bootstrap samples
#' @param strata String denoting stratifying variable for cluster bootstrap
#' defaulted to `NULL` which yields nonparametric resampling for bootstrap
#' @param par Boolean for whether bootstrapping should be parallelized through
#' the `parallel` library. Defaulted to `FALSE`.
#' @param ncpus Number of cpus used if `par` is set to `TRUE`. Defaulted to `NULL` and
#' must be supplied if `par` is set to true.
#' @param treatment A string indicating the column name of the treatment variable.
#' @param outcome A string indicating the column name of the outcome variable.
#' @param covars A character vector of column names of the relevant covariates in the dataframe.
#' @param df A dataframe containing `treatment`, `outcome`, and `covars`.
#'
#'
#' @import parallel
#' @import dplyr
#'
#'
#' @returns A Bx2 matrix with each row corresponding to the
#' wls estimate and weighted sd ratio of a bootstrap sample
#' @export


wls_bias_stats_boot <- function(weight_function, B, strata = NULL, par = F, ncpus = NULL,
                    treatment, outcome, df, covars){

  # Grab n-size for data
  n <- nrow(df)

  # Construct resampling function that inputs df and returns bootstrapped dataframe, depending on preset form of bootstrapping
  if(!is.null(strata)){
    levels <- unique(df[, strata])
    resample_mech <- function(df){
      boot_groups <- sample(levels, size = length(levels), replace = T)
      lapply(boot_groups,
             function(g){
               df[df[, strata] == g,]
            }) %>% do.call('rbind', .)
    }
  }else{
    resample_mech <- function(df){
      indicies_boot <- sample(seq_len(n),n,replace = T)
      df[indicies_boot, ]
    }
  }

   # Construct function to be ran multiple times
   # when evaluated, bootstraps df, calculated bootstrap weights, and calls `wls_bias_stats()`
   # returns matrix of bootstrapped adjusted wls estimated in matrix
   fn <- function(i){
     df_boot <- resample_mech(df)

     weights_boot <- weight_function(treatment, covars, df_boot)

     wls_bias_stats(outcome = outcome,
              treatment = treatment,
              df = df_boot,
              weights = weights_boot) # returns matrix of all values
   }

   # Bootstrap resampled wls_bias outputs in list object
   bootstrap_data <- if(par){
       parallel::mclapply(1:B, fn, mc.cores = ncpus)
   } else {
     lapply(1:B, fn)
   }

   # final output is Bx2 matrix of data
   return(do.call("rbind", bootstrap_data))
}





