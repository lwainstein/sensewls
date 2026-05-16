

#' Contour plot method for sensewls objects
#'
#' @description Forms a contour plot of the adjusted estimates for values of the unknown partial
#' correlations between a confounder and treatment and the confounder and outcome.
#'
#'
#' @param sensewls_obj An object of class `sensewls`, constructed from the `sensewls` function
#' @param xlim A vector providing the lower and upper bounds of the x-axis.
#' Must be a subinterval of the unit interval
#' @param ylim A vector providing the lower and upper bounds of the y-axis.
#' Must be a subinterval of the unit interval.
#' @param sig_contour A single logical, indicating whether or not to plot the contour line for the *first*
#' adjusted estimate that is statistically insignificant (at the `sig_alpha` level). Only allowed if bootstrap inference was used when creating `sensewls_obj`. Default value is `FALSE`.
#' @param sig_alpha Number between 0 and 1. Updated level of significance for the `sig_contour` line. If `NULL`, then assumes the significance level from provided `sensewls_obj`. Default value is `NULL`.
#' @param palette An argument which is supplied to `scale_color_brewer` to
#' color the scatterplot with a seqential palette.
#'
#'
#' @returns
#' Contour plot
#'
#'
#' @import ggrepel
#' @import ggplot2
#' @import dplyr
#'
#'
#' @export
#'
#'
#'
#'
plot_sensewls <- function(sensewls_obj, xlim = c(0,1), ylim = c(0,1), sig_contour = FALSE, sig_alpha = NULL,
                     palette = 7)
  {

  ##### ERROR/WARNING MESSAGES
  ### check basic argument quality
  # sensewls_obj
  if(
    class(sensewls_obj)!="sensewls"
  ){stop("`sensewls_obj` must be a `sensewls` object.\n")}

  # xlim
  if(
    !is.numeric(xlim)
  ){stop("`xlim` must be a numeric vector of length 2, and a subinterval of the unit interval (0, 1).\n")}
  if(
    is.numeric(xlim) & length(xlim)!=2
  ){stop("`xlim` must be a numeric vector of length 2, and a subinterval of the unit interval (0, 1).\n")}
  if(
    xlim[1]<0 | xlim[1]>1 | xlim[2]<0 | xlim[2]>1 | xlim[1]>=xlim[2]
  ){stop("`xlim` must be a subinterval of the unit interval (0, 1).\n")}

  # ylim
  if(
    !is.numeric(ylim)
  ){stop("`ylim` must be a numeric vector of length 2, and a subinterval of the unit interval (0, 1).\n")}
  if(
    is.numeric(ylim) & length(ylim)!=2
  ){stop("`ylim` must be a numeric vector of length 2, and a subinterval of the unit interval (0, 1).\n")}
  if(
    ylim[1]<0 | ylim[1]>1 | ylim[2]<0 | ylim[2]>1 | ylim[1]>=ylim[2]
  ){stop("`ylim` must be a subinterval of the unit interval (0, 1).\n")}

  # sig_contour
  if(
    !is.logical(sig_contour)
  ){stop("`sig_contour` must be a single logical.\n")}
  if(
    is.logical(sig_contour) & length(sig_contour)>1
  ){stop("`sig_contour` must be a single logical.\n")}
  if(
    sig_contour==TRUE & !(sensewls_obj$Inference_Method %in% c("boot", "fix_boot"))
  ){stop("`sig_contour=TRUE` is only allowed when bootstrap inference (`'boot'` or `'fix_boot'`) was used to create `sensewls_obj`.\n")}

  # sig_alpha
  if(!is.null(sig_alpha)){
    if(
      !is.numeric(sig_alpha)
    ){stop("`sig_alpha` must be a single number, or `NULL`.\n")}
    if(
      is.numeric(sig_alpha) & length(sig_alpha)>1
    ){stop("`sig_alpha` must be a single number, or `NULL`.\n")}
    if(sig_alpha<=0 | sig_alpha >=1){stop("`sig_alpha` must be a number between 0 and 1, or `NULL`.\n")}
  }





  ##### PLOTTING
  ### Pull out necessary information from sensewls object
  bias_stats <- sensewls_obj$plot_data

  covariates_data <- sensewls_obj$Covariate_Bound_Estimates
  strength_values <- covariates_data$Strength
  covariates_data <- covariates_data %>%
                     mutate(Strength = factor(Strength,
                                              levels = rev(sort(unique(strength_values)))))
  adj_ests <- covariates_data$Adjusted_Est
  covariates_data <- covariates_data %>%
    mutate(label = paste("(", round(Adjusted_Est, 3), ")"
                         , sep = ""))

  # Only if they want the significance contour
  actually_plot <- F # whether or not this actually gets plotted, based on RV
  if(sig_contour==TRUE){
    # Get correct alpha level
    if(is.null(sig_alpha)) alpha <- sensewls_obj$alpha
    if(!is.null(sig_alpha)){
      alpha <- sig_alpha
      sensewls_obj <- alter_sensewls(sensewls_obj, alpha=sig_alpha) # need to update object
    }

    # Message
    message_text <- paste0(
      "sig_contour=TRUE: Note that this only plots the contour for the *first* adjusted estimate that is statistically insignificant (at the ",
      alpha,
      " level). It is possible that an adjusted estimate is still statistically significant beyond this contour. Users should confirm this is not the case for any points beyond the contour."
    )
    message(message_text)

    # Find the adjusted estimate that is not stat sig
    temp_name <- paste0("RV_{alpha=", alpha, "}")
    RV_alpha <- sensewls_obj$RVs[temp_name]
    sign_wls <- sign(bias_stats[1])
    wls_insig  <- sign_wls * (abs(bias_stats[1]) - wls_bias_standardized(RV_alpha, RV_alpha) * bias_stats[2])

    # Determine whether or not to plot
    if(RV_alpha==0){
      actually_plot <- F
      message_text <- paste0(
        "Unadjusted estimate is already statistically insignificant at the ",
        alpha,
        " level. Not plotting contour for first statistically insignificant estimate. Consider choosing `sig_alpha`>",
        alpha,
        " if this contour is still desired."
      )
      warning(message_text)
    }
    if(RV_alpha>0) actually_plot <- T
  }

  ### Create function that will plot contours for a fixed adjusted estimate
  null_function <- function(rd, l){
    (sqrt(1 - rd) / sqrt(rd) * (bias_stats[1] - l) / bias_stats[2])^2
  }

  ### Get plot base with null effect contours plotted
  plot_base <- ggplot() +
    theme_minimal()  +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    geom_function(fun = null_function,
                  args = list(l = 0),
                  color = 'gray30',
                  aes(linetype = "Null Effect"))

  ### Add in Significance contour?
  if(actually_plot==T){
    plot_base <-
      plot_base +
      geom_function(fun = null_function,
                  args = list(l = wls_insig),
                  color = 'gray30',
                  aes(linetype = paste0("First Insignificant Effect (at ", alpha, " level)")))
  }

  ### Add in adjusted estimates from benchmarking
  # Get correct ordering of strengths
  temp <- covariates_data[, "Strength"]
  temp <- unique(temp)
  temp <- gsub("kd=", "", temp)
  temp <- gsub("ky=", "", temp)
  temp <- strsplit(temp, split=", ")
  temp <- t(data.frame(temp))
  temp <- cbind(as.numeric(temp[, 1]), as.numeric(temp[, 2]))
  if(nrow(temp)>1) temp <- temp[order(temp[, 1], temp[, 2]), ]
  strength_levels <- paste0("kd=", temp[, 1], ", ky=", temp[, 2])

  # Plot
    plot_base +
      geom_point(data = covariates_data,
                 aes(y = wR2.YZ.DX,
                     x = wR2.DZ.X,
                     shape = Bounding_Var),
                 inherit.aes = F, size = 2.75, color = 'black') +
    geom_point(data = covariates_data,
               aes(y = wR2.YZ.DX,
                   x = wR2.DZ.X,
                   shape = Bounding_Var,
                   color = factor(Strength, strength_levels)),
               inherit.aes = F, size = 2) +
    ggrepel::geom_text_repel(data = covariates_data,
              aes(label = label,
                  y = wR2.YZ.DX,
                  x = wR2.DZ.X,
                  ), vjust = "inward",
                      hjust = 'inward',
              size = 3.4) +
    labs(title = 'Adjusted Estimate Contours',
         color = 'Strength Ratio',
         shape = "Bounding Variables",
         x = 'Partial R^2 of confounder(s) with treatment',
         y = 'Partial R^2 of confounder(s) with outcome',
         linetype = "Estimate Contour") +
    scale_color_brewer(type = "seq",
                       direction = -1,
                       palette = palette) +
    theme(legend.position="right",
          plot.title = element_text(hjust = 0.5))
}




