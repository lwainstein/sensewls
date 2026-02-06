

#' Contour plot method for sensewls objects
#'
#' @description Forms a contour plot of the adjusted estimates for values of the unknown partial
#' correlations between a confounder and treatment and the confounder and outcome.
#'
#'
#' @param sensewls_obj An object of class `sensewls`, constructed from the `sensewls` function
#' @param contour_levels A vector containing values of adjusted estimates to plot the level sets of
#' @param xlim A vector providing the lower and upper bounds of the x-axis.
#' Must be a subinterval of the unit interval
#' @param ylim A vector providing the lower and upper bounds of the y-axis.
#' Must be a subinterval of the unit interval.
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
plot_sensewls <- function(sensewls_obj, contour_levels = NULL, xlim = c(0,1), ylim = c(0,1),
                     palette = 7)
  {

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

  ### Add in extra contour lines if necessary
  if(!is.null(contour_levels)){
    for (tau in contour_levels){
      plot_base <-  plot_base +
        geom_function(fun = null_function,
                      args = list(l = tau),
                      color = 'grey',
                      aes(linetype = "Additional Adjusted Estimate(s)"))
    }
  }

  ### Add in adjusted estimates
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
                   color = Strength),
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




