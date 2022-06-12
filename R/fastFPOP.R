library(tidyverse)
library(ggplot2)
library(cowplot)

#' @title fastFPOP
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm ( fast version with Gaussian loss)
#'
#' @param data is a vector of data points (a univariate time series)
#' @param cost is a number equal to the global cost
#' @param beta is a value of penalty (a non-negative real number)
#' @param affiche if it is TRUE Table V and the plots of functions q^{i}_{t+1} are displayed on the screen at each iteration
#' @param type  simple or mid type of version
#'
#' @return a vector of changepoints, global cost
#'
#' @examples
#' data <- data_gen(n = 10, chpts = 5, means = c(0,5), noise = 1, type = 'gauss')
#' fastFPOP(data, 2*log(10), affiche = FALSE, type = 'simple')
#' fastFPOP(data, 2*log(10), affiche = FALSE, type = 'mid')

fastFPOP <- function(data, beta, affiche = FALSE, type = 'simple')
{
  if (type == 'simple')
  {
    return (fastFPOPsimpleVersion(data, beta,  affiche))
  }
  else if (type == 'mid')
  {
    return (fastFPOPmidVersion(data, beta,  affiche))
  }
  else ({stop('type should be simple or mid')})
}
