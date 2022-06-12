#GAUSS
#' @title data_gen_gauss
#' @description Generation of gaussian data with a given values of means and changepoints
#' @param n number of data points
#' @param chpts a vector of increasing change-point indices. Last index is always less than 'n'.
#' By default, 'chpts = NULL' (the data without change-points).
#' @param means vector of successive means for data, by default 'means = 0'.
#' @param noise standard deviation of an additional normal noise, by default 'noise = 1'.
#' @return a vector of size n generated data with a given values of means by the segmentation.
#' @export
#' @examples
#' data <- data_gen_gauss(n = 100, chpts = 50, means = c(1,5), noise = 1)

data_gen <- function(n, chpts = NULL, means = 0, noise = 1, type  = 'gauss'){
  #---stop---#
  if (!is.null(chpts) && n <= chpts[length(chpts)]){stop('last element of changepoints is always less than n')}
  if(!is.null(chpts) && !is.numeric(chpts)){stop('changepoints are not all numeric')}
  if(is.unsorted(chpts)){stop('changepoints should be an increasing vector')}
  if(!is.numeric(means)){stop('means are not all numeric')}
  if ((length(chpts)+1) !=  length(means)){stop('The length of the means is always equal to the number of changepoints plus one')}
  if(!is.double(noise)){stop('noise is not a double')}
  if (type !='gauss') {stop('type should be gauss')}
  if(noise < 0){stop('noise must be non-negative')}
  #---function---#
  res <- vector(mode = "numeric", length = n)
  InttT<- diff(c(0,chpts,n))
  res <- rep(means, InttT) + rnorm(n, 0, noise)
  return(res)
}

#' @title q_abp1_theta
#'
#' @description Value of q^i_{t+1}(theta)-function  in theta( using cumsum of data)
#'
#' @param a  is a value of i
#' @param b is a value of t
#' @param theta is a function argument
#' @param csd_a is a value of cumsum value data(1:i)
#' @param csd2_a is a value of cumsum value data^2(1:i)
#' @param csd_b is a value of cumsum value data(1:t+1)
#' @param csd2_b is a value of cumsum value data^2(1:t+1)
#' @param beta is a value of penalty
#' @param m_am1 is the value of m_{i-1}
#'
#' @return a value of q^i_{t+1}-function  in theta
#' @export
#'
#' @examples
q_abp1_theta <- function(a, b, theta, csd_a, csd2_a, csd_b, csd2_b, beta, m_am1){
  res <- (b+2-a)*(theta-(csd_b-csd_a)/(b +2-a))^2+csd2_b-csd2_a-(csd_b-csd_a)^2/(b+2-a)+m_am1+beta
  return(res)
}


#' @title qin
#'
#' @description Value of q^i_{n}-function  in theta
#'
#' @param theta is a function argument
#' @param data is a vector of data points (a univariate time series)
#' @param min is the vector of the min values
#' @param beta is a value of penalty
#' @param i  is the left border value
#' @param n is the right border value
#'
#' @return a value of q^i_{n}-function in theta
#'
#' @examples
#'

qin <- function(theta, data, min, beta, i, n){
  res <- min[i] + beta + (n-i+1)*theta^2 - 2*sum(data[i:n])*theta + sum(data[i:n]^2)
  return(res)
}

#' @title plot_V in the time t+1
#'
#' @description Plots of the cost functions contained in the table V
#'
#' @param data is a vector of data points (a univariate time series)
#' @param V is a table of table of triplets and quadruplets
#' @param mi  is the vector of the min values
#' @param beta is a value of penalty
#' @param t is the value of t
#'
#' @return the plots of the cost functions contained in the table V
#' @export
#'
#' @examples
#'

plot_V <- function(data, V, mi, beta, t)
{
  data_q_val <- data.frame(data_val = seq(from = min(data, V[,4]), to = max(data, V[, 4]), by = 0.1))
  q_val <- rep(0, nrow(data_q_val))
  for ( r in 1 : nrow(V))
  {
    if(V[r, 1] == V[r, 2])
    { #q^lbl1_{t+1}(data_val[i])
      for (i in 1 : nrow(data_q_val))
      {
        q_val[i] = qin(data_q_val[i,1], data, mi, beta, V[r, 1], t+1)
      }
      data_q_val <- data.frame(data_q_val, q_val)
    }
  }
  plot_v <- ggplot()

  minpbeta <- data.frame(y = data_q_val[, 1], cost = rep((min(V[1, 3]) + beta), length(data_q_val[, 1])))
  plot_v <- plot_v + geom_line(data = minpbeta, aes(x = y, y = cost), col = 'red', size =0.5, linetype = 'dotted')

  dat <- vector("list", (ncol(data_q_val) - 1))
  for (cl in 1 : (ncol(data_q_val) - 1))
  {
    dat[[cl]] <- data.frame(y = data_q_val[, 1], cost = data_q_val[, cl+1])
    plot_v <- plot_v + geom_line(data = dat[[cl]], aes(x = y, y = cost), col = cl)
  }
  dat[[cl+1]] <- data.frame(y = V[, 4], cost = V[, 3])
  plot_v <- plot_v + geom_point(data = dat[[cl + 1]], aes(x = y , y = cost), shape = 4, size = 2, color = 'black')
  return (plot_v)
}

