library(tidyverse)
library(ggplot2)
library(cowplot)

#' @title fastFPOPmidVersion
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm ( fast version with Gaussian loss)
#'
#' @param data is a vector of data points (a univariate time series)
#' @param cost is a number equal to the global cost
#' @param beta is a value of penalty (a non-negative real number)
#' @param affiche if it is TRUE Table v and the plots of functions q^{i}_{t+1} are displayed on the screen at each iteration
#'
#' @return a vector of changepoints, global cost
#'
#' @examples
#' library(fastFPOP)
#' data <- data_gen(n = 10, chpts = 5, means = c(0,5), noise = 1)
#' fastFPOPmidVersion(data, 2*log(10),  affiche = TRUE)
#' fastFPOPmidVersion(data, 2*log(100),  affiche = FALSE)

fastFPOPmidVersion <- function(data, beta,  affiche = FALSE) {

  # DATA PLOT
  if (affiche)
  {
    List_Plot = list()
    data_frame <- data.frame(time = c(1 : length(data)), data)
    List_Plot[[1]] <- ggplot(data_frame) + geom_point(mapping = aes(x = time, y = data), size = 2, color = 'blue')
    print(List_Plot[[1]])
  }
  # 'n' is a data length,'tau'  is a vector of best last changepoints and 'mi' is vector of optimal cost
  n <- length(data)
  tau <- c(1, rep(0, (n-1))) ## tau <- rep(1,n)
  mi <- c(-beta, rep(0, n)) # m = (-beta, 0, ...)

  # 'csData' matrix of cumulative sums data  and data^2
  csData <-matrix(c(0,cumsum(data), 0, cumsum(data^2)), ncol = 2, byrow = FALSE)
  # matrix of minimum and intersection quadruplets
  v <- matrix(nrow = 1, ncol = 4, data = c(1, 1, 0, data[1]))
  colnames(v) <- c("lbl1", "lbl2", "value", "position")
  # loop1
  for (t in 1 : (n-1))
  {
    #"br_t" is a matrix of the branch indicators (l,r) (by default,  is always "TRUE")
    br_t <- matrix(ncol = 2, nrow = (t+1), data = TRUE)
    colnames(br_t) = c("left", "right")
    #add new "triplet", and update all the points
    mpbeta <- v[1, 3] + beta
    v <- rbind(v, c(t + 1, t + 1, mpbeta, data[t + 1]))
    for(i in 1: (nrow(v) - 1) )
    {
      v[i, 4] <- ifelse(v[i, 1] == v[i, 2], (csData[t + 2, 1] - csData[v[i, 1], 1]) / (t + 2 - v[i, 1]), v[i, 4])
      v[i, 3] <- ifelse(v[i, 1] == v[i, 2], csData[t+2,2]-csData[v[i,1],2]-(csData[t + 2, 1] - csData[v[i, 1], 1])^2 / (t + 2 - v[i, 1]) + mi[v[i, 1]] + beta, v[i, 3] + (data[t + 1] - v[i, 4])^2)
    }
    # Sort and division
    fl_level <- TRUE#(v[index,3] < (v[1, 3]+beta))
    v <- v[order(v[, 3]), ]
    id <- length(v[(v[, 3] < mpbeta), 3])
    id2 <- length(v[(v[, 3] < (v[1, 3] + beta)), 3])

    if (id2 <= id)
    {
      v <-  v[-c((id2 + 1) : nrow(v)), ]
    }
    else
    {
      #SET W, INTER, UPDATE BR_T,
      w <- NULL
      inter <- matrix(ncol = 4, nrow =0)
      intSet <- matrix(ncol = 4, nrow =0)

      if (id  == 0)
      {
        w <- c(w, (t + 1))
        index <- 2
      }
      else
      {
        for (i in 1 : id)
        {
          if (v[i, 1] == v[i, 2])
          {
            w <- c(w, v[i, 1])
          }
          else
          {
            if(v[i, 1] < v[i, 2])
            {
              br_t[v[i, 1], 1] <- FALSE
            }
            else{
              br_t[v[i, 2], 2] <- FALSE
            }
          }
        }
        index <- id + 1

        #ADD INTERSECTIONS INTO INTER
        for (j in w)
        {
          if ((br_t[j, 1]) | (br_t[j, 2]))
          {
            mlbltp1 <- (t + 2 - j) * ((csData[t + 2, 2] - csData[j, 2]) / (t + 2 - j) - ((csData[t + 2, 1] - csData[j, 1]) / (t + 2 - j))^2) + mi[j] + beta
            R <- (mi[t + 1] - mi[j]) / (t + 1 - j) - ((csData[t + 1, 2] - csData[j, 2]) / (t + 1 - j) - ((csData[t + 1, 1] - csData[j, 1]) / (t + 1 - j))^2)
            if (R > 0)
            {
              R <- sqrt(R)
              delta <- ((csData[t + 1, 1] - csData[j, 1]) / (t + 1 - j) - ((csData[t + 2, 1] - csData[j, 1]) / (t + 2 - j)))
              if (br_t[j, 1])
              {
                root <- (t + 2 - j) * (delta - R)^2 + mlbltp1
                if (root < (v[1, 3] + beta))
                {
                  inter <- rbind(inter, c(j, t + 1, root, (csData[t + 1, 1] - csData[j, 1]) / (t + 1 - j) - R))
                  br_t[j, 1] <- FALSE
                }
              }
              if (br_t[j, 2])
              {
                root <- (t + 2 - j) * (delta + R)^2 + mlbltp1
                if (root < (v[1, 3] + beta))
                {
                  inter <- rbind(inter, c(t + 1, j, root, (csData[t + 1, 1] - csData[j, 1])/(t + 1 - j) + R))
                  br_t[j, 2] <-FALSE
                }
              }
            }
          }
        }
      }
      ## q-PLOT BEFORE PRUNING
      if (affiche)
      {
        List_Plot[[2]] <- plot_V (data, v, mi, beta, t)
        print('V before :')
        print(v)
      }
      fl_level <- TRUE#(v[index,3] < (v[1, 3]+beta))
      while(fl_level)
      {
        ##IF MINIMUM POINT ("lbl1" = lbl2"),ADD NEW INTERSECTIONS INTO INTER IF NOT T+1
        if ((v[index, 1] == v[index, 2]) & (v[index, 1] != (t + 1)))
        {
          j <- v[index, 1]
          if ((br_t[j, 1]) | (br_t[j, 2]))
          {
            mlbltp1 <- (t + 2 - j) * ((csData[t + 2, 2] - csData[j, 2]) / (t + 2 - j) - ((csData[t + 2, 1] - csData[j, 1]) / (t + 2 - j))^2) + mi[j] + beta
            R <- (mi[t + 1] - mi[j]) / (t + 1 - j) - ((csData[t + 1, 2] - csData[j, 2]) / (t + 1 - j) - ((csData[t + 1, 1] - csData[j, 1]) / (t + 1 -j))^2)
            if (R > 0)
            {
              R <- sqrt(R)
              delta <- ((csData[t + 1, 1] - csData[j, 1]) / (t + 1 - j) - ((csData[t + 2, 1] - csData[j, 1]) / (t + 2 - j)))
              if (br_t[j, 1])
              {
                root <- (t + 2 - j) * (delta - R)^2 + mlbltp1
                if (root < (v[1, 3] + beta))
                {
                  inter <- rbind(inter, c(j, t + 1, root, (csData[t + 1, 1] - csData[j, 1]) / (t + 1 - j) - R))
                  br_t[j, 1] <-FALSE
                }
              }
              if (br_t[j, 2])
              {
                root <- (t + 2 - j) * (delta + R)^2 + mlbltp1
                if (root < (v[1, 3] + beta))
                {
                  inter <- rbind(inter, c(t + 1, j, root, (csData[t + 1, 1] - csData[j, 1])/(t + 1 - j) + R))
                  br_t[j, 2] <-FALSE
                }
              }
            }
          }
        }
        ################################
        #### Intersection :fl_vis = TRUE - the point is visible(not covered) => add "labels" in w fl_vis = FALSE the point isn't visible(covered) => remove the element
        #### Minimum : fl_vis = TRUE - the point is visible(not covered) => add "label" in w fl_vis = FALSE the point isn't visible(covered)
        fl_vis <- TRUE
        j <- 1
        eps <- 10^(-13)
        while((j <= length(w)) & fl_vis)
        {
          qwjtp1 <- q_abp1_theta(w[j], t, v[index, 4] , csData[w[j], 1], csData[w[j], 2], csData[t + 2, 1], csData[t + 2, 2], beta, mi[w[j]])
          fl_vis <- (!((v[index, 3] - (qwjtp1 + eps)) > 0)) # point is covered >=
          j <- j + 1
        }
        if (v[index, 1] != v[index, 2])
        {
          if (fl_vis)
          {
            w <- unique(c(w, v[index, 1], v[index, 2]))
          }
          else {
            v <- v[-index, ]
            index <-index - 1
          }
        }
        if ((v[index, 1] == v[index, 2]) & (fl_vis) )
        {
          w <- unique(c(w, v[index, 1]))
        }

        index <-index + 1

        if (nrow(inter) > 0)
        {
          if (index <= nrow(v))
          {
            intSet <- matrix(inter[inter[, 3] <= v[index, 3], ], ncol = 4)
            if (nrow(intSet) > 0)
            {
              inter <- matrix(inter[-c(inter[, 3] <= v[index, 3]), ], ncol = 4)
              v <- rbind(v[1 : (index - 1), ], intSet[order(intSet[ , 3]), ], v[index : nrow(v), ])
            }
          }
          else
          {
            v <- rbind(v, inter[order(inter[ , 3]), ])
            inter <- matrix(ncol = 4, nrow =0)
          }
        }
        if (index <= nrow(v))
        {
          fl_level <- (v[index, 3] < (v[1, 3] + beta))
        }
        else
        {
          fl_level <-FALSE
        }
      }
      if (index != (nrow(v) + 1))
      {
        v <- v[-c(index : nrow(v)), ]
      }
      ###step4 :
      v <- matrix(v[(v[ , 1] == v[ , 2]) & (is.element(v[ , 1], w)) | (v[ , 1] != v[ , 2]), ], ncol = 4)
    }
    ### q-PLOT AFTER PRUNING
    if (affiche)
    {
      print('V after :')
      print(v)
      List_Plot[[3]] <- plot_V(data, v, mi, beta, t)
      print(plot_grid(List_Plot[[2]], List_Plot[[3]]))
    }
    mi[t + 2] <- v[1, 3]
    tau[t + 1] <- v[1, 1]
  }
  ##BACKTRACKING
  ind <- tau[n]
  chpnt <- tau[n]
  while (ind > 1)
  {
    chpnt <- c(chpnt, tau[ind - 1])
    ind <- tau[ind - 1]
  }
  chpnt <- rev(chpnt)[-1] - 1
  glCost <- mi[n + 1] - length(chpnt) * beta
  return(list(changepoints = chpnt, globalCost = glCost))
}

