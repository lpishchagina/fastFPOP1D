library(tidyverse)
library(ggplot2)
library(cowplot)

#' @title fastFPOPsimpleVersion
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm ( fast version with Gaussian loss)
#'
#' @param data is a vector of data points (a univariate time series)
#' @param cost is a number equal to the global cost
#' @param beta is a value of penalty (a non-negative real number)
#' @param affiche if it is TRUE Table V and the plots of functions q^{i}_{t+1} are displayed on the screen at each iteration
#'
#' @return a vector of changepoints, global cost
#'
#' @examples
#' data <- data_gen(n = 10, chpts = 5, means = c(0,5), noise = 1, type = 'gauss')
#' fastFPOPsimpleVersion(data, 2*log(10),  affiche = TRUE)

fastFPOPsimpleVersion <- function(data, beta, affiche = FALSE) {
  #DATA PLOT
  if (affiche)
  {
    List_Plot = list()
    data_frame <- data.frame(time = c(1:length(data)), data)
    List_Plot[[1]] <- ggplot(data_frame) + geom_point(mapping = aes(x = time, y = data), size = 2, color = 'blue')
    print(List_Plot[[1]])
  }

  # 'n' is a data length,'tau'  is a vector of best last changepoints and 'mi' is vector of optimal cost
  n <- length(data)
  tau <- c(1, rep(0,(n-1))) ## tau <- rep(1,n)
  mi <- c(-beta, rep(0, n)) # m = (-beta, 0, ...)

  # 'csData' table of cumulative sums data  and data^2
  # mean(data[i:t+1]) = (csData[t+1 + 1,1] - csData[i,1])/(t+1 - i + 1)
  # Var_{i:t+1} = (csData[t+1 + 1,2] - csData[i,2])/(t+1 - i + 1) - (mean(data[i:t+1]))^2
  csData <-data.frame(sum_y = c(0, cumsum(data)), sum_y2 = c(0, cumsum(data^2)))

  # table of minimum quadruplets(lbl1, lbl1, value, position) and the intersection quadruplets  (lbl1, lbl2, value, position)
  V <- data.frame(lbl1 = 1, lbl2 = 1, value = 0, position = data[1]) #add (1, 1, 0, 0)

  # loop1
  for (t in 1 : (n-1)){

    ##ADD NEW POINT : add new "triplet"
    V[nrow(V) + 1, ] <- c(t + 1, t + 1, V[1, "value"] + beta, data[t + 1]) # add new triplet element

    ##UPDATE : update all lines except the last one (it is already correct)
    for(i in 1: (nrow(V)-1))
      {
      V[i, "position"] <- ifelse(V[i, "lbl1"] == V[i, "lbl2"], (csData[t + 2, 1] - csData[V[i, "lbl1"], 1]) / (t + 2 -V[i, "lbl1"]), V[i, "position"])
      V[i, "value"] <- ifelse(V[i, "lbl1"] == V[i, "lbl2"], csData[t + 2, 2] - csData[V[i, "lbl1"], 2] - (csData[t+2, 1] - csData[V[i, "lbl1"], 1])^2 / (t + 2 - V[i, "lbl1"]) + mi[V[i, "lbl1"]] + beta, V[i, "value"] + (data[t + 1] - V[i, "position"])^2)
    }

    ##ADD NEW INTERSECTIONS
    ### lblSet is the label set for intersection.(the last label "unique(V[,1])" is always t+1, so we don't consider it)
    lblSet <- unique(V[,"lbl1"])[-length(unique(V[,"lbl1"]))]
    ###loop3: finding the intersection points for given label (t+1) with labels from"lblSet"(for all "lbl" from "lblSet" label < t+1)
    for (lbl in lblSet){
      ####calculate m^lbl_{t+1} and R^2
      mlbltp1 <- (t + 2 - lbl) * ((csData[t + 2, 2]-csData[lbl, 2]) / (t + 2 - lbl) - ((csData[t + 2, 1] - csData[lbl, 1]) / (t + 2 - lbl))^2) + mi[lbl] + beta
      R <- (mi[t + 1] - mi[lbl]) / (t + 1 - lbl) - ((csData[t + 1, 2] - csData[lbl, 2]) / (t + 1 - lbl) - ((csData[t + 1, 1] - csData[lbl, 1]) / (t + 1 - lbl))^2)
      ####add new pair intersections
      if (R > 0)
      {
        R <- sqrt(R)
        delta <- ((csData[t + 1, 1] - csData[lbl, 1]) / (t + 1 - lbl) - ((csData[t + 2, 1] - csData[lbl, 1])/(t + 2 - lbl)))
        V[nrow(V) + 1,]  <- c(lbl, t + 1, (t + 2 - lbl) * (delta - R)^2 + mlbltp1, (csData[t + 1, 1] - csData[lbl, 1]) / (t + 1 - lbl) - R)
        V[nrow(V) + 1,]  <- c(t + 1, lbl, (t + 2 - lbl)*(delta + R)^2 + mlbltp1, (csData[t + 1, 1] - csData[lbl, 1])/(t + 1 - lbl) + R)
      }
    }

    ##SORT
    V <- V %>% arrange(value)

    ## q-PLOT BEFORE PRUNING
    if (affiche)
    {
      List_Plot[[2]] <- plot_V(data, V, mi, beta, t)
      print('V before :')
      print(V)
    }

    ##PRUNING
    ### step1 + step2 + step 3 : w is the vector of the labels of the available functions in functional cost
    w <- c(V[1, "lbl1"]) #by default
    index <- 2
    while((index <= nrow(V)) & (V[index,'value'] < (V[1, 'value'] + beta)))
    {
      #### Intersection :fl = TRUE - the point is visible(not covered) => add "labels" in w fl = FALSE the point isn't visible(covered) => remove the element
      #### Minimum : fl = TRUE - the point is visible(not covered) => add "label" in w fl = FALSE the point isn't visible(covered)
      fl <- TRUE
      j <- 1
      eps <- 10^(-13)
      while((j <= length(w)) & fl)
      {
        qwjtp1 <- q_abp1_theta(w[j], t, V[index, "position"] , csData[w[j], 1], csData[w[j], 2], csData[t + 2, 1], csData[t + 2, 2], beta, mi[w[j]])
        fl <- (!((V[index, "value"] - (qwjtp1 + eps)) > 0))# point is covered >=
        j <- j + 1
      }
      if (V[index,"lbl1"] != V[index,"lbl2"])
      {
        if (fl)
        {
          w <- unique(c(w, V[index, "lbl1"], V[index, "lbl2"]))
        }
        else
        {
          V <- V[-index, ]
          index <- index - 1
        }
      }
      if ((V[index, "lbl1"] == V[index, "lbl2"]) & (fl))
      {
        w <- unique(c(w, V[index, "lbl1"]))
      }
      index <-index + 1
    }
    if (index != (nrow(V) + 1))
    {
      V <- V[-c(index:nrow(V)), ]
    }

    ###step4 :
    V <- dplyr::filter(V, ((lbl1 == lbl2) & is.element(V[, "lbl1"], w)) | (lbl1 != lbl2))

    ### q-PLOT AFTER PRUNING
    if (affiche)
    {
      print('V after :')
      print(V)
      List_Plot[[3]] <- plot_V(data, V, mi, beta, t)
      print(plot_grid(List_Plot[[2]], List_Plot[[3]]))
    }

    mi[t + 2] <- V[1, "value"]
    tau[t + 1] <- V[1, "lbl1"]
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
  glCost = mi[n + 1] - length(chpnt) * beta
  return(list(changepoints = chpnt, globalCost = glCost))
}
