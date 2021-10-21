# Written by Kimlee Chea
# Date: 2021-09-27
# STAT 128-01 F21

install.packages("tidyverse")
install.packages("lattice")
install.packages("gridExtra")

library(tidyverse)
library(lattice)
library(gridExtra)

set.seed(1095)


#'Create a randomly infected matrix.
#'
#' @param x number of rows in grid
#' @param z number of columns in grid
#' @param w probability of initial infection
#' @export
step_infection = function(x=50,z=50,w=0.1){
  if(x >=1 && z >=1 && w >= 0 && w <= 1){
    y = sample(c(1, 0), size = x * z, replace = TRUE, prob = c(w, (1-w)))
    dd <- matrix(y, nrow = x, ncol = z)
    class(dd) = c("SIR", class(dd))
    return(dd)
  }else{
    stop("x must be a positive integer; z must be a positive integer; w must be between 0 and 1")
  }
}

#'Create a matrix with previous infections removed/recovered and calculate new infection spread.
#'
#' @param x a matrix with susceptible, infected, and removed values
#' @param y probability of infecting a neighbor cell
#' @export
step_susceptible = function(x, y = 0.2){
  new_infections = x
  new_infections = step_removed(new_infections)
  infcoord <- which(new_infections == 2, arr.ind = T)
  n = nrow(x)
  p = ncol(x)

  # TOP CELL
  cc = sample(c(1, 0), size = (n)*(p-1), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n, ncol = p-1)
  infected = x[, -p] == 1
  new_infections[,-1][infected & cc] = 1

  # BOTTOM CELL
  cc = sample(c(1, 0), size = (n)*(p-1), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n, ncol = p-1)
  infected = x[,-1] == 1
  new_infections[,-p][infected & cc] = 1

  # LEFT CELL
  cc = sample(c(1, 0), size = (n-1)*(p), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n-1, ncol = p)
  infected = x[-1,] == 1
  new_infections[-n,][infected & cc] = 1

  # RIGHT CELL
  cc = sample(c(1, 0), size = (n-1)*(p), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n-1, ncol = p)
  infected = x[-n,] == 1
  new_infections[-1,][infected & cc] = 1

  # DIAGONAL TOP-RIGHT CELL
  cc = sample(c(1, 0), size = (n-1)*(p-1), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n-1, ncol = p-1)
  infected = x[-n, -p] == 1
  new_infections[-1,-1][infected & cc] = 1

  # DIAGONAL TOP-LEFT CELL
  cc = sample(c(1, 0), size = (n-1)*(p-1), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n-1, ncol = p-1)
  infected = x[-1, -p] == 1
  new_infections[-n,-1][infected & cc] = 1

  # DIAGONAL BOTTOM-RIGHT CELL
  cc = sample(c(1, 0), size = (n-1)*(p-1), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n-1, ncol = p-1)
  infected = x[-n, -1] == 1
  new_infections[-1,-p][infected & cc] = 1

  # DIAGONAL BOTTOM-LEFT CELL
  cc = sample(c(1, 0), size = (n-1)*(p-1), replace = TRUE, prob = c(y,(1-y)))
  cc = matrix(cc, nrow = n-1, ncol = p-1)
  infected = x[-1,-1] == 1
  new_infections[-n,-p][infected & cc] = 1

  # SET PREVIOUSLY INFECTED/REMOVED GREY CELLS
  for(i in 1:nrow(infcoord)){
      new_infections[infcoord[i,1],infcoord[i,2]] <- 2
  }

  return(new_infections)
}

#'Return a matrix with infected converted to removed/recovered cells.
#'
#' @param x a matrix with susceptible, infected, and removed values
#' @export
step_removed = function(x){
  x[x==1] <- 2
  return(x)
}

#'Returns true if current matrix contains infections; otherwise false.
#'
#' @param x a matrix with susceptible, infected, and removed values
#' @export
is_infectious = function(x){
  hasinfecting <- any(x==1)
  if(hasinfecting==TRUE){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#'Plot a matrix
#'
#' @param dd a matrix with susceptible, infected, and removed values
#' @export
plot_infection = function(dd){
  sir_colors = c("white","red","grey")
  lattice::levelplot( dd,
             col.regions=sir_colors, colorkey = FALSE,
             at = c(-Inf, seq(0.0, 1, by = 1), Inf))
}

#'Create an SIR simulation with given arguments.
#'
#' @param nrow number of rows in grid
#' @param ncol number of columns in grid
#' @param p probability of infecting a neighbor cell
#' @param p0 initial probability of population to be infected
#' @param plot set to TRUE to display plot of matrix
#' @param time_per_frame time interval between each infection step
#' @export
sir = function( x = step_infection(nrow, ncol, p0),
    nrow = 50, ncol = 50, p = 0.2, p0 = 0.1, plot = FALSE, time_per_frame = 0.2, ...){
  if(plot == TRUE){
    infection <- step_infection(nrow, ncol, p0)
    plot_infection(infection)
    Sys.sleep(time_per_frame)
    repeat{
      infection <- step_susceptible(infection,p)
      p1 <- plot_infection(infection)
      gridExtra::grid.arrange(p1, nrow=1)
      Sys.sleep(time_per_frame)
      if(is_infectious(infection) == FALSE){
        break
      }
    }
  }
}

sir(plot = TRUE)

