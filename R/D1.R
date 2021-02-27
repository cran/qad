#' Calculate the D1 distance between two dependence structures
#'
#' @description Computation of the D1 distance between two checkerboard copulas A and B, corresponding to the random vectors (X1,Y1) and (X2,Y2), respectively.
#' The function \code{D1()} computes the difference between the dependence structures of two random vectors. The function \code{D1.ECBC()} computes the D1-distance between two checkerboard copulas with the same resolution.
#' The function \code{zeta1()} is defined as 3D1(A,Pi), where Pi denotes the independence copula and returns the dependence measure computed in qad.
#'
#' @param x1 a (non-empty) numeric vector of data values for the first random vector (first coordinate)
#' @param y1 a (non-empty) numeric vector of data values for the first random vector (second coordinate)
#' @param x2 a (non-empty) numeric vector of data values for the second random vector (first coordinate)
#' @param y2 a (non-empty) numeric vector of data values for the second random vector (second coordinate)
#'
#' @return \code{D1()} returns the D1 distance, introduced in (Trutschnig, 2011).
#'
#' @references
#' Trutschnig, W. (2011). On a strong metric on the space of copulas and its induced dependence measure. Journal of Mathematical Analysis and Applications. 384 (2), 690-705.
#'
#' Junker, R.R., Griessenberger, F. and Trutschnig, W. (2021). Estimating scale-invariant directed dependence of bivariate distributions. Computational Statistics and Data Analysis, 153, 107058.
#'
#' @examples
#'
#' n <- 100
#' x1 <- runif(n)
#' y1 <- x1
#' x2 <- runif(n)
#' y2 <- 1-x2
#' D1(x1,y1,x2,y2)
#'
#' n <- 1000
#' x <- runif(n, 0, 1)
#' y1 <- ifelse(x < 0.5, runif(length(x < 0.5), 0,0.5), runif(length(x >= 0.5), 0.5, 1))
#' y2 <- ifelse(x > 0.5, runif(length(x < 0.5), 0,0.5), runif(length(x >= 0.5), 0.5, 1))
#' A <- ECBC(x,y1, resolution = 50)
#' B <- ECBC(x,y2, resolution = 50)
#' #plot_density(A)
#' #plot_density(B)
#' D1.ECBC(A,B)


D1 <- function(x1,y1,x2,y2, resolution = NULL){
  if(is.null(resolution)){
    resolution <- floor(min(c(length(unique(x1)), length(unique(y1)),
                            length(unique(x2)), length(unique(y2))))^0.5)
  }else{
    resolution <- resolution
  }

  .X1 <- rank(x1, ties.method="max")
  .Y1 <- rank(y1, ties.method="max")
  .X2 <- rank(x2, ties.method="max")
  .Y2 <- rank(y2, ties.method="max")
  A1 <- round(build_checkerboard_weights(.X1, .Y1, resolution),15)
  A2 <- round(build_checkerboard_weights(.X2, .Y2, resolution),15)

  res = 0
  for (i in 1:resolution){
    yA_sum1 = 0
    yB_sum1 = 0
    for (j in 1:resolution){
      output = .r_local_kernel_integral_AB(A1,A2,i,j,yA_sum1, yB_sum1, resolution)
      res = res + output["yD"]
      yA_sum1 = output["yA_sum"]
      yB_sum1 = output["yB_sum"]
    }
  }
  return(unname(res/resolution))
}


#' @param A Numeric matrix of dimension NxN indicating the mass of the first N-checkerboad copula
#' @param B Numeric matrix of dimension NxN indicating the mass of the second N-checkerboad copula
#' @rdname D1
#' @return \code{D1.ECBC()} returns the D1-distance between to checkerboard copulas A and B with same resolution
#'

D1.ECBC <- function(A,B){
  #Check if the dimensions of the checkerboard copulas are the same
  if(any(dim(A) != dim(B))){
    warning("The dimensions of the checkerboard copulas have to coincide!")
  }else{
    N = NROW(A)
    res = 0
    for (i in 1:N) {
      yA_sum1 = 0
      yB_sum1 = 0
      for (j in 1:N) {
        output = .r_local_kernel_integral_AB(A,B,i,j,yA_sum1, yB_sum1, N)
        res = res + output["yD"]
        yA_sum1 = output["yA_sum"]
        yB_sum1 = output["yB_sum"]
      }
    }
    return(unname(res/N))
  }
}


#' @param X Numeric vector of values in the first coordinate
#' @param Y Numeric vector of values in the second coordinate
#' @param resolution integer indicating the resolution of the checkerboard copula. (default = NULL)
#' @rdname D1
#' @return \code{zeta1()} returns the directed dependence from x to y.
#'
zeta1 <- function(X,Y, resolution = NULL){
  #N <- floor(length(x)^0.5)
  #An <- ECBC(x,y, resolution = N)
  #zeta1 <- 3*D1(An, matrix(1/N^2, nrow = N, ncol = N))
  #return(zeta1)
  if(is.null(resolution)){
    sample_size <- min(length(unique(X)), length(unique(Y)))
    resolution <- floor(sample_size ^ (1/2))    #s \in [0,1/2)
  }
  .X <- rank(X, ties.method="max")
  .Y <- rank(Y, ties.method="max")
  A <- build_checkerboard_weights(.X, .Y, resolution)
  return(3*D1_Pi(A, resolution))
}





