#' Calculate empirical checkerboard copula
#'
#' The function \code{ECBC()} computes the mass distribution of the
#' empirical (checkerboard) copula, given a bi-variate sample.
#' If resolution equals sample size, the bi-linearly extended empirical copula is returned.
#' Note, if there are ties in the sample an adjusted empirical copula is calculated. The function
#' \code{ECBC.eval()} evaluates the checkerboard copula at given points.
#'
#' @param X Numeric vector of values in the first coordinate
#' @param Y Numeric vector of values in the second coordinate
#' @param resolution Integer indicating the resolution of the checkerboard aggregation, i.e.
#' the number of vertical/horizontal strips of the checkerboard copula. (default = NULL, sets resolution = floor(sqrt(sample size)))
#'
#' @return \code{ECBC()} returns a matrix with the mass distribution of the empirical
#' (checkerboard) copula.
#'
#' @details If the observations are drawn from a continuous distribution (no ties in the sample),
#' the function \code{ECBC()} returns the commonly used empirical checkerboard copula.
#' If there are ties in the sample, the empirical copula is adjusted and calculated in the following way: \cr
#' Let \eqn{(u_i,v_i) := (F_n(x_i),G_n(y_i))} be the pseudo-observations for \eqn{i \in \{1,\ldots,n\}} and \eqn{(u_1',v_1'),\ldots, (u_m',v_m')} the distinct pairs of pseudo-observations with m leq n. Moreover set \eqn{S_1:=\{0, u_1, \ldots, u_{m_1}\}} and \eqn{S_2:=\{0, v_1,\ldots, v_{m_2}\}} and define the quantities \eqn{t_i,r_i} and
#' \eqn{s_i} for \eqn{i=1,\ldots, m} by
#' \deqn{t_i := \sum_{j=1}^n 1_{(u_i',v_i')}(u_j,v_j)}
#' \deqn{r_i := \sum_{j=1}^n 1_{u_i}(u_j)}
#' \deqn{s_i := \sum_{j=1}^n 1_{v_i}(v_j)}
#' where 1 defines the indicator function.
#' Define the empirical subcopula \eqn{A'_n: S_1 x S_2 \to \{0,1/n, \ldots, (n-1)/n,1\}} by
#' \deqn{A'_n(s_1,s_2)= 1/n  \sum_{i=1}^m t_i * 1_{[0,s_1] x [0,s_2]} (u_i', v_i')=1/n \sum_{i=1}^n 1_{[0,s_1] x [0,s_2]} (u_i, v_i)}
#' for all \eqn{s_1 \in S_1} and \eqn{s_2 in S_2}. \cr
#' We extend the subcopula \eqn{A'_n} to a copula by defining the transformations \eqn{w_i:[0,1]^2 \to [u_i'-r_i/n,u_i'] x [v_i'-s_i/n,v_i']} by
#' \deqn{w_i(x,y)=(u_i'-r_i/n+r_i*x/n, v_i'-s_i/n + s_iy/n)}
#' and set the measure of the empirical copula \eqn{mu_{A_n}^B := 1/n \sum_{i=1}^m t_i mu_B^{w_i}}, where B denotes the product copula.
#'
#'
#' @references
#' Deheuvels, P. (1979). La fonction de dependance empirique et ses proprietas: un test non parametrique d'independance, Acad. Roy. Belg. Bull. Cl. Sci., 5th Ser. 65, 274-292.
#'
#' Li, X., Mikusinski, P. and Taylor, M.D. (1998). Strong approximation of copulas, Journal of Mathematical Analysis and Applications, 255, 608-623.
#'
#' Genest, C., Neshlehova J.G. and Remillard, B. (2014). On the empirical multilinear copula process for count data. Bernoulli, 20 (3), 1344-1371.
#'
#' Junker, R.R., Griessenberger, F. and Trutschnig, W. (2021). Estimating scale-invariant directed dependence of bivariate distributions. Computational Statistics and Data Analysis, 153, 107058.
#'
#' @examples
#' ##Generate data drawn from the product copula and compute the empirical (checkerboard) copula
#' n <- 100
#' x <- runif(n, 0, 1)
#' y <- runif(n, 0, 1)
#' mass <- ECBC(x,y, resolution = 10)
#' plot_density(mass)
#' mass <- ECBC(x,y, resolution = n)
#' plot_density(mass)
#'
#' ## Compute empirical checkerboard copula of a sample with ties and plot density
#' n <- 100
#' x <- sample(runif(n, -1, 1), n, replace=TRUE)
#' y <- x^2 + rnorm(n, 0, 1)
#' mass <- ECBC(x,y)
#' plot_density(mass)


ECBC <- function(X,Y, resolution = NULL){
  if(is.null(resolution)){
    resolution <- floor(min(c(length(unique(X)), length(unique(Y))))^0.5)
  }else if(resolution > length(X)){
    resolution <- length(X)
    warning(paste0("Resolution cannot exceed sample size. Automatically set to ", length(X)))
  }else if(resolution < 1){
    resolution <- 1
    warning(paste0("Resolution has to be positive! Automatically set to 1, which returns independence copula!"))
  }else{
    resolution <- resolution
  }
  X <- rank(X, ties.method = "max")
  Y <- rank(Y, ties.method = "max")
  sample_size = length(X)
  result = matrix(0, nrow=resolution, ncol=resolution)
  x_range = .r_range(X)
  y_range = .r_range(Y)
  for (i in 1:sample_size) {
    rx = x_range[X[i]]
    ry = y_range[Y[i]]
    x_upper = ceiling(X[i] / sample_size * resolution)
    x_lower = max(ceiling((X[i] - rx) / sample_size * resolution), 1)
    y_upper = ceiling(Y[i] / sample_size * resolution)
    y_lower = max(ceiling((Y[i] - ry) / sample_size * resolution), 1)
    for (x in x_lower:x_upper) {
      lambda_x = min(X[i], x / resolution * sample_size) - max(X[i] - rx, (x-1) / resolution * sample_size)
      for (y in y_lower:y_upper) {
        lambda_y = min(Y[i], y / resolution * sample_size) - max(Y[i] - ry, (y-1) / resolution * sample_size)
        result[x,y] = result[x,y] + lambda_x*lambda_y/(sample_size*rx*ry)
      }
    }
  }
  return(round(result,15))
}


#' @param CB A numeric mass matrix of a checkerboard copula (ECBC)
#' @param eval.points A numeric matrix or data.frame indicating the eval.points (x,y)
#' @rdname ECBC
ECBC.eval <- function(CB, eval.points){
  # N = NROW(CB)
  # x_rect <- as.numeric(cut(eval.points[,1], breaks = seq(0,1,length.out = N+1), labels = 1:N, right = T, include.lowest = T))
  # y_rect <- as.numeric(cut(eval.points[,2], breaks = seq(0,1,length.out = N+1), labels = 1:N, right = T, include.lowest = T))
  # M = rbind(0,cbind(0,CB))
  # eval <- rep(0, NROW(eval.points))
  # for(i in 1:NROW(eval.points)){
  #   Sum1 = sum(M[1:x_rect[i], 1:y_rect[i]])
  #   Sum2 = sum(M[x_rect[i] + 1, 1:y_rect[i]] * (eval.points[i,1] - (x_rect[i]-1)/N)*N)
  #   Sum3 = sum(M[1:x_rect[i], y_rect[i] + 1] * (eval.points[i,2] - (y_rect[i]-1)/N)*N)
  #   Sum4 = sum(M[x_rect[i] + 1, y_rect[i] + 1] * (eval.points[i,2] - (y_rect[i]-1)/N) * (eval.points[i,1] - (x_rect[i]-1)/N)*N^2)
  #   eval[i] = Sum1 + Sum2 + Sum3 + Sum4
  # }
  # return(eval)

  res <- NROW(CB)
  g <- seq(0,1,length.out = res + 1)
  mass <- rbind(0,cbind(0,CB,0),0)
  mass <- t(apply(apply(mass,2,cumsum),1,cumsum))

  bilinear_interpolation <- function(grid, z, x0, y0){
    index_NA <- which(x0 < 0 | x0 > 1 | y0 < 0 | y0 >1)
    x0 <- ifelse(x0 < 0 | x0 > 1, 0, x0)
    y0 <- ifelse(y0 < 0 | y0 > 1, 0, y0)

    N <- length(grid) - 1
    lbx_nr <- sapply(x0, function(x) max(which(grid <= x)))
    lby_nr <- sapply(y0, function(x) max(which(grid <= x)))
    grid <- c(grid,2*grid[length(grid)]-grid[length(grid)-1])
    lbx <- x0 - grid[lbx_nr]
    rbx <- grid[lbx_nr + 1] - x0

    lby <- y0 - grid[lby_nr]
    rby <- grid[lby_nr + 1] - y0

    output <- z[cbind(lbx_nr, lby_nr)] * rbx * rby +
      z[cbind(lbx_nr + 1, lby_nr)] * lbx * rby +
      z[cbind(lbx_nr, lby_nr + 1)] * rbx * lby +
      z[cbind(lbx_nr + 1, lby_nr + 1)] * lbx * lby
    output <- output*N^2
    output[index_NA] <- NA
    return(output)
  }
  return(bilinear_interpolation(g, mass, eval.points[,1], eval.points[,2]))
}

