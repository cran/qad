#' The empirical checkerboard copula
#'
#' The function \code{emp_c_copula()} computes the mass distribution of the
#' empirical (checkerboard) copula, given a bivariate sample. \code{emp_c_copula_eval()} evaluates the
#' the empirical (checkerboard) copula at given points.
#' If \code{smoothing} = FALSE, the empirical copula is computed (if there are ties in the sample an adjusted empirical copula is computed),
#' otherwise the empirical checkerboard copula - a smoothed version of the empirical copula - is computed. For more information of the calculations, see details.
#'
#' @param X a data frame with two columns containing the observations of the sample. Each row contains
#' one observation.
#' @param smoothing a logial indicating whether the checkerboard aggregation is computed (default = TRUE).
#' @param resolution an integer indicating the resolution of the checkerboard aggregation, i.e.
#' the number of vertical/horizontal strips of the checkerboard copula.
#'
#' @return \code{emp_c_copula()} returns a matrix with the mass distribution of the empirical
#' (checkerboard) copula.
#'
#' @details If the observations come from a distribution with continuous margins,
#' i.e. there are no ties in the sample, the function \code{emp_c_copula()} gives the same result
#' as the function \code{C.n()} in the \code{copula} package.
#' If there are ties in the sample, the empirical copula is adjusted and calculated in the following way: \cr
#' Let (u_i,v_i) := (F_n(x_i),G_n(y_i)) be the pseudo-observations for i in \{1,\ldots,n\} and (u_1',v_1'),\ldots, (u_m',v_m') the distinct pairs of pseudo-observations with m leq n. Moreover set S_1:=\{0, u_1, \ldots, u_{m_1}\} and S_2:=\{0, v_1,\ldots, v_{m_2}\} and define the quantities t_i,r_i and s_i for i=1,\ldots, m by
#' \deqn{t_i := sum_{j=1}^n 1_{(u_i',v_i')}(u_j,v_j)}
#' \deqn{r_i := sum_{j=1}^n 1_{u_i}(u_j)}
#' \deqn{s_i := sum_{j=1}^n 1_{v_i}(v_j)}
#' where 1 defines the indicator function.
#' Define the empirical subcopula A'_n: S_1 x S_2 to \{0,1/n, \ldots, (n-1)/n,1\} by
#' \deqn{A'_n(s_1,s_2)= 1/n  sum_{i=1}^m t_i * 1_{[0,s_1] x [0,s_2]} (u_i', v_i')=1/n sum_{i=1}^n 1_{[0,s_1] x [0,s_2]} (u_i, v_i)}
#' for all s_1 in S_1 and s_2 in S_2. \cr
#' We extend the subcopula A'_n to a copula by defining the transformations w_i:[0,1]^2 to [u_i'-r_i/n,u_i'] x [v_i'-s_i/n,v_i'] by
#' \deqn{w_i(x,y)=(u_i'-r_i/n+r_i*x/n, v_i'-s_i/n + s_iy/n)}
#' and set the measure of the empirical copula mu_{A_n}^B := 1/n sum_{i=1}^m t_i mu_B^{w_i}, where B denotes the product copula.
#'
#' The checkerboard aggregation is computed as usual (see references).
#'
#' @note The calculation of the empirical copula with a high sample size (and resolution rate) can take time.
#'
#' @references
#' Deheuvels, P. (1979). La fonction de dependance empirique et ses proprietas: un test non parametrique d'independance, Acad. Roy. Belg. Bull. Cl. Sci., 5th Ser. 65, 274-292.
#'
#' Li, X., Mikusinski, P. and Taylor, M.D. (1998). Strong approximation of copulas, Journal of Mathematical Analysis and Applications, 255, 608-623.
#'
#' Genest, C., Neshlehova J.G. and Remillard, B. (2014). On the empirical multilinear copula process for count data. Bernoulli, 20 (3), 1344-1371.

emp_c_copula <- function(X, smoothing = TRUE, resolution) {
  .Deprecated("ECBC")
  x <- as.numeric(data.frame(X[, 1:2])[, 1])
  y <- as.numeric(data.frame(X[, 1:2])[, 2])

  X <- na.omit(X[,1:2])

  if(smoothing == FALSE){
    N <- NROW(X)
  }else{
    N <- floor(resolution)
  }

  #Determine the grid wrt the resolution
  grid <- seq(0, 1, length.out = N + 1)
  grid <- expand.grid(grid,grid)
  names(grid) <- c('x','y')

  #if there are no ties we can use the function C.n from the copula package
  if (length(x) == length(unique(x)) &
      length(y) == length(unique(y))) {
    z <- copula::C.n(as.matrix(grid), as.data.frame(X), smoothing = 'checkerboard')
    M <- matrix(z, nrow = N + 1, ncol = N + 1)
    mass <-
      M[-1,-1] + M[-(N + 1),-(N + 1)] - M[-1,-(N + 1)] - M[-(N + 1),-1]
  } else{
    #Calculate the ranks and the frequency of each rank of the sample
    rank_x <- data.frame(rank_x = rank(x, ties.method = 'max'))
    tab_x <- rank_x %>%
      dplyr::group_by(rank_x) %>%
      dplyr::summarise(freq_rank_x = n())

    #tab_x <- rank_x[, .(freq_rank_x = .N), by = rank_x]
    rank_y <- data.frame(rank_y = rank(y, ties.method = 'max'))
    #tab_y <- rank_y[, .(freq_rank_y = .N), by = rank_y]
    tab_y <- rank_y %>%
      dplyr::group_by(rank_y) %>%
      dplyr::summarise(freq_rank_y = n())

    #Create data.frame with rank informations
    rank_df <- data.frame(rank_x, rank_y)
    #rank_df <- rank_df[, .(freq = .N), by = .(rank_x, rank_y)]
    rank_df <- rank_df %>%
      dplyr::group_by(rank_x, rank_y) %>%
      dplyr::summarise(freq = n())

    rank_df <- dplyr::left_join(rank_df, tab_x, by = 'rank_x')
    rank_df <- dplyr::left_join(rank_df, tab_y, by = 'rank_y')

    #Calculate the mass and the rectangles with min and max
    n <- length(x)
    rank_df$xmax <- rank_df$rank_x / n
    rank_df$xmin <- rank_df$xmax - rank_df$freq_rank_x / n
    rank_df$ymax <- rank_df$rank_y / n
    rank_df$ymin <- rank_df$ymax - rank_df$freq_rank_y / n
    rank_df$mass <- (1 / n) * rank_df$freq
    #setkey(rank_df, xmin, ymin)   #This line does not go through check --> note: global variable

    #mass function for the checkerboard copula
    mass_function <- function(u, v, rank_df, n) {
      #A <- rank_df[xmin <= u & ymin <= v,]   #This line does not go through check --> note: global variable
      A <- rank_df[rank_df$xmin <= u & rank_df$ymin <= v, ]
      output <- sum(A$mass * (pmin(A$xmax,u)-A$xmin)/(A$xmax-A$xmin)*(pmin(A$ymax,v)-A$ymin)/(A$ymax-A$ymin))
      return(output)
    }

    z <-
      mapply(
        mass_function,
        u = grid$x,
        v = grid$y,
        MoreArgs = list(rank_df = rank_df, n = n)
      )

    M <- matrix(z, nrow = N+1, ncol = N+1)
    mass <- M[-1, -1] + M[-(N+1), -(N+1)] - M[-1, -(N+1)] - M[-(N+1), -1]
  }
  return(round(mass,14))
}

#' @param u a data.frame with two columns containing the evaluation points. Each row consists of a x and y value.
#' @rdname emp_c_copula
#' @return \code{emp_c_copula_eval()} returns a vector of evaluations of the empirical
#' (checkerboard) copula.

emp_c_copula_eval <- function(X, u, smoothing = TRUE, resolution){
  .Deprecated("ECBC.eval")

  if(inherits(u, 'data.frame')){
    mass <- emp_c_copula(X, smoothing = smoothing, resolution)
    N <- NROW(mass)
    u <- data.frame(u)
    #Set grid for evaluation
    grid <- seq(0,1,length.out = N+1)
    mass_new <- matrix(0, nrow= N+2, ncol = N+2)
    mass_new[2 : (N+1), 2 : (N+1)] <- mass

    mass <- t(apply(apply(mass_new,2,cumsum),1,cumsum))

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
    return(bilinear_interpolation(grid, mass, u[,1], u[,2]))
  }else{
    warning('u must be a data.frame')
  }

}
