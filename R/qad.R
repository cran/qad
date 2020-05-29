#' Measure of (asymmetric and directed) dependence
#'
#' @description Quantification of (asymmetric and directed) dependence structures between two random variables X and Y.
#'
#' @rdname qad
#'
#' @param x a data.frame containing two columns with the observations of the bivariate sample or a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values.
#' @param resolution an integer indicating the number of strips for the checkerboard aggregation (see \link{emp_c_copula}).
#' Default = NULL uses the optimal resolution.
#' @param permutation a logical indicating whether a p-value (based on permutations) is computed; otherwise a p-value is computed on MC-simulations (see pqad()).
#' @param nperm an integer indicating the number of permutation runs.
#' @param DoParallel a logical value indicating whether the repetitions in the permutation test is computed parallel.
#' @param registerC function to register the parallel environment. It is recommended to use registerDoParallel(), contained in the doParallel package (default). Another option, especially for a linux based system, is to install the
#' doMC package and use registerDoMC
#' @param ncores an integer indicating the number of cores used for parallel computation. (Default = NULL, which is defined by max(cores)-1)
#' @param print a logical indicating whether the result of qad is printed.
#' @param remove.00 a logical indicating whether double 0 entries should be excluded (default = FALSE)
#' @param ... Further arguments passed to 'qad' will be ignored
#'
#' @details qad is the implementation of a strongly consistent estimator of the copula based dependence measure zeta_1 introduced in Trutschnig 2011.
#' We first compute the empirical copula of a two-dimensional sample, aggregate it to the so called empirical checkerboard copula (ECB), and
#' calculate zeta_1 of the ECB copula and its transpose. In order to test for independence (in both directions), the distribution (and hence the p-value)
#' of a Monte-Carlo simulation is provided (default). Alternatively, a permutation test can be used to obtain p-values for the direction and asymmetry.
#'
#'
#' @note The computation of the p-values (aggregated by permutations) take some time.
#'
#' @return qad returns an object of class qad containing the following components:
#' \item{data}{ a data.frame containing the input data.}
#' \item{results}{ a data.frame containing the results of the dependence measures.}
#' \item{mass_matrix}{ a matrix containing the mass distribution of the empirical checkerboard copula.}
#' \item{resolution}{an integer containing the used resolution of the checkerboard aggregation.}
#'
#' @references Trutschnig, W. (2011). On a strong metric on the space of copulas and its induced dependence measure, Journal of Mathematical Analysis and Applications 384, 690-705.
#'
#' @examples
#' #Example 1 (independence)
#'
#' n <- 1000
#' x <- runif(n,0,1)
#' y <- runif(n,0,1)
#' sample <- data.frame(x,y)
#' qad(sample)
#'
#' ###
#'
#' #Example 2 (mutual complete dependence)
#'
#' n <- 1000
#' x <- runif(n,0,1)
#' y <- x^2
#' sample <- data.frame(x,y)
#' qad(sample)
#'
#' #Example 3 (complete dependence)
#'
#' n <- 1000
#' x <- runif(n,-10,10)
#' y <- sin(x)
#' sample <- data.frame(x,y)
#' qad(sample)



qad <- function(x, ...){
  UseMethod("qad")
}

#' @rdname qad
#' @method qad data.frame
qad.data.frame <- function(x, resolution = NULL, permutation = FALSE, nperm=1000, DoParallel=TRUE, registerC = registerDoParallel, ncores=NULL, print=TRUE, remove.00 = FALSE,...){
  X <- x

  if(remove.00){
    X <- filter(X, rowSums(abs(X)) > 0)
  }

  x <- as.numeric(data.frame(na.omit(X[,1:2]))[,1])
  y <- as.numeric(data.frame(na.omit(X[,1:2]))[,2])

  if(NROW(X) > NROW(na.omit(X[,1:2]))){
    warning(paste(NROW(X) - NROW(na.omit(X)), ' observation(s) containing NAs were removed!'))
  }

  X <- na.omit(X[,1:2])
  X_size <- NROW(X)

  check_pvalue <- resolution
  # Calculate the default resolution
  if (is.null(resolution)) {
    unique_x <- length(unique(x))
    unique_y <- length(unique(y))
    sample_size <- min(unique_x, unique_y)
    resolution <- floor(sample_size ^ (1/2))    #s \in [0,1/2)
  } else{
    resolution <- floor(resolution)
  }

  dist <- .zeta1(X,resolution=resolution)
  ##y ~ x
  zeta1 <- dist$zeta1
  ##x ~ y
  zeta1.t <- dist$zeta1.t
  ##mean dependence
  mean.dependence <- mean(c(zeta1, zeta1.t))
  ##asymmetry
  asym <- zeta1 - zeta1.t
  ##Resolution and mass matrix
  resolution <- dist$resolution.checkerboard
  mass_matrix <- dist$mass_matrix


  p_asym_perm <- p_zeta1 <- p_zeta1.t <- p_mean_dep <- as.numeric(NA)

  #MC simulation to get the p-value
  if(is.null(check_pvalue)){
    p_zeta1 <- ifelse(round(zeta1,8) == 0, 1, 1 - .ppqad(zeta1, sample_size))
    p_zeta1.t <- ifelse(round(zeta1.t,8) == 0, 1, 1 - .ppqad(zeta1.t, sample_size))
    p_mean_dep <- ifelse(round(mean.dependence,8) == 0, 1, 1 - .ppmqad(mean.dependence, sample_size))
  }

  #Permutation test to obtain a p-value
  if(permutation == TRUE){
    if(DoParallel){
      #Detect cores for parallelisation
      if(is.null(ncores)){
        cores <- parallel::detectCores()-1
      }else{
        cores <- ncores
      }

      #Register parallel backend
      registerC(cores)


      xAll <- c(x,y)
      permutation_loop <- foreach(i=1:nperm, .combine='rbind', .packages = c('data.table','copula'), .export=c('emp_c_copula', '.markov_kernel','.zeta1_checkerboard_strip','.zeta1')) %dopar% {
        xPerm <- sample(xAll, size = length(xAll), replace = FALSE)
        x1Perm <- xPerm[1:(length(xAll)/2)]
        x2Perm <- xPerm[((length(xAll)/2)+1):length(xAll)]

        sample <-data.table(x1Perm,x2Perm)
        dist <- .zeta1(X=sample, resolution=resolution)
        ##y ~ x
        zeta1.perm <- dist$zeta1
        ##x ~ y
        zeta1.t.perm <- dist$zeta1.t
        ##asymmetry
        asym.perm <- abs(zeta1.perm - zeta1.t.perm)
        ##mean dependence
        mean_dep.perm <- mean(c(zeta1.perm, zeta1.t.perm), na.rm = TRUE)

        return(c(zeta1.perm,zeta1.t.perm, asym.perm, mean_dep.perm))
      }

      if(.Platform$OS.type == "windows"){
        stopImplicitCluster()
      }

      zeta1.perm <- permutation_loop[,1]
      zeta1.t.perm <- permutation_loop[,2]
      asym.perm <- permutation_loop[,3]
      mean_dep.perm <- permutation_loop[,4]

    }else{
      xAll <- c(x,y)
      zeta1.perm <- zeta1.t.perm <- asym.perm <- mean_dep.perm <- rep(NA, nperm)
      for(i in 1:nperm){
        xPerm <- sample(xAll, size = length(xAll), replace = FALSE)
        x1Perm <- xPerm[1:(length(xAll)/2)]
        x2Perm <- xPerm[((length(xAll)/2)+1):length(xAll)]

        sample <-data.table(x1Perm,x2Perm)
        dist <- .zeta1(X=sample, resolution=resolution)
        ##y ~ x
        zeta1.perm[i] <- dist$zeta1
        ##x ~ y
        zeta1.t.perm[i] <- dist$zeta1.t
      }
      ##asymmetry
      asym.perm <- abs(zeta1.perm - zeta1.t.perm)
      ##mean dependence
      mean_dep.perm <- apply(cbind(zeta1.perm, zeta1.t.perm),1,mean, na.rm = TRUE)
    }

    p_asym_perm <- mean(ifelse(asym.perm >= abs(asym),1,0))
    p_zeta1 <- mean(ifelse(zeta1.perm >= abs(zeta1),1,0))
    p_zeta1.t <- mean(ifelse(zeta1.t.perm >= abs(zeta1.t),1,0))
    p_mean_dep <- mean(ifelse(mean_dep.perm >= abs(mean.dependence),1,0))
  }


  #output
  names <- c('q(x1,x2)', 'q(x2,x1)','mean.dependence','asymmetry      ')
  q.values <- c(zeta1, zeta1.t, mean.dependence, asym)
  p.values <- c(p_zeta1, p_zeta1.t, p_mean_dep, p_asym_perm)
  output <- data.frame(names, q.values, p.values)
  names(output) <- c('','coef','p.values')

  output_q <- output[1:3,]
  names(output_q) <- c('','q','p.values')
  output_q[,2:3] <- round(output_q[,2:3],3)
  output_a <- output[4,]
  names(output_a) <- c('','a','p.values')
  output_a[,2:3] <- round(output_a[,2:3],3)

  if(print){
    cat("\n")
    cat("quantification of asymmetric dependence:", "\n")
    cat("\nData: x1 :=", paste(colnames(X)[1]))
    cat("\n      x2 :=", paste(colnames(X)[2]))
    cat("\n")
    cat(paste("\nSample Size:"),X_size)
    cat(paste("\nNumber of unique ranks:", "x1:", length(unique(x))))
    cat(paste("\n                        x2:", length(unique(y))))
    cat(paste("\n                   (x1,x2):", NROW(unique(X))))
    cat(paste("\nResolution:",resolution,'x',resolution))
    cat("\n\nDependence measures:")
    cat("\n")
    if(all(is.na(output$p.values))){
      print.data.frame(format(output_q[,1:2], justify='left', digits=3), row.names = FALSE)
      cat("\n")
      print.data.frame(format(output_a[,1:2], justify='left', digits=3), row.names = FALSE)
    }else{
      print.data.frame(format(output_q, justify='left', digits=3), row.names = FALSE)
      cat("\n")
      print.data.frame(format(output_a, justify='left', digits=3), row.names = FALSE)
    }
  }

  output_qad <- list(data = data.frame(x1=x,x2=y),
                     results = output,
                     mass_matrix = mass_matrix,
                     resolution = resolution)
  class(output_qad) <- 'qad'
  invisible(output_qad)
}


#' @rdname qad
#' @method qad numeric
qad.numeric <- function(x, y , resolution = NULL, permutation=FALSE, nperm = 1000, DoParallel = TRUE, registerC = registerDoParallel, ncores = NULL, print = TRUE,remove.00 = FALSE,...){
  X <- data.frame(x,y)
  names(X) <- c(deparse(substitute(x)),deparse(substitute(y)))
  return(qad.data.frame(X, resolution = resolution, permutation = permutation, nperm = nperm, DoParallel = DoParallel, registerC = registerC, ncores = ncores, print = print, remove.00 = remove.00))
}



#' Pairwise quantification of (asymmetric and directed) dependencies
#'
#' Pairwise computation of the function \code{qad}(). \code{qad}() is applied on each pair of variables of a numeric data.frame.
#'
#'
#' @param data_df a data frame containing numeric columns with the observations of the sample.
#' @param resolution an integer indicating the number of strips for the checkerboard aggregation (see \link{emp_c_copula}()).
#' Default (NULL) uses the optimal resolution, computed out of the sample size.
#' @param remove.00 a logical indicating whether double 0 entries should be excluded (default = FALSE)
#' @param min.res an integer indicating the necessary minimum resolution of the checkerboard grid to compute qad, otherwise the result is NA (default = 3).
#' @param permutation a logical indicating whether a p-value (based on permutations) is computed; (otherwise the p-value is computed by MC-simulation - see pqad()).
#' @param nperm an integer indicating the number of permutation runs.
#' @param DoParallel a logical value indicating whether the permutation test is computed parallelized.
#' @param registerC function to register the parallel backend. It is recommended to use registerDoParallel() of the doParallel package (default). Other option is for example on a linux based system to install the
#' doMC package and use registerDoMC
#' @param ncores an integer indicating the number of cores used for parallelization. Default (NULL) uses the maximum number of cores minus 1.
#'
#' @return a list, containing 8 data.frames with the dependence measures, corresponding p.values, the resolution of the checkerboard aggregation and the number of removed double zero entries (only if remove.00 = TRUE).
#' The output of pairwise.qad() can be illustrated using the function \code{heatmap.qad()}.
#'
#' @examples
#' n <- 100
#' x <- runif(n, 0, 1)
#' y <- runif(n, 0, 1)
#' z <- runif(n, 0, 1)
#' sample_df <- data.frame(x,y,z)
#'
#' #qad
#' model <- pairwise.qad(sample_df, permutation = FALSE)
#' heatmap.qad(model, select = "dependence", fontsize = 5, significance = TRUE, sign.level = 0.05)


pairwise.qad <- function(data_df, resolution = NULL, remove.00 = FALSE, min.res = 3,
                         permutation = FALSE, nperm = 1000, DoParallel = FALSE, registerC = registerDoParallel, ncores = NULL){



  #==== Data preparations ====#
  data_df <- data.frame(data_df)

  var_names <- colnames(data_df)
  n_var <- length(var_names)

  M <- data.frame(matrix(as.numeric(NA), nrow = n_var, ncol = n_var))
  colnames(M) <- row.names(M) <- var_names
  qM <- mdM <- aM <- M
  qM.pvalue <- mdM.pvalue <- aM.pvalue <- M
  resM <- uniqueranksM <- M
  N_00 <- M



  #==== Start with the pairwise computations ====#
  print('computation process...')

  for(i in 1:(n_var-1)){
    for(j in (i+1):n_var){

      pw_data <- data_df[,c(i,j)]
      u <- dplyr::n_distinct(pw_data)
      u1 <- dplyr::n_distinct(pw_data[[1]], na.rm = TRUE)
      u2 <- dplyr::n_distinct(pw_data[[2]], na.rm = TRUE)
      res <- floor(sqrt(min(u,u1,u2)))
      ties <- NROW(pw_data) - min(u,u1,u2)

      #==== Remove columns ====#
      if(remove.00){
        N_00[i,j] <- N_00[j,i] <- NROW(filter(pw_data, rowSums(abs(pw_data)) == 0))
        pw_data <- filter(pw_data, rowSums(abs(pw_data)) > 0)

        u <- dplyr::n_distinct(pw_data)
        u1 <- dplyr::n_distinct(pw_data[[1]], na.rm = TRUE)
        u2 <- dplyr::n_distinct(pw_data[[2]], na.rm = TRUE)
        res <- floor(sqrt(min(u,u1,u2)))
      }

      n_ranks <- min(u,u1,u2)

      #==== Calculate the qad fits ====#

      if(res < min.res){
        qad_coefficients <- rep(as.numeric(NA), 8)
        names(qad_coefficients) <- c('q(x1,x2)', 'q(x2,x1)','mean.dependence','asymmetry',
                                     'p.q(x1,x2)', 'p.q(x2,x1)','p.mean.dependence','p.asymmetry')
        qad_res <- res
      }else{
        qad_fit <- qad(pw_data, resolution = resolution,
                       permutation = permutation, nperm = nperm,
                       DoParallel = DoParallel, registerC = registerC,
                       ncores = ncores, print = FALSE)
        qad_coefficients <- coef(qad_fit)
        qad_res <- qad_fit$resolution
      }

      #q values
      qM[i,j] <- qad_coefficients[c('q(x1,x2)')]
      qM[j,i] <- qad_coefficients[c('q(x2,x1)')]
      #mean dependence values
      mdM[i,j] <- qad_coefficients[c('mean.dependence')]
      mdM[j,i] <- qad_coefficients[c('mean.dependence')]
      #asymmetry values
      aM[i,j] <- qad_coefficients[c('asymmetry')]
      aM[j,i] <- -qad_coefficients[('asymmetry')]
      #q p.values
      qM.pvalue[i,j] <- qad_coefficients[c('p.q(x1,x2)')]
      qM.pvalue[j,i] <- qad_coefficients[c('p.q(x2,x1)')]
      #mean dependence p.values
      mdM.pvalue[i,j] <- qad_coefficients[c('p.mean.dependence')]
      mdM.pvalue[j,i] <- qad_coefficients[c('p.mean.dependence')]
      #asymmetry p.values
      aM.pvalue[i,j] <- qad_coefficients[c('p.asymmetry')]
      aM.pvalue[j,i] <- qad_coefficients[c('p.asymmetry')]

      #Resolution
      resM[i,j] <- resM[j,i] <- qad_res
      uniqueranksM[i,j] <- uniqueranksM[j,i] <- n_ranks
    }
    print(paste('computation process:', i,'/',n_var))
  }



  return(list(q = qM,
              mean.dependence = mdM,
              asymmetry = aM,
              q_p.values = qM.pvalue,
              mean.dependence_p.values = mdM.pvalue,
              asymmetry_p.values = aM.pvalue,
              resolution = resM,
              #n_distinct_ranks = uniqueranksM,
              n_removed_00 = N_00))
}










#' Distribution of qad (H0: independence)
#'
#' Distribution function - P_H0(qad <= q) - and quantile function for the qad distribution with regard
#' to the null hypthesis (H0) stating independence between X and Y.
#'
#' @name qad_distribution
#' @rdname qad_distribution
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#'
#' @details The distribution of qad was computed in the setting of independence
#' between the random variables X and Y in the following way:
#'
#' For n < 1000, Monte Carlo (MC) simulation of H0 with 20.000 repetitions were executed
#' for each sample size. According to these values the empirical cumulative distribution
#' functions and the quantile functions were computed and then approximated on a coarser grid.
#'
#' For n >= 1000, MC simulations were executed again, but this time on a coarser sample size grid (steps of 100) until the size of 10.000.
#' The so obtained quantiles were approximated using the parametric function a*n^b+c, whereby
#' the parameters a,b,c were estimated using the R-function nls. Using the so calculated quantiles,
#' the empirical distribution function and the quantile functions were approximated.
#'
#'
#' @return \code{pqad} gives the distribution function, i.e. P(qad <= q). \code{qqad} gives the quantile function.
#' The length of the result is determined by the length of q or p, respectively.
#'
#' @examples
#' pqad(0.3, 45)
#' qqad(0.5, 30)



pqad <- function(q, n){
  data <- mcData_independence$q
  data_values <- data$qad_values
  data_param <- data$param

  if(n < 4){
    return(rep(1,length(q)))
  }else if(n < 1000){
    grid <- seq(data_values[[as.character(n)]]$min_grid, data_values[[as.character(n)]]$max_grid, length.out = data_values[[as.character(n)]]$grid_length)
    result <- approxfun(x = grid,
                        y = data_values[[as.character(n)]]$q_n_approx_values,
                        yleft = 0, yright = 1)(q)
    return(result)
  }else{
    #Sort is included, otherwise there are wrong quantiles for very big n
    result <- sort(as.numeric(apply(data_param, 1, function(p) p[1]*n^p[2] + p[3])))
    return(approxfun(result, seq(0, 1, length.out = NROW(data_param)),
                     yleft = 0, yright = 1)(q))
  }
}

#' @name qad_distribution
#' @rdname qad_distribution

qqad <- function(p, n){
  data <- mcData_independence$q
  data_values <- data$qad_values
  data_param <- data$param

  if(n < 4){
    return(rep(0,length(p)))
  }else if(n < 1000){
    grid <- seq(data_values[[as.character(n)]]$min_grid, data_values[[as.character(n)]]$max_grid, length.out = data_values[[as.character(n)]]$grid_length)
    q_ecdf <- approxfun(grid,
                        data_values[[as.character(n)]]$q_n_approx_values,
                        yleft = 0, yright = 1)
    pseudoinverse <- function(y){
      stats::uniroot(function(x) {
        q_ecdf(x) - max(y, q_ecdf(0))
      }, lower = 0, upper = 1)$root
    }
    result <- unlist(sapply(p, pseudoinverse))
    return(result)
  }else{
    #Sort is included, otherwise there are wrong quantiles for very big n
    values <- sort(as.numeric(apply(data_param, 1, function(p) p[1]*n^p[2] + p[3])))
    result <- approxfun(seq(0, 1, length.out = NROW(data_param)), values, yleft = 0, yright = 1)(p)
    return(result)
  }
}



#' Conditional confidence interval
#'
#' Conditioned on the sample size n, approximated confidence intervals of the dependence measure qad(x,y) for independent random variables are computed.
#' \code{cci()} can be used to test two random variables for independence.
#'
#' @param n an integer indicating the sample size.
#' @param alternative character string specifying the type of the confidence interval; must be one of "one.sided" (default) or "two.sided".
#'
#' @details \code{alternative = "one.sided"} provides a one-sided confidence interval, which can be interpreted that
#' in 95 of 100 realizations of two independet random variables X and Y
#' the calculated dependence measure qad(x,y) is less than the upper interval boundary.
#' If \code{alternative = "two.sided"} 95 of 100 realizations lie in between the interval boundaries.
#'
#'
#' @return a named vector indicating the lower and upper boundary of the confidence interval.


cci <- function(n, alternative = c("one.sided", "two.sided")){
  .Deprecated("qqad")
  if(alternative[1] == 'one.sided'){
    #lower: 0 - upper: 0.95
    #Estimated coefficients
    a <-  1.1056475
    b <- -0.2617012
    ci_upper <- ifelse(n <= 1, 0, a * n^b)
    ci <- c('"0"' = 0, '"0.95"' = ci_upper)
    return(ci)
  } else if(alternative[1] == 'two.sided'){
    #lower: 0.025  -  upper: 0.975
    #Estimated coefficients
    a_u <- 1.1406510
    b_u <- -0.2623318
    ci_upper <- ifelse(n <= 1, 0, a_u * n^b_u)

    a_l <- 0.5586709
    b_l <- -0.2046032
    ci_lower <- ifelse(n <= 1, 0, a_l * n^b_l)

    ci <- c('"0.025"' = ci_lower, '"0.975"' = ci_upper)
    return(ci)
  }else{
    warning('Choose an appropriate alternative: c("one.sided","two.sided")')
  }
}
