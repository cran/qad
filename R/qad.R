#' Measure of (asymmetric and directed) dependence
#'
#' @description Quantification of (asymmetric and directed) dependence structures between two random variables X and Y.
#'
#' @rdname qad
#'
#' @param x a data.frame containing two columns with the observations of the bi-variate sample or a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values.
#' @param resolution an integer indicating the number of strips for the checkerboard aggregation (see \link{ECBC}). We recommend to use the default value (resolution = NULL)
#' @param p.value a logical indicating whether to return a p-value of rejecting independence (based on permutation).
#' @param nperm an integer indicating the number of permutation runs (if p.value = TRUE)
#' @param p.value_asymmetry a logical indicating whether to return a p-value for the measure of asymmetry (based on bootstrap).
#' @param nboot an integer indicating the number of runs for the bootstrap.
#' @param print a logical indicating whether the result of qad is printed.
#' @param remove.00 a logical indicating whether double 0 entries should be excluded (default = FALSE)
#' @param ... Further arguments passed to 'qad' will be ignored
#'
#' @details qad is the implementation of a strongly consistent estimator of the copula based dependence measure zeta_1 introduced in Trutschnig 2011.
#' We first compute the empirical copula of a two-dimensional sample, aggregate it to the so called empirical checkerboard copula (ECBC), and
#' calculate zeta_1 of the ECBC and its transpose. In order to test for independence (in both directions), a built-in p-value
#' is implemented (a permutation test with nperm permutation runs to estimate the p-value).
#' Furthermore, a bootstrap test with nboot runs can be applied to estimate a p-value for the measure of asymmetry a.
#'
#'
#' @return qad returns an object of class qad containing the following components:
#' \item{data}{ a data.frame containing the input data.}
#' \item{q(X,Y)}{influence of X on Y}
#' \item{q(Y,X)}{influence of Y on X}
#' \item{max.dependence}{maximal dependence}
#' \item{results}{ a data.frame containing the results of the dependence measures.}
#' \item{mass_matrix}{ a matrix containing the mass distribution of the empirical checkerboard copula.}
#' \item{resolution}{an integer containing the used resolution of the checkerboard aggregation.}
#' \item{n}{Sample size.}
#'
#' @references Trutschnig, W. (2011). On a strong metric on the space of copulas and its induced dependence measure, Journal of Mathematical Analysis and Applications 384, 690-705.
#'
#' Junker, R., Griessenberger, F. and Trutschnig, W. (2021). Estimating scale-invariant directed dependence of bivariate distributions.  Computational Statistics and Data Analysis, 153.
#'
#' @seealso
#' A tutorial can be found at \url{http://www.trutschnig.net/software.html}.
#'
#' @examples
#' #Example 1 (independence)
#'
#' n <- 100
#' x <- runif(n,0,1)
#' y <- runif(n,0,1)
#' sample <- data.frame(x,y)
#' qad(sample)
#'
#' ###
#'
#' #Example 2 (mutual complete dependence)
#'
#' n <- 500
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
#'
#' #Example 4 (Asymmetry)
#'
#' n <- 100
#' x <- runif(n,0,1)
#' y <- (2*x) %% 1
#' qad(x, y, p.value_asymmetry = TRUE)




qad <- function(x, ...){
  UseMethod("qad")
}

#' @rdname qad
#' @method qad data.frame
qad.data.frame <- function(x, resolution = NULL,
                           p.value = TRUE,
                           nperm = 1000,
                           p.value_asymmetry = FALSE,
                           nboot = 1000,
                           print=TRUE, remove.00 = FALSE,...){
  X <- x

  #Remove double zero entries if desired
  if(remove.00){
    X <- filter(X, rowSums(abs(X)) > 0)
  }

  #Remove NAs
  Z <- na.omit(X[,1:2])
  x <- as.numeric(Z[,1])
  y <- as.numeric(Z[,2])
  if(NROW(X) > NROW(Z)){
    warning(paste(NROW(X) - NROW(Z), ' observation(s) containing NAs were removed!'))
  }

  X_size <- NROW(Z)

  # Calculate the default resolution
  if(is.null(resolution)) {
    sample_size <- min(length(unique(x)), length(unique(y)))
    resolution <- floor(sample_size ^ (1/2))    #s \in [0,1/2)
  } else{
    sample_size <- NULL
    resolution <- floor(resolution)
  }

  .X <- rank(x, ties.method="max")
  .Y <- rank(y, ties.method="max")
  mass_matrix <- round(build_checkerboard_weights(.X, .Y, resolution),15) #Round to avoid computer errors

  #Influence of x on y
  zeta1 <- 3*D1_Pi(mass_matrix, resolution)
  #Influence of y on x
  zeta1.t <- 3*D1_Pi(t(mass_matrix), resolution)
  #Max influence
  max.dependence <- max(c(zeta1, zeta1.t))
  #Asymmetry
  asym <- zeta1 - zeta1.t

  #_____________________________________________________________________________
  #Compute the p-values:
  p_asym_perm <- p_zeta1 <- p_zeta1.t <- p_max_dep <- as.numeric(NA)

  # #Distribution of zeta1 under H0: independence
  # if(p.value && !permutation && !is.null(sample_size)){
  #   mc1 <- mc2 <- rep(0,nperm)
  #   for(i in 1:nperm){
  #     n <- length(x)#min(length(unique(x)), length(unique(y)))
  #     xx <- runif(n)
  #     yy <- runif(n)
  #     mc1[i] <- zeta1(xx,yy, resolution)
  #     mc2[i] <- zeta1(yy,xx, resolution)
  #   }
  #   p_zeta1 <- 1-ecdf(mc1)(zeta1)
  #   p_zeta1.t <- 1-ecdf(mc2)(zeta1.t)
  #   p_max_dep <- 1-ecdf(pmax(mc1,mc2))(max(c(zeta1.t,zeta1)))
  # }

  mass_matrix0 <- mass_matrix
  #Alternative pvalue: Permutation test for zeta1, zeta1.t and mean.dependence
  if(p.value){
    n <- length(.X)
    rx <- unname(table(.X)[as.character(.X)])
    ry <- unname(table(.Y)[as.character(.Y)])

    zeta1.perm <- zeta1.t.perm <- max.dependence.perm <- rep(0,nperm)
    for(i in 1:nperm){
      u_new <- as.numeric((.X-rx*runif(n,0,1))/n)
      v_new <- as.numeric((.Y-ry*runif(n,0,1))/n)
      samplePerm <- sample(c(u_new,v_new), size = 2*n)
      x1Perm <- samplePerm[1:n]
      x2Perm <- samplePerm[(n+1):(2*n)]

      .x1 <- rank(x1Perm, ties.method="max")
      .x2 <- rank(x2Perm, ties.method="max")
      mass_matrix <- build_checkerboard_weights(.x1, .x2, resolution)

      #Influence of x on y
      zeta1.perm[i] <- 3*D1_Pi(mass_matrix, resolution)
      #Influence of y on x
      zeta1.t.perm[i] <- 3*D1_Pi(t(mass_matrix), resolution)
      #Mean influence
      max.dependence.perm[i] <- max(c(zeta1.perm[i], zeta1.t.perm[i]))
    }

    p_zeta1 <- mean(ifelse(zeta1.perm >= zeta1, 1, 0))
    p_zeta1.t <- mean(ifelse(zeta1.t.perm >= zeta1.t, 1, 0))
    p_max_dep <- mean(ifelse(max.dependence.perm >= max.dependence,1,0))
  }

  #optional a p-value for asymmetry in dependence
  if(p.value_asymmetry){
    #Draw from empirical checkerboard copula
    n <- length(.X)
    zeta1.boot <- zeta1.t.boot <- asymmetry.boot <- rep(0,nboot)
    N <- resolution

    for(i in 1:nboot){
      ru_new <- sample(1:N, n, replace = T)
      tb_ru <- table(ru_new)
      rv_new <- sapply(names(tb_ru), function(k) sample(1:N, unname(tb_ru[k]), prob = N*mass_matrix0[as.numeric(k),], replace = T))
      ru_new <- sort(ru_new)
      rv_new <- unname(unlist(rv_new))
      u_new <- runif(n, 0, 1/N) + (ru_new-1)/N
      v_new <- runif(n, 0, 1/N) + (rv_new-1)/N

      .x1 <- rank(u_new, ties.method="max")
      .x2 <- rank(v_new, ties.method="max")
      mass_matrix <- build_checkerboard_weights(.x1, .x2, resolution)

      #Influence of x on y
      zeta1.boot[i] <- 3*D1_Pi(mass_matrix, resolution)
      #Influence of y on x
      zeta1.t.boot[i] <- 3*D1_Pi(t(mass_matrix), resolution)
      #Asymmetry
      asymmetry.boot[i] <- zeta1.boot[i]-zeta1.t.boot[i]
    }
    critical_value <- ecdf(asymmetry.boot)(0)
    if(critical_value < 0.5){
      critical_value <- critical_value
    }else if(critical_value == 1 & median(asymmetry.boot) == 0){
      critical_value <- 1/2
    }else{
      critical_value <- 1-critical_value
    }
    p_asym_perm <- 2*critical_value
  }


  #output
  names <- c('q(x1,x2)', 'q(x2,x1)','max.dependence','asymmetry      ')
  q.values <- c(zeta1, zeta1.t, max.dependence, asym)
  p.values <- c(p_zeta1, p_zeta1.t, p_max_dep, p_asym_perm)
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
  if(resolution <= 3){
    warning("Resolution is less or equal to 3. Results must be interpreted with caution!")
  }

  output_qad <- list(data = data.frame(x1=x,x2=y),
                     `q(X,Y)` = zeta1,
                     `q(Y,X)` = zeta1.t,
                     max.dependence = max.dependence,
                     results = output,
                     mass_matrix = mass_matrix0,
                     resolution = resolution,
                     n = sample_size)
  class(output_qad) <- 'qad'
  invisible(output_qad)
}


#' @rdname qad
#' @method qad numeric
qad.numeric <- function(x, y , resolution = NULL,
                        p.value = TRUE,
                        nperm = 1000,
                        p.value_asymmetry = FALSE, nboot = 1000,
                        print=TRUE, remove.00 = FALSE,...){
  X <- data.frame(x,y)
  names(X) <- c(deparse(substitute(x)),deparse(substitute(y)))
  return(qad.data.frame(X, resolution = resolution,
                        p.value = p.value,nperm = nperm,
                        p.value_asymmetry = p.value_asymmetry, nboot = nboot,
                        print = print, remove.00 = remove.00))
}



#' Pairwise quantification of (asymmetric and directed) dependencies
#'
#' Pairwise computation of the function \code{qad}(). \code{qad}() is applied on each pair of variables of a numeric data.frame.
#'
#'
#' @param data_df a data frame containing numeric columns with the observations of the sample.
#' @param remove.00 a logical indicating whether double 0 entries should be excluded (default = FALSE)
#' @param min.res an integer indicating the necessary minimum resolution of the checkerboard grid to compute qad, otherwise the result is NA (default = 3).
#' @param p.value a logical indicating whether to return a p-value of rejecting independence (based on permutation).
#' @param nperm an integer indicating the number of permutation runs.
#' @param p.value_asymmetry a logical indicating whether a p-value (based on bootstrap) is computed for the measure of asymmetry.
#' @param nboot an integer indicating the number of bootstrapping runs.
#'
#' @return a list, containing 8 data.frames with the dependence measures, corresponding p.values, the resolution of the checkerboard aggregation and the number of removed double zero entries (only if remove.00 = TRUE).
#' The output of pairwise.qad() can be illustrated using the function \code{heatmap.qad()}.
#'
#' @examples
#' n <- 100
#' x1 <- runif(n, 0, 1)
#' x2 <- x1^2 + rnorm(n, 0, 0.1)
#' x3 <- runif(n, 0, 1)
#' x4 <- x3 - x2 + rnorm(n, 0, 0.1)
#' sample_df <- data.frame(x1,x2,x3,x4)
#' #Fit qad
#' model <- pairwise.qad(sample_df, p.value = FALSE)
#' heatmap.qad(model, select = "dependence", fontsize = 6)


pairwise.qad <- function(data_df, remove.00 = FALSE, min.res = 3,
                         p.value = TRUE, nperm = 1000,
                         p.value_asymmetry = FALSE, nboot = 1000){



  #==== Data preparations ====#
  data_df <- data.frame(data_df)

  var_names <- colnames(data_df)
  n_var <- length(var_names)

  M <- data.frame(matrix(as.numeric(NA), nrow = n_var, ncol = n_var))
  colnames(M) <- row.names(M) <- var_names
  qM <- maxM <- aM <- M
  qM.pvalue <- maxM.pvalue <- aM.pvalue <- M
  resM <- uniqueranksM <- M
  N_00 <- M



  #==== Start with the pairwise computations ====#
  print('computation process...')

  for(i in 1:(n_var-1)){
    for(j in (i+1):n_var){

      pw_data <- na.omit(data_df[,c(i,j)])
      u1 <- length(unique(pw_data[[1]]))
      u2 <- length(unique(pw_data[[2]]))
      res <- floor(sqrt(min(u1,u2)))
      ties <- NROW(pw_data) - min(u1,u2)

      #==== Remove columns ====#
      if(remove.00){
        N_00[i,j] <- N_00[j,i] <- NROW(filter(pw_data, rowSums(abs(pw_data)) == 0))
        pw_data <- filter(pw_data, rowSums(abs(pw_data)) > 0)

        u1 <- length(unique(pw_data[[1]]))
        u2 <- length(unique(pw_data[[2]]))
        res <- floor(sqrt(min(u1,u2)))
      }

      n_ranks <- min(u1,u2)

      #==== Calculate the qad fits ====#

      if(res < min.res){
        qad_coefficients <- rep(as.numeric(NA), 8)
        names(qad_coefficients) <- c('q(x1,x2)', 'q(x2,x1)','max.dependence','asymmetry',
                                     'p.q(x1,x2)', 'p.q(x2,x1)','p.max.dependence','p.asymmetry')
        qad_res <- res
      }else{
        qad_fit <- qad(pw_data, p.value = p.value, nperm = nperm,
                       p.value_asymmetry = p.value_asymmetry, nboot = nboot,
                       print = FALSE)
        qad_coefficients <- coef(qad_fit)
        qad_res <- qad_fit$resolution
      }

      #q values
      qM[i,j] <- qad_coefficients[c('q(x1,x2)')]
      qM[j,i] <- qad_coefficients[c('q(x2,x1)')]
      #mean dependence values
      maxM[i,j] <- qad_coefficients[c('max.dependence')]
      maxM[j,i] <- qad_coefficients[c('max.dependence')]
      #asymmetry values
      aM[i,j] <- qad_coefficients[c('asymmetry')]
      aM[j,i] <- -qad_coefficients[('asymmetry')]
      #q p.values
      qM.pvalue[i,j] <- qad_coefficients[c('p.q(x1,x2)')]
      qM.pvalue[j,i] <- qad_coefficients[c('p.q(x2,x1)')]
      #mean dependence p.values
      maxM.pvalue[i,j] <- qad_coefficients[c('p.max.dependence')]
      maxM.pvalue[j,i] <- qad_coefficients[c('p.max.dependence')]
      #asymmetry p.values
      aM.pvalue[i,j] <- qad_coefficients[c('p.asymmetry')]
      aM.pvalue[j,i] <- qad_coefficients[c('p.asymmetry')]

      #Resolution
      resM[i,j] <- resM[j,i] <- qad_res
      uniqueranksM[i,j] <- uniqueranksM[j,i] <- n_ranks
    }
    print(paste('computation process:', i+1,'/',n_var))
  }



  return(list(q = qM,
              max.dependence = maxM,
              asymmetry = aM,
              q_p.values = qM.pvalue,
              max.dependence_p.values = maxM.pvalue,
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
#' @param n number of observations (or minimum of unique values, if ties occur).
#' @param R number of repetitions (default R = 1000)
#' @param resolution resolution of checkerboard copula (default  = NULL)
#'
#' @details The distribution of qad in the setting of independence, i.e., the random variables
#' X and Y are independent. The distribution is calculated in the following way: Samples of size n
#' are drawn from independent random variables. Then qad is calculated. The procedure is repeated R times. #'
#'
#'
#' @return \code{pqad} gives the distribution function, i.e. P(qad <= q). \code{qqad} gives the quantile function.
#' The length of the result is determined by the length of q or p, respectively.
#'
#' @examples
#' pqad(0.3, 45)
#' qqad(0.5, 30)



pqad <- function(q, n, R = 1000, resolution = NULL){
  if(is.null(resolution)){
    resolution <- floor(n^(1/2))
  }
  zeta1.values <- rep(0,R)
  for(i in 1:R){
    x <- runif(n)
    y <- runif(n)
    zeta1.values[i] <- zeta1(x,y,resolution = resolution)
  }
  return(ecdf(zeta1.values)(q))
}

#' @name qad_distribution
#' @rdname qad_distribution

qqad <- function(p, n, R = 1000, resolution = NULL){
  if(is.null(resolution)){
    resolution <- floor(n^(1/2))
  }
  zeta1.values <- rep(0,R)
  for(i in 1:R){
    x <- runif(n)
    y <- runif(n)
    zeta1.values[i] <- zeta1(x,y, resolution = resolution)
  }
  return(quantile(zeta1.values, probs = p, type = 1))
}



