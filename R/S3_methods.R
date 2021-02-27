####################################################################################################################################
### Filename:    S3methods.R
### Description: S3 Generic Methods providing different ways to print the results
###
####################################################################################################################################

#' Summarize a qad object
#' @description Summary and coefficients of a \code{qad} output. The function
#' \code{summary()} prints the dependence measures, sample size and resolution of
#' the checkerboard copula and returns a list with the mentioned values.
#' The function \code{coef()} returns a named vector with the selected values.
#'
#' @param object an object of class 'qad'
#' @param ... some methods for this generic require additional arguments. None are used in this method.
#'
#' @return an object containing the calculated values of a \code{qad} object.
#'
#' @examples
#' n <- 100
#' x <- runif(n, 0, 1)
#' y <- x^2 + rnorm(n, 0, 0.1)
#' sample <- data.frame(x, y)
#' ##(Not Run)
#' fit <- qad(sample)
#' summary(fit)
#' coef(fit)
#' coef(fit, select = c('q(x1,x2)','p.q(x1,x2)'))
#'
#' @method summary qad
summary.qad <- function(object, ...) {
  x <- object
  output <- x$results
  X <- x$data
  resolution <- x$resolution

  output_q <- output[1:3,]
  names(output_q) <- c('','q','p.values')
  output_a <- output[4,]
  names(output_a) <- c('','a','p.values')

  cat("\n")
  cat("quantification of asymmetric dependencies:", "\n")
  cat("\n")
  cat(paste("\nSample Size:"),NROW(X))
  cat(paste("\nNumber of unique ranks:", "x1:", length(unique(X$x1))))
  cat(paste("\n                        x2:", length(unique(X$x2))))
  cat(paste("\n                   (x1,x2):", NROW(unique(X))))
  cat(paste("\nResolution:",resolution,'x',resolution))
  cat("\n\nDependence measures:")
  cat("\n")
  if(all(is.na(output$p.values))){
    print.data.frame(format(output_q[,1:2],justify='left', digits=4), row.names = FALSE)
    cat("\n")
    print.data.frame(format(output_a[,1:2],justify='left', digits=3), row.names = FALSE)
  }else{
    print.data.frame(format(output_q,justify='left', digits=4), row.names = FALSE)
    cat("\n")
    print.data.frame(format(output_a,justify='left', digits=3), row.names = FALSE)
  }
  invisible(list(SampleSize=NROW(x$data), resolution=x$resolution, dependence_values=x$results))
}


#' @rdname summary.qad
#' @method coef qad
#' @param select a vector of strings indicating which dependence measure should be returned. Options are c('q(x1,x2)', 'q(x2,x1)', 'max.dependence', 'asymmetry')
coef.qad <- function(object, select = c('q(x1,x2)', 'q(x2,x1)', 'max.dependence', 'asymmetry',
                                        'p.q(x1,x2)','p.q(x2,x1)','p.max.dependence','p.asymmetry'), ...){
  results <- object$results
  coef_values <- c('q(x1,x2)' = results[1,2],
                   'q(x2,x1)' = results[2,2],
                   'max.dependence' = results[3,2],
                   'asymmetry' = results[4,2],
                   'p.q(x1,x2)' = results[1,3],
                   'p.q(x2,x1)' = results[2,3],
                   'p.max.dependence' = results[3,3],
                   'p.asymmetry' = results[4,3])
  return(coef_values[select])
}










