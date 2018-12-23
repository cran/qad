#' Predict conditional probabilities
#'
#' @description The function \code{predict.qad()} predicts the probabilities
#' to end up in specific intervals given x or y values and plots the conditional probabilities.
#' The prediction can be computed in the copula setting
#' or in the data setting.
#'
#' @param object an object of class 'qad', which determines the underlying checkerboard aggregation.
#' @param values a vector containing the x or the y values for which the conditional probabilities should be predicted.
#' @param conditioned a character specifying on which variable is conditioned. Options are "x1" (default) or "x2".
#' @param nr_intervals an integer, which determines the number of intervals for the prediction. Note, that in the copula setting
#' the intervals are equidistant, in the data setting the retransformed intervals have different lengths. (default = NULL: the number of intervals is the
#' resolution of the checkerboard copula)
#' @param copula a logical (default =FALSE) determining whether the empirical checkerboard copula is used or the retransformed data.
#' @param prediction_interval a vector specifying the interval boundaries for which the conditional probability is computed. Options are NULL (default) to predict the conditional probabilites for all intervals or a vector c(lower_boundary, upper_boundary) indicating the boundaries.
#' @param pred_plot a logical indicating if the conditional probabilites are plotted.
#' @param panel.grid a logical indicating whether the panel.grid is plotted.
#' @param ... some methods for this generic require additional arguments.  None are used in this method.
#'
#' @return a named data.frame and a plot (optional). Each row stands for an evaluation point and the columns contain the conditional probabilities of the intervals.
#'
#' @note Predictions are only possible for values within the range of the sample (or between 0 and 1 in the copula setting). Values exceeding the range are removed.
#'
#' @examples
#' n <- 1000
#' x <- runif(n, -1 ,1)
#' y <- x^2 + rnorm(n, 0, 1)
#' sample <- data.frame(x, y)
#'
#' ##(Not Run)
#' #mod <- qad(sample)
#' #val <- c(-0.5, 0,1)
#' #predict(mod, values = val, conditioned = "x1", copula = FALSE, pred_plot = TRUE)
#' #predict(mod, values = val, conditioned = "x1", copula = TRUE)
#' #predict(mod, values = val, conditioned = "x1", copula = TRUE, pred_plot = TRUE)
#'
#'
#' @export predict qad
predict.qad <- function(object, values, conditioned = c("x1","x2"), nr_intervals = NULL, prediction_interval = NULL, copula = FALSE,  pred_plot = FALSE, panel.grid = TRUE,  ...) {
  qad_output <- object
    #Consider the copula setting
    if(copula){
      xpobs <- values[values>=0 & values <= 1]
      prediction <- .predict_qad_pseudoobservations(values = xpobs, conditioned, qad_output, nr_intervals, prediction_interval)

      #Prediction plot
      if(pred_plot){
        grid <- seq(0 ,1 ,length.out = qad_output$resolution + 1)
        index <- sapply(rownames(prediction), function(x){return(min(which(1:qad_output$resolution/qad_output$resolution >= x)))})
        p <- plot.qad(qad_output, copula = TRUE, panel.grid = panel.grid)
        if(conditioned == 'x1'){
          df_rect <- data.frame(xmin = grid[index],xmax = grid[index+1],ymin = 0,ymax = 1)
          p <- p + geom_rect(data = df_rect, aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), color = 'red', fill = NA)
        }else{
          df_rect <- data.frame(ymin = grid[index],ymax = grid[index+1],xmin = 0,xmax = 1)
          p <- p + geom_rect(data = df_rect, aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), color = 'red', fill = NA)
        }
        print(p)
      }

      return(prediction)

    #Consider the data setting
    }else{

      #x := unique values to predict
      x <- unique(values)
      #conditioned on which variable
      if(conditioned[1] == 'x2'){
        x2 <- qad_output$data$x1
        x1 <- qad_output$data$x2
      }else{
        x1 <- qad_output$data$x1
        x2 <- qad_output$data$x2
      }
      x <- x[x >= min(x1) & x <= max(x1)]

      resolution <- qad_output$resolution
      #Consider possible duplicated x values before we transform to pseudoobservations
      duplicated_x_values <- which(duplicated(ecdf(x1)(x)))
      x_pobs <- unique(ecdf(x1)(x))

      if(!is.null(prediction_interval)){
        prediction_interval_pobs <- c(ecdf(x2)(prediction_interval[1]),ecdf(x2)(prediction_interval[2]))
      }else{
        prediction_interval_pobs <- NULL
      }
      #Calculate the predictions for the pseudoobservations
      df_predict <- .predict_qad_pseudoobservations(x_pobs, conditioned[1], qad_output, nr_intervals, prediction_interval = prediction_interval_pobs)

      #Define the row names
      if(length(duplicated_x_values) > 0){
        rownames(df_predict) <- x[-duplicated_x_values]
      }else{
        rownames(df_predict) <- x
      }

      if(!is.null(prediction_interval)){
        colnames(df_predict) <- paste('[',prediction_interval[1],',',prediction_interval[2],']',sep='')
      }else{
      #Consider equidistant case
        if(is.null(nr_intervals)){
          N <- resolution
        } else{
          N <- nr_intervals
        }
        colnames(df_predict) <- paste('(',round(quantile(x2, 0:(N-1)/N),2),',',round(quantile(x2,1:N/N),2),']', sep='')
      }

      #Prediction plot
      if(pred_plot){
        grid <- seq(0 ,1 ,length.out = qad_output$resolution + 1)
        p <- plot.qad(qad_output, copula = FALSE, panel.grid = panel.grid)
        if(conditioned == 'x1'){
          lower <- round(quantile(x1, grid[-length(grid)]),2)
          upper <- round(quantile(x1, grid[-1]),2)
          index <- sapply(as.numeric(row.names(df_predict)), function(x){return(min(which(upper >= x)))})
          df_rect <- data.frame(xmin = lower[index],xmax = upper[index],ymin = min(qad_output$data$x2),ymax = max(qad_output$data$x2))
          p <- p + geom_rect(data = df_rect, aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), color = 'red', fill = NA)
        }else{
          lower <- round(quantile(x1, grid[-length(grid)]),2)
          upper <- round(quantile(x1, grid[-1]),2)
          index <- sapply(as.numeric(row.names(df_predict)), function(x){return(min(which(upper >= x)))})
          df_rect <- data.frame(ymin = lower[index],ymax = upper[index],xmin = min(qad_output$data$x1),xmax = max(qad_output$data$x1))
          p <- p + geom_rect(data = df_rect, aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), color = 'red', fill = NA)
        }
        print(p)
      }

      return(df_predict)
    }
}
