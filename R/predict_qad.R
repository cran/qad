#' Predict conditional probabilities
#'
#' @description The function \code{predict.qad()} can be utilized to predict the probabilities of the event that Y lies in
#' specific intervals given X=x, or vice versa. Additionally, the mass of the conditional distribution functions are plotted.
#' The prediction can be computed in the sample setting as well as in the copula setting (pseudo-observation in the unit square).
#'
#' @param object an object of class 'qad', which determines the underlying checkerboard aggregation.
#' @param values a vector containing the x or the y values for which the conditional probabilities should be predicted.
#' @param conditioned a character specifying on which variable is conditioned. Options are "x1" (default) or "x2".
#' @param nr_intervals an integer, which determines a different number of intervals for the prediction (only possible in the copula setting).
#' @param copula a logical (default =FALSE) determining whether the empirical checkerboard copula is used or the retransformed data.
#' @param prediction_interval a vector specifying the interval boundaries for which the conditional probability is computed. Options are NULL (default) to predict the conditional probabilites for all intervals or a vector c(lower_boundary, upper_boundary) indicating the boundaries.
#' @param pred_plot a logical indicating if the conditional probabilites are plotted.
#' @param panel.grid a logical indicating whether the panel.grid is plotted.
#' @param ... some methods for this generic require additional arguments.  None are used in this method.
#'
#' @return a list containing a data.frame with the interval boundaries and the prediction probabilities and a plot depicting the mass of the conditional distributions functions.
#'
#' @note Predictions are only possible for values within the range of the sample (or between 0 and 1 in the copula setting). For given values exceeding the range NA is returned.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- runif(n, -4 ,4)
#' y <- x^2 + rnorm(n, 0, 1)
#' sample <- data.frame(x, y)
#'
#' ##(Not Run)
#' qad.fit <- qad(sample)
#' predict.qad(qad.fit, values = c(-2,0.6), conditioned = "x1", pred_plot = TRUE)
#' predict.qad(qad.fit, values = c(1,9), conditioned = "x2", pred_plot = TRUE)
#' predict.qad(qad.fit, values = c(-2,0.6), conditioned = "x1", pred_plot = FALSE,
#'         nr_intervals = 4)
#' predict.qad(qad.fit, values = c(-2,0.6), conditioned = "x1", pred_plot = FALSE,
#'             prediction_interval = c(4,6))
#' predict.qad(qad.fit, values = c(4,0.6), conditioned = "x2", pred_plot = FALSE,
#'             prediction_interval = c(2,3))
#'
#' qad.pred <- predict.qad(qad.fit, values = c(-2,0.6), conditioned = "x1", pred_plot = FALSE)
#' qad.pred$prediction
#' qad.pred$plot
#'
#'
#' @export predict qad
predict.qad <- function(object,
                        values,
                        conditioned = "x1",
                        nr_intervals = NULL,
                        prediction_interval = NULL,
                        copula = FALSE,
                        pred_plot = FALSE,
                        panel.grid = TRUE,  ...) {


  if(copula){
    pred_cop <- .predict_qad_copula(values = values, conditioned = conditioned, qad_output = object,
                       nr_intervals = nr_intervals, prediction_interval = prediction_interval)
    pred_prob <- pred_cop$pred_prob
    grid <- seq(0,1,length.out = object$resolution+1)
    index <- pred_cop$indices

    #Prediction plot
    p <- plot.qad(object, copula = TRUE, panel.grid = panel.grid)

    if(conditioned == 'x1'){
      df_rect <- data.frame(xmin = grid[index],
                            xmax = grid[index+1],
                            ymin = min(pred_prob$lowerBound),ymax = max(pred_prob$upperBound))
    }else if(conditioned == "x2"){
      df_rect <- data.frame(ymin = grid[index],
                            ymax = grid[index+1],
                            xmin = min(pred_prob$lowerBound), xmax = max(pred_prob$upperBound))
    }
    p <- p + geom_rect(data = df_rect, aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), color = 'red', fill = NA, size = 1)

  }else{
      if(!is.null(nr_intervals)){
        warning("nr_intervals only works in the copula setting. Checkerboard resolution is used instead.")
      }
      grid <- seq(0,1,length.out = object$resolution + 1)
      GridX <- quantile(object$data$x1, grid, type = 1)
      GridY <- quantile(object$data$x2, grid, type = 1)
      if(conditioned == "x1"){
        x_index <- as.numeric(cut(values, GridX, labels = 1:(length(GridX)-1), include.lowest = T))
        new_values <- (grid[x_index] + grid[x_index + 1])/2
        cInv <- "x2"
        Grid <- GridX
      }else if(conditioned == "x2"){
        y_index <- as.numeric(cut(values, GridY, labels = 1:(length(GridY)-1), include.lowest = T))
        new_values <- (grid[y_index] + grid[y_index + 1])/2
        cInv <- "x1"
        Grid <- GridY
      }
      pred_cop <- .predict_qad_copula(values = new_values, conditioned = conditioned, qad_output = object,
                                            nr_intervals = NULL, prediction_interval = NULL)

      pred_prob <- pred_cop$pred_prob
      index <- pred_cop$indices

      if(length(which(values > max(object$data[conditioned]) | values < min(object$data[conditioned])))>0){
        pred_prob[,3+which(values > max(object$data[conditioned]) | values < min(object$data[conditioned]))] <- NA
        exceeding <- which(values > max(object$data[conditioned]) | values < min(object$data[conditioned]))
      }

      names(pred_prob) <- c("Interval", "lowerBound", "upperBound", paste0(conditioned, "=",values))
      if(conditioned == "x1"){
        pred_prob$lowerBound <- GridY[1:NROW(pred_prob)]
        pred_prob$upperBound <- GridY[2:(NROW(pred_prob)+1)]
      }else if(conditioned == "x2"){
        pred_prob$lowerBound <- GridX[1:NROW(pred_prob)]
        pred_prob$upperBound <- GridX[2:(NROW(pred_prob)+1)]
      }

      if(!is.null(prediction_interval)){
        pred_prob_new <- data.frame(Interval = "I", lowerBound = prediction_interval[1], upperBound = prediction_interval[2])
        for(i in 4:(NCOL(pred_prob))){
          f_int <- approxfun(c(pred_prob$lowerBound[1],pred_prob$upperBound), c(0,cumsum(pred_prob[,i])), method = "linear", yleft = 0, yright = 1)
          pred_prob_new[1,names(pred_prob)[i]] <- f_int(prediction_interval[2])-f_int(prediction_interval[1])
        }
        pred_prob <- pred_prob_new
      }


    #Prediction plot
    p <- plot.qad(object, copula = FALSE, panel.grid = panel.grid)
    if(conditioned == 'x1'){
      df_rect <- data.frame(xmin = Grid[index],
                            xmax = Grid[index+1],
                            ymin = min(pred_prob$lowerBound),ymax = max(pred_prob$upperBound))
      if(length(which(values > max(object$data[conditioned]) | values < min(object$data[conditioned])))>0){
        df_rect[exceeding, ] <- rep(NA, 4)
      }
    }else if(conditioned == "x2"){
      df_rect <- data.frame(ymin = Grid[index],
                            ymax = Grid[index+1],
                            xmin = min(pred_prob$lowerBound), xmax = max(pred_prob$upperBound))
      if(length(which(values > max(object$data[conditioned]) | values < min(object$data[conditioned])))>0){
        df_rect[exceeding, ] <- rep(NA, 4)
      }
    }
    p <- p + geom_rect(data = df_rect, aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax), color = 'red', fill = NA, size = 1)

  }

  print_prob <- data.frame(t(pred_prob[,4:NCOL(pred_prob)]))
  names(print_prob) <- pred_prob$Interval
  row.names(print_prob) <- names(pred_prob)[-c(1:3)]
  cat("Intervall definition:\n")
  print.data.frame(pred_prob[,1:3])
  cat("\nPrediction probabilities:\n")
  print.data.frame(print_prob)

  if(pred_plot){
    print(p)
  }

  return(invisible(list(prediction = pred_prob,
              plot = p)))

}



