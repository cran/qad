#' Plot conditional probabilites
#'
#' @description  Plots conditional probabilities for each strip of the checkerboard copula in the copula setting or in the retransformed data setting.
#'
#'
#' @param x an object of class qad.
#' @param addSample a logical determining whether the observations (or pseudo-observations) are added in the plot (default = FALSE).
#' @param copula a logical indicating whether the plot depicts the conditional probabilities
#' of the empirical checkerboard copula or of the retransformed data setting (default = FALSE).
#' @param density a logical indicating whether the density should be plotted instead of the conditional probabilites (default = FALSE).
#' Only works in the copula setting, i.e. if copula = TRUE.
#' @param margins a logical indicating whether the margin distribution is added in form of a rug plot.
#' @param point.size a numeric specifying the point size of the sample (relevant if addSample = TRUE).
#' @param panel.grid a logical indicating whether the panel grid is plotted. (default = TRUE)
#' @param ... some methods for this generic require additional arguments.  None are used in this method.
#'
#' @note The conditional probabilities are constant at squares in the copula setting. If the squares are retransformed in the data setting, the resulting objects are rectangles.
#'
#' @examples
#' ## Example 1
#' n <- 1000
#' x <- runif(n, 0, 1)
#' y <- runif(n, 0, 1)
#' sample <- data.frame(x, y)
#'
#' #qad (Not Run)
#' # mod <- qad(sample)
#' # plot(mod, addSample = TRUE, copula = FALSE)
#'
#' ## Example 2
#' n <- 1000
#' x <- runif(n, -1, 1)
#' y <- x^2 + rnorm(n, 0, 0.1)
#' sample <- data.frame(x, y)
#'
#' #qad (Not Run)
#' # mod <- qad(sample)
#' # plot(mod, addSample = TRUE, copula = TRUE)
#' # plot(mod, addSample = TRUE, copula = FALSE)
#'
#' @method plot qad

plot.qad <- function(x, addSample = FALSE, copula = FALSE, density = FALSE, margins = FALSE, point.size = 0.8, panel.grid = TRUE, ...){
  qad_output <- x
  if(class(qad_output)=='qad'){
#Copula setting
    if(copula){
      mass_matrix <- qad_output$mass_matrix

      if(density){
        density_matrix <- mass_matrix*NROW(mass_matrix)*NCOL(mass_matrix)
        density_matrix <- ifelse(density_matrix == 0, NA, density_matrix)
        data <- melt(density_matrix)
        n <- NROW(density_matrix)
      }else{
        cond_prob_matrix <- mass_matrix*NROW(mass_matrix)
        cond_prob_matrix <- ifelse(cond_prob_matrix==0, NA, cond_prob_matrix)
        data <- melt(cond_prob_matrix)
        n <- NROW(cond_prob_matrix)
      }

      names(data) <- c('x1','x2','value')
      sample_size <- length(qad_output$data$x1)

      p <- ggplot()
      p <- p + geom_raster(data=data, aes_(x=~x1/n-1/n/2,y=~x2/n-1/n/2 , fill=~value), alpha=0.95, na.rm = TRUE)
      if(density){
        p <- p + scale_fill_gradient(low="#cde7f0",high='#006699',guide = guide_colorbar(title='Density'), na.value=NA)
      }else{
        p <- p + scale_fill_gradient(low="#cde7f0",high='#006699',guide = guide_colorbar(title='Conditional probabilities'), na.value=NA)
      }
      p <- p + theme_bw() + xlab('x1') + ylab('x2') + ggtitle('Empirical checkerboard copula')
      if(addSample){
        p <- p + geom_point(data=qad_output$data, aes(x=rank(x1, ties.method = 'max')/sample_size,
                                                      y=rank(x2, ties.method = 'max')/sample_size), size=point.size)
      }
      if(margins){
        p <- p + geom_rug(data=qad_output$data, aes(x=rank(x1, ties.method = 'max')/sample_size,
                                                    y=rank(x2, ties.method = 'max')/sample_size))
      }
      if(!panel.grid){
        p <- p + theme(panel.grid.minor = element_line(colour = NA),
                       panel.grid.major = element_line(colour = NA))
      }
      p <- p + theme(panel.grid = element_line(linetype = 'dashed'))
      p <- p + scale_x_continuous(breaks = function(x) pretty(x, 15))
      p <- p + scale_y_continuous(breaks = function(x) pretty(x, 15))

      return(p)
#Data setting
    } else{
      x1 <- qad_output$data$x1
      x2 <- qad_output$data$x2
      N <- qad_output$resolution

      #middle point of each strip (checkerboard)
      x_pobs <- seq(0,1-(1/N),length.out = N) + 1/(2*N)

      #prediction of each segment (for conditional probabilites)
      df_predict <- .predict_qad_pseudoobservations(values = x_pobs, conditioned = 'x1', qad_output = qad_output, nr_intervals = N)

      #Define the rectangles
      raster_x <- data.table(xmin=round(quantile(x1, 0:(N-1)/N),2), xmax=round(quantile(x1,1:N/N),2), id=1)
      raster_y <- data.table(ymin=round(quantile(x2, 0:(N-1)/N),2), ymax=round(quantile(x2,1:N/N),2), id=1)
      raster_grid <- plyr::join(raster_x,raster_y, by='id')

      #Calculate the conditional probabilites for each rectangle
      raster_grid$prob <- as.vector(t(as.matrix(df_predict)))
      raster_grid$prob <- ifelse(raster_grid$prob == 0, NA, raster_grid$prob)

      p <- ggplot()
      p <- p + geom_rect(data=raster_grid, aes_(xmin=~xmin, xmax=~xmax, ymin=~ymin, ymax=~ymax,fill=~prob), alpha=0.95, na.rm = TRUE)
      if(margins){
        p <- p + geom_rug(data = qad_output$data, aes_(x=~x1, y=~x2))
      }
      p <- p + scale_fill_gradient(low="#cde7f0",high='#006699',guide = guide_colorbar(title='Conditional probabilites'), na.value=NA)
      if(addSample){
        p <- p + geom_point(data=qad_output$data, aes_(x=~x1,y=~x2), size=point.size)
      }
      p <- p + theme_bw() + xlab('x1') + ylab('x2') + ggtitle('Conditional probabilities')
      if(!panel.grid){
        p <- p + theme(panel.grid.minor = element_line(colour = NA),
                       panel.grid.major = element_line(colour = NA))
      }
      p <- p + theme(panel.grid = element_line(linetype = 'dashed'))
      p <- p + scale_x_continuous(breaks = function(x) pretty(x, 15))
      p <- p + scale_y_continuous(breaks = function(x) pretty(x, 15))
      return(p)
    }
  }
}

#_____________________________________________________________________________________________________________________

#' Plot density of empirical checkerboard copula
#'
#' @description  Plots the density/mass of the empirical checkerboard copula.
#'
#' @param mass_matrix a squared matrix containing the mass distribution, e.g. output of the function \code{emp_c_copula()}.
#' @param density a logical (TRUE = default) whether the density or the mass is plotted.
#'
#' @return a density plot (or mass distribution)
#' @examples
#' n <- 1000
#' x <- runif(n,0,1)
#' y <- runif(n,0,1)
#' sample <- data.frame(x,y)
#' plot(sample)
#'
#' mass <- emp_c_copula(sample, resolution=8)
#' plot_density(mass, density=TRUE)
#' plot_density(mass, density=FALSE)

plot_density <- function(mass_matrix, density=TRUE){
  #input: matrix with mass distribution
  if(density){
    mass_matrix <- mass_matrix*NROW(mass_matrix)*NCOL(mass_matrix)
  }
  mass_matrix <- ifelse(mass_matrix==0, NA, mass_matrix)
  data <- melt(mass_matrix)
  names(data) <- c('x1','x2','value')
  n <- NROW(mass_matrix)
  p <- ggplot(data=data, aes_(x=~x1/n-1/n/2,y=~x2/n-1/n/2))
  p <- p + geom_raster(aes_(fill=~value), alpha=0.95, na.rm = TRUE)
  if(density){
    p <- p + scale_fill_gradient(low="#cde7f0",high='#006699',guide = guide_colorbar(title='Density'), na.value=NA)
  }else{
    p <- p + scale_fill_gradient(low="#cde7f0",high='#006699',guide = guide_colorbar(title='Mass'), na.value=NA)
  }
  p <- p + theme_bw() + xlab('x1') + ylab('x2') + ggtitle('Empirical checkerboard copula')
  p
}



#' Heatmap of dependence measures
#'
#' @description  The pairwise computed dependence measures (output of the function \code{pairwise.qad()}) are illustrated by a heatmap.
#'
#' @param pw_qad output of the function \code{pairwise.qad}().
#' @param select a character indicating which dependence value is plotted.
#' Options are c("dependence", "mean.dependence", "asymmetry").
#' @param fontsize a numeric specifying the font size of the values.
#' @param significance a logical indicating whether significant values - with respect to the permutated p.values - are marked with a star.
#' @param sign.level numeric value indicating the significance level.
#' @param scale character indicating whether the heatmap uses a relative or absolute scale. Options are "rel" or "abs" (default).
#'
#' @details If the output of \code{pairwise.qad}() contains p-values, significant values can be highlighted by stars by setting significance=TRUE.
#'
#' @return a heatmap
#' @examples
#' n <- 1000
#' x <- runif(n, 0, 1)
#' y <- x^2 + rnorm(n, 0, 1)
#' z <- runif(n, 0, 1)
#' sample_df <- data.frame(x, y, z)
#'
#' #qad (Not Run)
#' #model <- pairwise.qad(sample_df, permutation = TRUE, nperm = 10, DoParallel = TRUE)
#' #heatmap.qad(model, select = "dependence", fontsize = 10, significance = TRUE)

heatmap.qad <- function(pw_qad, select = c('dependence','mean.dependence','asymmetry'), fontsize = 4, significance = TRUE, sign.level = 0.05, scale = "abs"){
  #Dependence measures
  if(select == 'dependence'){
    melt_df <- melt(as.matrix(pw_qad$q))
    melt_df$value <- round(melt_df$value, 2)
    p <- ggplot(data = melt_df, aes_(x = ~Var2, y = ~Var1, fill = ~value))
    p <- p + geom_tile(color = 'white')
    if(scale == "rel"){
      p <- p + scale_fill_gradient2(low = '#006699', high = '#006699', mid = "white",
                                    midpoint = 0, space = "Lab",
                                    name = "Dependence: qad(x,y)", na.value = "lightgrey")
    }else{
      p <- p + scale_fill_gradient2(low = '#006699', high = '#006699', mid = "white",
                                    midpoint = 0, limit = c(0, 1), space = "Lab",
                                    name = "Dependence: qad(x,y)", na.value = "lightgrey")
    }
    p <- p + theme_bw() + coord_fixed() + xlab('Variable 2') + ylab("Variable 1")
    p <- p + scale_x_discrete(position = 'top')
    p <- p + scale_y_discrete(limits = rev(unique(melt_df$Var1)))
    if(significance & any(!is.na(pw_qad$q_p.values))){
      melt_df_sign <- melt(as.matrix(pw_qad$q_p.values))
      melt_df_sign$value <- ifelse(melt_df_sign$value < sign.level & !is.na(melt_df_sign$value), '*','')
      melt_df$sign <- melt_df_sign$value
      melt_df$paste_value <- ifelse(is.na(melt_df$value), NA, paste(melt_df$value,melt_df$sign, sep=''))
      p <- p + geom_text(data = melt_df, aes_(x = ~Var2, y = ~Var1, label = ~paste_value), color = "black", size = fontsize, na.rm=TRUE)
    }else{
      p <- p + geom_text(aes_(x = ~Var2, y = ~Var1, label = ~value), color = "black", size = fontsize, na.rm = TRUE)
    }
    p



  }else if(select == 'mean.dependence'){
    melt_df <- melt(as.matrix(pw_qad$mean.dependence))
    melt_df$value <- round(melt_df$value, 2)
    p <- ggplot(data = melt_df, aes_(x = ~Var2, y = ~Var1, fill = ~value))
    p <- p + geom_tile(color = 'white')
    if(scale == "rel"){
      p <- p + scale_fill_gradient2(low = '#006699', high = '#006699', mid = "white",
                                    midpoint = 0, space = "Lab",
                                    name = "Mean dependence", na.value = "lightgrey")
    }else{
      p <- p + scale_fill_gradient2(low = '#006699', high = '#006699', mid = "white",
                                    midpoint = 0, limit = c(0, 1), space = "Lab",
                                    name = "Mean dependence", na.value = "lightgrey")
    }
    p <- p + theme_bw() + coord_fixed() + xlab('Variable 2') + ylab("Variable 1")
    p <- p + scale_x_discrete(position = 'top')
    p <- p + scale_y_discrete(limits = rev(unique(melt_df$Var1)))
    if(significance & any(!is.na(pw_qad$mean.dependence_p.values))){
      melt_df_sign <- melt(as.matrix(pw_qad$mean.dependence_p.values))
      melt_df_sign$value <- ifelse(melt_df_sign$value < sign.level & !is.na(melt_df_sign$value), '*','')
      melt_df$sign <- melt_df_sign$value
      melt_df$paste_value <- ifelse(is.na(melt_df$value), NA, paste(melt_df$value,melt_df$sign, sep=''))
      p <- p + geom_text(data = melt_df, aes_(x = ~Var2, y = ~Var1, label = ~paste_value), color = "black", size = fontsize, na.rm=TRUE)
    }else{
      p <- p + geom_text(aes_(x = ~Var2, y = ~Var1, label = ~value), color = "black", size = fontsize, na.rm=TRUE)
    }
    p



  }else if(select == 'asymmetry'){
    melt_df <- melt(as.matrix(pw_qad$asymmetry))
    melt_df$value <- round(melt_df$value, 2)
    p <- ggplot(data = melt_df, aes_(x = ~Var2, y = ~Var1, fill = ~value))
    p <- p + geom_tile(color = 'white')
    if(scale == "rel"){
      p <- p + scale_fill_gradient2(low = '#006699', high = '#006699', mid = "white",
                                    midpoint = 0, space = "Lab",
                                    name = "Asymmetry", na.value = "lightgrey")
    }else{
      p <- p + scale_fill_gradient2(low = '#006699', high = '#006699', mid = "white",
                                    midpoint = 0, limit = c(-1, 1), space = "Lab",
                                    name = "Asymmetry", na.value = "lightgrey")
    }
    p <- p + theme_bw() + coord_fixed() + xlab('Variable 2') + ylab("Variable 1")
    p <- p + scale_x_discrete(position = 'top')
    p <- p + scale_y_discrete(limits = rev(unique(melt_df$Var1)))
    if(significance & any(!is.na(pw_qad$asymmetry_p.values))){
      melt_df_sign <- melt(as.matrix(pw_qad$asymmetry_p.values))
      melt_df_sign$value <- ifelse(melt_df_sign$value < sign.level & !is.na(melt_df_sign$value), '*','')
      melt_df$sign <- melt_df_sign$value
      melt_df$paste_value <- ifelse(is.na(melt_df$value), NA, paste(melt_df$value,melt_df$sign, sep=''))
      p <- p + geom_text(data = melt_df, aes_(x = ~Var2, y = ~Var1, label = ~paste_value), color = "black", size = fontsize, na.rm=TRUE)
    }else{
      p <- p + geom_text(aes_(x = ~Var2, y = ~Var1, label = ~value), color = "black", size = fontsize, na.rm=TRUE)
    }
    p



  }else{
    warning('Select an appropriate variable. Options are c("dependence","mean.dependence","asymmetry")')
  }
}

