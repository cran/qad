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
#' @param color a color palette of the viridis package or rainbow. options are c("viridis", "magma", "inferno", "plasma", "cividis", "rainbow")
#' @param rb_values a vector of size 3 with number of values, start value and end value in the rainbow colors space.
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

plot.qad <- function(x, addSample = FALSE, copula = FALSE, density = FALSE, margins = FALSE, point.size = 0.8, panel.grid = TRUE,
                     color = "plasma", rb_values = c(10, 0.315, 0.15), ...){
  qad_output <- x
  if(class(qad_output)=='qad'){
#Copula setting
    if(copula){
      mass_matrix <- qad_output$mass_matrix

      if(density){
        density_matrix <- mass_matrix*NROW(mass_matrix)*NCOL(mass_matrix)
        density_matrix <- ifelse(density_matrix == 0, NA, density_matrix)

        density_matrix <- data.table::as.data.table(density_matrix)
        names(density_matrix) <- as.character(1:NCOL(density_matrix))
        density_matrix$Var1 <- colnames(density_matrix)
        data <- data.table::melt(density_matrix, id.vars = "Var1", variable.name = "Var2")
        data$Var1 <- as.numeric(as.character(data$Var1))
        data$Var2 <- as.numeric(as.character(data$Var2))

        #data <- melt(density_matrix)
        n <- NROW(density_matrix)
      }else{
        cond_prob_matrix <- mass_matrix*NROW(mass_matrix)
        cond_prob_matrix <- ifelse(cond_prob_matrix==0, NA, cond_prob_matrix)

        cond_prob_matrix <- data.table::as.data.table(cond_prob_matrix)
        names(cond_prob_matrix) <- as.character(1:NCOL(cond_prob_matrix))
        cond_prob_matrix$Var1 <- colnames(cond_prob_matrix)
        data <- data.table::melt(cond_prob_matrix, id.vars = "Var1", variable.name = "Var2")
        data$Var1 <- as.numeric(as.character(data$Var1))
        data$Var2 <- as.numeric(as.character(data$Var2))

        #data <- melt(cond_prob_matrix)
        n <- NROW(cond_prob_matrix)
      }

      names(data) <- c('x1','x2','value')
      sample_size <- length(qad_output$data$x1)

      p <- ggplot()
      p <- p + geom_raster(data=data, aes_(x=~x1/n-1/n/2,y=~x2/n-1/n/2 , fill=~value), alpha=0.95, na.rm = TRUE)
      if(density){
        if(color == "rainbow"){
          p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                        guide = guide_colorbar(title='Density'), na.value=NA)
        }else{
          p <- p + scale_fill_viridis_c(na.value = NA, space = "Lab", name = 'Density', option = color, direction = -1)
        }
      }else{
        if(color == "rainbow"){
          p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                        guide = guide_colorbar(title='Cond. prob.'), na.value=NA)
        }else{
          p <- p + scale_fill_viridis_c(na.value = NA, space = "Lab", name = 'Cond. prob.', option = color, direction = -1)
        }
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
      if(color == "rainbow"){
        p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                      guide = guide_colorbar(title='Cond. prob.'), na.value=NA)
      }else{
        p <- p + scale_fill_viridis_c(na.value = NA, space = "Lab", name = 'Cond. prob.', option = color, direction = -1)
      }
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
#' @param color Select the color palette. Options are c("plasma" (default), "viridis", "inferno", "magma", "cividis").
#' @param rb_values a vector of size 3 with number of values, start value and end value in the rainbow colors space.
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

plot_density <- function(mass_matrix, density=TRUE, color = "plasma", rb_values = c(10, 0.315, 0.15)){
  #input: matrix with mass distribution
  if(density){
    mass_matrix <- mass_matrix*NROW(mass_matrix)*NCOL(mass_matrix)
  }
  mass_matrix <- ifelse(mass_matrix==0, NA, mass_matrix)

  mass_matrix <- data.table::as.data.table(mass_matrix)
  names(mass_matrix) <- as.character(1:NCOL(mass_matrix))
  mass_matrix$Var1 <- colnames(mass_matrix)
  data <- data.table::melt(mass_matrix, id.vars = "Var1", variable.name = "Var2")
  data$Var1 <- as.numeric(as.character(data$Var1))
  data$Var2 <- as.numeric(as.character(data$Var2))

  #data <- melt(mass_matrix)
  names(data) <- c('x1','x2','value')
  data <- as.data.frame(data)
  n <- NROW(mass_matrix)

  p <- ggplot(data=data, aes_(x=~x1/n - 1/n/2,y=~x2/n - 1/n/2))
  p <- p + geom_raster(aes_(fill=~value), alpha=0.95, na.rm = TRUE)
  if(density){
    if(color == "rainbow"){
      p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                    guide = guide_colorbar(title='Density'), na.value=NA)
    }else{
      p <- p + scale_fill_viridis_c(guide = guide_colorbar(title = "Density"), na.value = NA, option = color, direction = -1)
    }
    #p <- p + scale_fill_gradient(low="#cde7f0",high='#006699',guide = guide_colorbar(title='Density'), na.value=NA)
  }else{
    if(color == "rainbow"){
      p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                    guide = guide_colorbar(title='Mass'), na.value=NA)
    }else{
      p <- p + scale_fill_viridis_c(guide = guide_colorbar(title = "Mass"), na.value = NA, option = color, direction = -1)
    }
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
#' @param color Select the color palette. Options are c("plasma" (default), "viridis", "inferno", "magma", "cividis").
#' @param rb_values a vector of size 3 with number of values, start value and end value in the rainbow colors space.
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
#' model <- pairwise.qad(sample_df, permutation = FALSE)
#' heatmap.qad(model, select = "dependence", fontsize = 6)

heatmap.qad <- function(pw_qad, select = c('dependence','mean.dependence','asymmetry'), fontsize = 4, significance = FALSE,
                        sign.level = 0.05, scale = "abs", color = "plasma", rb_values = c(10, 0.315, 0.15)){
  #Select dependence measure
  if(select == 'dependence'){
    ma <- as.matrix(pw_qad$q)
    ma <- data.table::as.data.table(ma)
    names(ma) <- as.character(1:NCOL(ma))
    ma$Var1 <- colnames(ma)
    melt_df <- data.table::melt(ma, id.vars = "Var1", variable.name = "Var2")
    melt_df$Var1 <- as.numeric(as.character(melt_df$Var1))
    melt_df$Var2 <- as.numeric(as.character(melt_df$Var2))

    #melt_df <- melt(as.matrix(pw_qad$q))
    melt_df$value <- round(melt_df$value, 2)
    if(significance & any(!is.na(pw_qad$q_p.values))){
      ma_s <- as.matrix(pw_qad$q_p.values)
      ma_s <- data.table::as.data.table(ma_s)
      names(ma_s) <- as.character(1:NCOL(ma_s))
      ma_s$Var1 <- colnames(ma_s)
      melt_df_sign <- data.table::melt(ma_s, id.vars = "Var1", variable.name = "Var2")
      melt_df_sign$Var1 <- as.numeric(as.character(melt_df_sign$Var1))
      melt_df_sign$Var2 <- as.numeric(as.character(melt_df_sign$Var2))

      #melt_df_sign <- melt(as.matrix(pw_qad$q_p.values))
      melt_df_sign$value <- ifelse(melt_df_sign$value < sign.level & !is.na(melt_df_sign$value), '*','')
      melt_df$sign <- melt_df_sign$value
      melt_df$paste_value <- ifelse(is.na(melt_df$value), NA, paste(melt_df$value,melt_df$sign, sep=''))
    }
    limit_plot <- c(0,1)
    legend_title <- "Dependence: \nq:=q(x,y)"
    #color_pal <- c("white",  "#365C8DFF", "#1FA187FF", "#FDE725FF")
    if(color == "rainbow"){
      color_pal <- rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])
    }else{
      color_pal <- viridis(5, option = color, direction = -1)
    }
  }else if(select == 'mean.dependence'){
    ma <- as.matrix(pw_qad$mean.dependence)
    ma <- data.table::as.data.table(ma)
    names(ma) <- as.character(1:NCOL(ma))
    ma$Var1 <- colnames(ma)
    melt_df <- data.table::melt(ma, id.vars = "Var1", variable.name = "Var2")
    melt_df$Var1 <- as.numeric(as.character(melt_df$Var1))
    melt_df$Var2 <- as.numeric(as.character(melt_df$Var2))

    #melt_df <- melt(as.matrix(pw_qad$mean.dependence))
    melt_df$value <- round(melt_df$value, 2)
    if(significance & any(!is.na(pw_qad$mean.dependence_p.values))){
      ma_s <- as.matrix(pw_qad$mean.dependence_p.values)
      ma_s <- data.table::as.data.table(ma_s)
      names(ma_s) <- as.character(1:NCOL(ma_s))
      ma_s$Var1 <- colnames(ma_s)
      melt_df_sign <- data.table::melt(ma_s, id.vars = "Var1", variable.name = "Var2")
      melt_df_sign$Var1 <- as.numeric(as.character(melt_df_sign$Var1))
      melt_df_sign$Var2 <- as.numeric(as.character(melt_df_sign$Var2))

      #melt_df_sign <- melt(as.matrix(pw_qad$mean.dependence_p.values))
      melt_df_sign$value <- ifelse(melt_df_sign$value < sign.level & !is.na(melt_df_sign$value), '*','')
      melt_df$sign <- melt_df_sign$value
      melt_df$paste_value <- ifelse(is.na(melt_df$value), NA, paste(melt_df$value,melt_df$sign, sep=''))
    }
    limit_plot <- c(0,1)
    legend_title <- "Mean dependence: \nmd:=(q(x,y)+q(y,x))/2"
    #color_pal <- c("white",  "#365C8DFF", "#1FA187FF", "#FDE725FF")
    if(color == "rainbow"){
      color_pal <- rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])
    }else{
      color_pal <- viridis(5, option = color, direction = -1)
    }
  }else if(select == 'asymmetry'){
    ma <- as.matrix(pw_qad$asymmetry)
    ma <- data.table::as.data.table(ma)
    names(ma) <- as.character(1:NCOL(ma))
    ma$Var1 <- colnames(ma)
    melt_df <- data.table::melt(ma, id.vars = "Var1", variable.name = "Var2")
    melt_df$Var1 <- as.numeric(as.character(melt_df$Var1))
    melt_df$Var2 <- as.numeric(as.character(melt_df$Var2))

    #melt_df <- melt(as.matrix(pw_qad$asymmetry))
    melt_df$value <- round(melt_df$value, 2)
    if(significance & any(!is.na(pw_qad$asymmetry_p.values))){
      ma_s <- as.matrix(pw_qad$asymmetry_p.values)
      ma_s <- data.table::as.data.table(ma_s)
      names(ma_s) <- as.character(1:NCOL(ma_s))
      ma_s$Var1 <- colnames(ma_s)
      melt_df_sign <- data.table::melt(ma_s, id.vars = "Var1", variable.name = "Var2")
      melt_df_sign$Var1 <- as.numeric(as.character(melt_df_sign$Var1))
      melt_df_sign$Var2 <- as.numeric(as.character(melt_df_sign$Var2))

      #melt_df_sign <- melt(as.matrix(pw_qad$asymmetry_p.values))
      melt_df_sign$value <- ifelse(melt_df_sign$value < sign.level & !is.na(melt_df_sign$value), '*','')
      melt_df$sign <- melt_df_sign$value
      melt_df$paste_value <- ifelse(is.na(melt_df$value), NA, paste(melt_df$value,melt_df$sign, sep=''))
    }
    limit_plot <- c(-1,1)
    legend_title <- "Asymmetry: \na:=q(x,y)-q(y,x)"
    #color_pal <- c("#FDE725FF","#1FA187FF","#365C8DFF", "white",  "#365C8DFF", "#1FA187FF", "#FDE725FF")
    if(color == "rainbow"){
      color_pal <- c(rev(rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])),rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])[-1])
    }else{
      color_pal <- c(rev(viridis(5, option = color, direction = -1)),viridis(5, option = color, direction = -1)[-1])
    }
  }else{
    stop('Select an appropriate select variable. Options are c("dependence","mean.dependence","asymmetry")')
  }

  #Plot heatmap
  p <- ggplot(data = melt_df, aes_(x = ~Var2, y = ~Var1, fill = ~value))
  p <- p + geom_tile(color = 'white')
  if(scale == "rel"){
    p <- p + scale_fill_gradientn(colors = color_pal, na.value = "lightgrey", space = "Lab",
                                  name = paste(legend_title))
  }else{
    p <- p + scale_fill_gradientn(colors = color_pal,
                                  limits = limit_plot, na.value = "lightgrey", space = "Lab",
                                  name = paste(legend_title))
  }
  p <- p + theme_bw() + coord_fixed() + xlab('Variable 2 (Y)') + ylab("Variable 1 (X)")
  p <- p + scale_x_discrete(position = 'top')
  p <- p + scale_y_discrete(limits = rev(unique(melt_df$Var1)))

  if(significance & any(!is.na(pw_qad$q_p.values))){
    p <- p + geom_text(data = melt_df, aes_(x = ~Var2, y = ~Var1, label = ~paste_value), color = "black", size = fontsize, na.rm=TRUE)
  }else{
    p <- p + geom_text(aes_(x = ~Var2, y = ~Var1, label = ~value), color = "black", size = fontsize, na.rm = TRUE)
  }
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = -0.5))
  p
  return(p)
}

