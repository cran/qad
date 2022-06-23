#' Plot conditional probabilites
#'
#' @description  Visualizes the conditional probabilities for each strip of the checkerboard copula in the copula setting or in the retransformed sample setting.
#'
#'
#' @param x an object of class qad.
#' @param addSample a logical indicating whether the observations are returned. In the copula setting the mass squares of the empirical copula are added too. (default = FALSE).
#' @param copula a logical indicating whether the plot depicts the conditional probabilities
#' of the empirical checkerboard copula or of the retransformed data setting (default = FALSE).
#' @param density a logical indicating whether the density should be plotted instead of the conditional probabilites (default = FALSE).
#' Only works in the copula setting, i.e. if copula = TRUE.
#' @param margins a logical indicating whether the margin distribution is added in form of a rug plot.
#' @param point.size a numeric specifying the point size of the sample (relevant if addSample = TRUE).
#' @param panel.grid a logical indicating whether the panel grid is plotted. (default = TRUE)
#' @param color a color palette of the viridis package or rainbow. options are c("viridis", "magma", "inferno", "plasma", "cividis", "rainbow")
#' @param rb_values a vector of size 3 with number of values, start value and end value in the rainbow colors space.
#' @param title The text for the title
#' @param x.axis The text for the x-axis
#' @param y.axis The text for the y-axis
#' @param ... some methods for this generic require additional arguments.  None are used in this method.
#'
#' @note The conditional probabilities are constant at squares in the copula setting. If the squares are retransformed in the data setting, the resulting objects are rectangles.
#'
#' @examples
#' ## Example 1
#' n <- 100
#' x <- runif(n, 0, 1)
#' y <- runif(n, 0, 1)
#' sample <- data.frame(x, y)
#'
#' #qad
#' fit <- qad(sample)
#' plot(fit, addSample = TRUE, copula = FALSE)
#'
#' ## Example 2
#' n <- 100
#' x <- runif(n, -1, 1)
#' y <- x^2 + rnorm(n, 0, 0.1)
#' sample <- data.frame(x, y)
#'
#' #qad
#' fit <- qad(sample)
#' plot(fit, addSample = TRUE, copula = TRUE)
#' plot(fit, addSample = TRUE, copula = FALSE)
#'
#' @method plot qad

plot.qad <- function(x,
                     addSample = FALSE,
                     copula = FALSE,
                     density = FALSE,
                     margins = FALSE,
                     title = "",
                     x.axis = "X1",
                     y.axis = "X2",
                     point.size = 0.9,
                     panel.grid = TRUE,
                     color = "plasma",
                     rb_values = c(10, 0.315, 0.15), ...){
  qad_output <- x
  if(inherits(qad_output, 'qad')){
    #Mass matrix of the copula
    mass_matrix <- qad_output$mass_matrix
    if(density){
      density_matrix <- mass_matrix*NROW(mass_matrix)*NCOL(mass_matrix)
      density_matrix <- ifelse(density_matrix == 0, NA, density_matrix)
      data <- as.data.frame(as.table(density_matrix))
      n <- NROW(density_matrix)
      l.title <- "Density"
    }else{
      cond_prob_matrix <- mass_matrix*NROW(mass_matrix)
      cond_prob_matrix <- ifelse(cond_prob_matrix==0, NA, cond_prob_matrix)
      data <- as.data.frame(as.table(cond_prob_matrix))
      n <- NROW(cond_prob_matrix)
      l.title <- "Cond. prob."
    }

    names(data) <- c('x1','x2','value')
    data$x1 <- as.numeric(data$x1)
    data$x2 <- as.numeric(data$x2)
    sample_size <- length(qad_output$data$x1)


    if(copula){
      p <- ggplot()
      p <- p + geom_raster(data=data, aes_(x=~x1/n-1/n/2,y=~x2/n-1/n/2 , fill=~value), alpha=0.95, na.rm = TRUE)
      if(color == "rainbow"){
        p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                      guide = guide_colorbar(title=l.title), na.value=NA)
      }else{
        p <- p + scale_fill_viridis_c(na.value = NA, space = "Lab", name = l.title, option = color, direction = -1)
      }
      p <- p + theme_bw() + xlab(x.axis) + ylab(y.axis)
      if(addSample){
        if(n < 20){
          M <- round(ECBC(qad_output$data$x1, qad_output$data$x2, resolution = length(qad_output$data$x1)),15)
          M <- ifelse(M==0, NA, M)
          m <- NROW(M)

          Mdf <- as.data.frame(as.table(M), dnn =  1:m)
          names(Mdf) <- c('x1','x2','value')
          Mdf$x1 <- as.numeric(Mdf$x1)
          Mdf$x2 <- as.numeric(Mdf$x2)
          Mdf <- Mdf[!is.na(Mdf$value),]
          p <- p + geom_raster(data = Mdf, aes_(x=~x1/m - 1/m/2,y=~x2/m - 1/m/2, fill = ~value), alpha = 0.35, fill = "grey20", na.rm = T)
          p <- p + geom_point(data=qad_output$data, aes(x=rank(x1, ties.method = 'max')/sample_size,
                                                        y=rank(x2, ties.method = 'max')/sample_size), size=point.size)
        }else{
          p <- p + geom_point(data=qad_output$data, aes(x=rank(x1, ties.method = 'max')/sample_size,
                                                        y=rank(x2, ties.method = 'max')/sample_size), size=point.size)
        }

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
      p <- p + ggtitle(title) + xlab(x.axis) + ylab(y.axis)

      return(p)
#Data setting
    } else{
      x1 <- qad_output$data$x1
      x2 <- qad_output$data$x2
      N <- qad_output$resolution

      #New
      grid <- seq(0,1,length.out = N+1)
      x_grid <- quantile(x1, grid, type = 1)
      y_grid <- quantile(x2, grid, type = 1)
      raster_grid1 <- expand.grid(xmin = x_grid[1:N], ymin = y_grid[1:N])
      raster_grid2 <- expand.grid(xmax = x_grid[2:(N+1)], ymax = y_grid[2:(N+1)])
      raster_grid <- dplyr::bind_cols(raster_grid1, raster_grid2)
      raster_grid <- dplyr::bind_cols(raster_grid, data.frame(prob = data$value))


      p <- ggplot()
      p <- p + geom_rect(data=raster_grid, aes_(xmin=~xmin, xmax=~xmax, ymin=~ymin, ymax=~ymax,fill=~prob), alpha=0.95, na.rm = TRUE)
      if(margins){
        p <- p + geom_rug(data = qad_output$data, aes_(x=~x1, y=~x2))
      }
      if(color == "rainbow"){
        p <- p + scale_fill_gradientn(colours = rainbow(rb_values[1], start = rb_values[2], end = rb_values[3]),
                                      guide = guide_colorbar(title= l.title), na.value=NA)
      }else{
        p <- p + scale_fill_viridis_c(na.value = NA, space = "Lab", name = l.title, option = color, direction = -1)
      }
      if(addSample){
          p <- p + geom_point(data=qad_output$data, aes_(x=~x1,y=~x2), size=point.size)
      }
      p <- p + theme_bw() + xlab('x1') + ylab('x2')
      if(!panel.grid){
        p <- p + theme(panel.grid.minor = element_line(colour = NA),
                       panel.grid.major = element_line(colour = NA))
      }
      p <- p + theme(panel.grid = element_line(linetype = 'dashed'))
      p <- p + scale_x_continuous(breaks = function(x) pretty(x, 15))
      p <- p + scale_y_continuous(breaks = function(x) pretty(x, 15))
      p <- p + ggtitle(title) + xlab(x.axis) + ylab(y.axis)
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
#' plot(x,y,pch = 16)
#'
#' mass <- ECBC(x,y,resolution = 10)
#' plot_density(mass, density=TRUE)
#' plot_density(mass, density=FALSE)

plot_density <- function(mass_matrix, density=TRUE, color = "plasma", rb_values = c(10, 0.315, 0.15)){
  #input: matrix with mass distribution
  if(density){
    mass_matrix <- mass_matrix*NROW(mass_matrix)*NCOL(mass_matrix)
  }
  mass_matrix <- ifelse(mass_matrix==0, NA, mass_matrix)

  data <- as.data.frame(as.table(mass_matrix), dnn =  1:NROW(mass_matrix))
  names(data) <- c('x1','x2','value')
  data$x1 <- as.numeric(data$x1)
  data$x2 <- as.numeric(data$x2)

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
#' Options are c("dependence", "max.dependence", "asymmetry").
#' @param fontsize a numeric specifying the font size of the values.
#' @param significance a logical indicating whether significant values with respect to the (adjusted) qad p.values are denoted by a star.
#' @param use_p.adjust a logical indicating if the adjusted p.values are used (default = TRUE).
#' @param sign.level numeric value indicating the significance level.
#' @param scale character indicating whether the heatmap uses a relative or absolute scale. Options are 'rel' or 'abs' (default).
#' @param color Select the color palette. Options are c("plasma" (default), "viridis", "inferno", "magma", "cividis", "rainbow").
#' @param white_font numeric between 0 and 1 denoting the start value for white text font (default = 0.7)
#' @param rb_values a vector of size 3 with number of values, start value and end value in the rainbow colors space (if color = 'rainbow').
#' @param title The text for the title
#' @details If the output of \code{pairwise.qad}() contains p-values, significant values can be highlighted by stars by setting significance=TRUE.
#'
#' @return a heatmap
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

heatmap.qad <- function(pw_qad,
                        select = c('dependence','max.dependence','asymmetry'),
                        fontsize = 4,
                        significance = FALSE,
                        use_p.adjust = TRUE,
                        sign.level = 0.05,
                        scale = "abs",
                        color = "plasma",
                        white_font = 0.7,
                        rb_values = c(10, 0.315, 0.15),
                        title = ""){

  prepare_data_long <- function(matr, matr.p, n_round = 3){
    df_matr <- as.data.frame(as.table(round(as.matrix(matr), n_round)))
    names(df_matr) <- c("Var1", "Var2", "value")

    df_p <- as.data.frame(as.table(as.matrix(matr.p)))
    names(df_p) <- c("Var1", "Var2", "p.value")

    df_p$sign <- with(df_p, ifelse(p.value < sign.level & !is.na(p.value), "*", ""))
    df_matr <- dplyr::left_join(df_matr, df_p, by = c("Var1", "Var2"))
    df_matr$sign_label <- with(df_matr, ifelse(is.na(value), NA, paste0(value,sign)))

    return(df_matr)
  }


  if(color == "rainbow"){
    color_pal <- rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])
  }else{
    color_pal <- viridis(5, option = color, direction = -1)
  }
  limit_plot <- c(0,1)

  #Select dependence measure
  if(select == 'dependence'){
    if(use_p.adjust){
      df_long <- prepare_data_long(pw_qad$q, pw_qad$q_p.values.adjusted)
    }else{
      df_long <- prepare_data_long(pw_qad$q, pw_qad$q_p.values)
    }
    legend_title <- "Dependence: \nq:=q(x,y)"

  }else if(select == 'max.dependence'){
    if(use_p.adjust){
      df_long <- prepare_data_long(pw_qad$max.dependence, pw_qad$max.dependence_p.values.adjusted)
    }else{
      df_long <- prepare_data_long(pw_qad$max.dependence, pw_qad$max.dependence_p.values)
    }
    legend_title <- "Max dependence:=\nmax(q(x,y)+q(y,x))"

  }else if(select == 'asymmetry'){
    if(use_p.adjust){
      df_long <- prepare_data_long(pw_qad$asymmetry, pw_qad$asymmetry_p.values.adjusted)
    }else{
      df_long <- prepare_data_long(pw_qad$asymmetry, pw_qad$asymmetry_p.values)
    }
    limit_plot <- c(-1,1)
    legend_title <- "Asymmetry: \na:=q(x,y)-q(y,x)"
    if(color == "rainbow"){
      color_pal <- c(rev(rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])),rainbow(rb_values[1], start = rb_values[2], end = rb_values[3])[-1])
    }else{
      color_pal <- c(rev(viridis(5, option = color, direction = -1)),viridis(5, option = color, direction = -1)[-1])
    }
  }else{
    stop('Select an appropriate select variable. Options are c("dependence","max.dependence","asymmetry")')
  }

  df_long$text_color <- ifelse(df_long$value > white_font, "high", "low")



  #Plot heatmap
  p <- ggplot(data = df_long, aes_(x = ~Var2, y = ~Var1, fill = ~value))
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

  if(significance){
    p <- p + geom_text(aes_(label = ~sign_label, color = ~text_color), size = fontsize, na.rm=TRUE)
  }else{
    p <- p + geom_text(aes_(label = ~value, color = ~text_color), size = fontsize, na.rm = TRUE)
  }
  p <- p + scale_color_manual(values = c("high" = "white", "low" = "black")) + guides(color = "none")

  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  p <- p + scale_y_discrete(limits = rev(levels(df_long$Var1)))
  if(significance){
    if(all(is.na(df_long$p.value))){
      warning(paste0("No p.values available"))
    }else{
      if(use_p.adjust){
        p <- p + labs(caption = paste0("Signif. codes:"," '*' ","p.value <", sign.level, " (p.value correction method: ", pw_qad$p.adjust.method, ")"))
        message(paste0("p.values are adjusted by using the method ", pw_qad$p.adjust.method))
      }else{
        p <- p + labs(caption = paste0("Signif. codes:"," '*' ","p.value <", sign.level))
      }
    }
  }
  p <- p + ggtitle(title)
  return(p)
}





