#.markov_kernel calculates the markov kernel of a checkerboard copula
.markov_kernel <- function(x, mass) {
  #x = point at which value the markov kernel is evaluated
  resolution <- NROW(mass)
  if (resolution == 1)
    mass <- as.matrix(mass)
  #Select the appropriate row according to x
  index <- max(which(0:resolution / resolution < x))
  kernel <-
    data.table(
      x = rep(x, resolution + 1),
      y = seq(0, 1, length.out = resolution + 1),
      Markov_kernel = cumsum(c(0, mass[index,] * resolution ^ 2 * (1 / resolution)))
    )
  return(kernel)
}

#.D1_checkerboard_strip calculates the D1 distance between a piecewise linear kernel and Pi for one strip
.zeta1_checkerboard_strip <- function(kernel_values_x){
  #input: data.frame of 3 columns
  #       constant x values
  #       y values defining the checkerboard grid
  #       Markov_kernel defining the Markov kernel at x for [0,y]
  resolution <- NROW(kernel_values_x) - 1
  N <- resolution + 1
  kernel_values_x$Markov_kernel_Pi <- kernel_values_x$y
  df_kernel <- kernel_values_x

  #Calculate intersections of the piecewise linear functions with the kernel of PI
  k <- (df_kernel$Markov_kernel[-1] - df_kernel$Markov_kernel[-N]) / (df_kernel$y[-1] - df_kernel$y[-N]) #slope
  d <- df_kernel$Markov_kernel[-1] - k * df_kernel$y[-1]   #intercept
  intersection_points <- round(d / (1 - k), 14)
  #Remove intersection points outside the relevant interval
  intersection_points <- na.omit(
    ifelse(intersection_points < round(df_kernel$y[-1], 14) &
             intersection_points > round(df_kernel$y[-N], 14),
           intersection_points,NA)
  )
  n <- length(intersection_points)

  #Append the additional values for exact caluclation of the integral
  if(n >0){
    df_kernel <- rbind(df_kernel, data.table(x=rep(df_kernel$x[1],n),
                                             y=intersection_points,
                                             Markov_kernel = intersection_points,
                                             Markov_kernel_Pi = intersection_points))
  }
  df_kernel <- df_kernel[order(df_kernel$y),]
  df_kernel$difference <- abs(df_kernel$Markov_kernel - df_kernel$Markov_kernel_Pi)

  #Calculate the integral
  N <- NROW(df_kernel)
  zeta1_x <- 3*sum(0.5*(df_kernel$y[-1] - df_kernel$y[-N]) * (df_kernel$difference[-1] + df_kernel$difference[-N]))
  return(zeta1_x)
}

#.zeta1 calculates the dependence measure zeta_1
.zeta1 <- function(X, resolution = NA) {
  x <- as.numeric(data.frame(X)[, 1])
  y <- as.numeric(data.frame(X)[, 2])

  # Calculate the default resolution
  if (is.na(resolution)) {
    unique_x <- length(unique(x))
    unique_y <- length(unique(y))
    sample_size <- min(unique_x, unique_y)
    resolution <- floor(sample_size ^ 0.5)
  } else{
    resolution <- floor(resolution)
  }

  # Calculate the empirical checkerboard mass
  checkerboard_mass <- emp_c_copula(X = X , resolution = resolution)
  mass.xy <- checkerboard_mass
  mass.yx <- t(checkerboard_mass)

  # Calculate the markov kernels of the empirical checkerboard copula A
  x.values <- seq(0, 1 - 1 / resolution, length.out = resolution) +
    (1 / (2 * resolution))
  # x.values = middle points of the grid on the x-axis
  kernel_values <- lapply(x.values, .markov_kernel, mass.xy)
  zeta1_values_x <- unlist(lapply(kernel_values, .zeta1_checkerboard_strip))
  zeta1 <- sum(zeta1_values_x)/resolution

  # Calculate the markov kernels of the empirical checkerboard copula A^t
  kernel_values <- lapply(x.values, .markov_kernel, mass.yx)
  zeta1_values_x <- unlist(lapply(kernel_values, .zeta1_checkerboard_strip))
  zeta1_transposed <- sum(zeta1_values_x)/resolution

  return(list(zeta1=zeta1, zeta1.t=zeta1_transposed, resolution.checkerboard=resolution, mass_matrix=mass.xy))
}


#Predict function to predict y values for defined x values wrt pseudobservations
.predict_qad_pseudoobservations <- function(values, conditioned = 'x1', qad_output, nr_intervals = NULL, prediction_interval = NULL) {
  x <- unique(values)      #delete duplicate entries
  if(conditioned == 'x2'){
    mass_matrix <- t(qad_output$mass_matrix)    #mass matrix
    resolution <- qad_output$resolution
    mass_cond <- resolution*mass_matrix      #conditional mass matrix
  }else{
    mass_matrix <- qad_output$mass_matrix    #mass matrix
    resolution <- qad_output$resolution
    mass_cond <- resolution*mass_matrix      #conditional mass matrix
  }

  #If prediction_interval is not null consider the prediction_interval
  if(!is.null(prediction_interval)){
    #Calculate conditional mass distribution for new intervall lengths (new y resolution)
    cumulated_mass_cond <- t(apply(mass_cond, 1, cumsum))
    cumulated_mass_cond <- cbind(0,cumulated_mass_cond)
    grid <- seq(0,1,length.out = resolution+1)
    #Calculate the index of each x
    index <- sapply(x, function(x){return(min(which(1:resolution/resolution >= x)))})
    #Calculate the conditional distribution function
    prob_functions <- apply(cumulated_mass_cond, 1, function(x) stats::approxfun(grid, x))
    cond_prob <- rep(NA, length(x))
    for(i in 1:length(index)){
      cond_prob[i] <- prob_functions[[index[i]]](prediction_interval[2]) - prob_functions[[index[i]]](prediction_interval[1])
    }
    df_predict <- data.frame(prob=cond_prob)
    names(df_predict) <- c(paste('[',prediction_interval[1],',',prediction_interval[2],']', sep=''))
    rownames(df_predict) <- x

  #If there is no prediction interval consider the other cases
  }else{
    if(is.null(nr_intervals)) {
      #Calculate conditional mass distribution for intervalls with length 1/resolution
      df_predict <- data.frame(matrix(NA, nrow = length(x), ncol = resolution))
      index <- sapply(x, function(x){return(min(which(1:resolution/resolution >= x)))})
      df_predict[1:length(x), ] <- mass_cond[index, ]
      rownames(df_predict) <- x
      colnames(df_predict) <- paste('(', round(0:(resolution - 1) / resolution, 2), ',', round(1:resolution / resolution, 2), ']', sep ='')
    } else{

      #Calculate conditional mass distribution for new intervall lengths (new y resolution)
      mass_cond_new <- matrix(NA, nrow = resolution, ncol = nr_intervals * resolution)
      for(i in 1:NROW(mass_cond)){
        mass_cond_new[i,] <- as.vector(sapply(mass_cond[i,],function(x) rep(x/nr_intervals,nr_intervals)))
      }
      mass_cond <- matrix(NA, nrow=resolution,ncol=nr_intervals)
      indices <- matrix(1:(nr_intervals*resolution), nrow=nr_intervals, byrow = TRUE)
      for(i in 1:nr_intervals){
        mass_cond[,i] <- apply(mass_cond_new[,as.vector(indices[i,])],1,sum)
      }
      #Output data.frame (rows: x values; columns: intervalls with probabilites)
      df_predict <- data.frame(matrix(NA, nrow=length(x), ncol=nr_intervals))
      #Select for each x value the right strip
      index <- sapply(x,function(x){return(min(which(1:resolution/resolution>=x)))})
      df_predict[1:length(x),] <- mass_cond[index,]
      rownames(df_predict) <- x
      colnames(df_predict) <- paste('(',round(0:(nr_intervals-1)/nr_intervals,2),',',round(1:nr_intervals/nr_intervals,2),']',sep='')
    }
  }

  return(df_predict)
}



#Distribution function of qad (right continuous) P(qad < x)
.ppqad <- function(q, n){
  data <- mcData_independence$q
  data_values <- data$qad_values
  data_param <- data$param

  if(n < 4){
    return(rep(1,length(q)))
  }else if(n < 1000){
    grid <- seq(data_values[[as.character(n)]]$min_grid, data_values[[as.character(n)]]$max_grid, length.out = data_values[[as.character(n)]]$grid_length)
    result <- approxfun(x = grid,
                        y = data_values[[as.character(n)]]$q_n_approx_values_l,
                        yleft = 0, yright = 1)(q)
    return(result)
  }else{
    #Sort is included, otherwise there are wrong quantiles for very big n
    result <- sort(as.numeric(apply(data_param, 1, function(p) p[1]*n^p[2] + p[3])))
    return(approxfun(result, seq(0, 1, length.out = NROW(data_param)),
                     yleft = 0, yright = 1)(q))
  }
}


#Distribution function of mean qad (right continuous) P(meanqad < x)
.ppmqad <- function(q, n){
  data <- mcData_independence$mq
  data_values <- data$qad_values
  data_param <- data$param

  if(n < 4){
    return(rep(1,length(q)))
  }else if(n < 1000){
    grid <- seq(data_values[[as.character(n)]]$min_grid, data_values[[as.character(n)]]$max_grid, length.out = data_values[[as.character(n)]]$grid_length)
    result <- approxfun(x = grid,
                        y = data_values[[as.character(n)]]$q_n_approx_values_l,
                        yleft = 0, yright = 1)(q)
    return(result)
  }else{
    #Sort is included, otherwise there are wrong quantiles for very big n
    result <- sort(as.numeric(apply(data_param, 1, function(p) p[1]*n^p[2] + p[3])))
    return(approxfun(result, seq(0, 1, length.out = NROW(data_param)),
                     yleft = 0, yright = 1)(q))
  }
}


