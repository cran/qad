#Utility function for the computation of the empirical checkerboard copula
#Computes the number of ties
.r_range = function(X) {
  r = rep(0, length(X))
  for (x in X) {
    r[x] = r[x] + 1
  }
  return(r)
}

#Compute the local kernel integral of two checkerboard copulas
.r_local_kernel_integral_AB <- function(A,B,i,j,yA_sum1, yB_sum1, N){
  yA_sum2 <- N*A[i,j] + yA_sum1
  yB_sum2 <- N*B[i,j] + yB_sum1

  #Separate two cases: No intersection and intersection
  if((yA_sum2-yB_sum2)*(yA_sum1-yB_sum1) >= 0){
    I1 <- (yA_sum2-yA_sum1)/(2*N) + yA_sum1/N
    I2 <- (yB_sum2-yB_sum1)/(2*N) + yB_sum1/N
    yD <- abs(I1-I2)
  }else{
    x_intersect <- (yA_sum1-yB_sum1)/(((yB_sum2-yB_sum1)-(yA_sum2-yA_sum1))/(1/N))
    y_intersect <- yA_sum1 + N*(yA_sum2-yA_sum1)*x_intersect
    I11 <- ((y_intersect-yA_sum1)*x_intersect)/2 + yA_sum1*x_intersect
    I12 <- ((y_intersect-yB_sum1)*x_intersect)/2 + yB_sum1*x_intersect
    I21 <- ((yA_sum2 - y_intersect)*(1/N-x_intersect))/2 + y_intersect*(1/N-x_intersect)
    I22 <- ((yB_sum2 - y_intersect)*(1/N-x_intersect))/2 + y_intersect*(1/N-x_intersect)
    yD <- abs(I11-I12)+abs(I21-I22)
  }
  return(c(yA_sum = unname(yA_sum2), yB_sum = unname(yB_sum2), yD = unname(yD)))
}



#.markov_kernel calculates the markov kernel of a checkerboard copula
.r_markov_kernel <- function(x, CB) {
  N <- NROW(CB)
  grid <- seq(0,1,length.out = N+1)
  index <- as.numeric(cut(x, grid,labels = 1:N, right = T, include.lowest = T))
  MK <- sapply(index, function(x) approxfun(grid, N*cumsum(c(0,CB[x,])), method = "linear", yleft = 0, yright = 1))
  names(MK) <- x
  return(MK)
}



#Computes the directed dependence measures
.r_zeta1 <- function(X,Y, resolution){
  N <- resolution
  Pi <- matrix(1/N^2, nrow = N, ncol = N)
  An <- ECBC(X,Y,resolution = N)
  An_t <- t(An)

  return(list(zeta1 = 3*D1(An,Pi),
              zeta1.t = 3*D1(An_t,Pi),
              resolution = resolution,
              mass_matrix = An))
}


.predict_qad_copula <- function(values,
                                conditioned = 'x1',
                                qad_output,
                                nr_intervals = NULL,
                                prediction_interval = NULL) {

  if(conditioned == "x2"){
    mass_matrix <- t(qad_output$mass_matrix)    #mass matrix
    resolution <- qad_output$resolution
    mass_cond <- resolution*mass_matrix      #conditional mass matrix
  }else if(conditioned == "x1"){
    mass_matrix <- as.matrix(qad_output$mass_matrix)    #mass matrix
    resolution <- qad_output$resolution
    mass_cond <- resolution*mass_matrix      #conditional mass matrix
  }

  #Compute markov kernel --> MK_f
  MK <- cbind(0,t(apply(mass_cond, 1, cumsum)))
  grid <- seq(0,1,length.out = resolution + 1)
  MK_f <- apply(MK, 1, function(x) stats::approxfun(grid, x, method = "linear", yleft = 0, yright = 1))


  #Define intervals
  N <- resolution
  if(is.null(nr_intervals) & is.null(prediction_interval)){
    lb <- seq(0,1-1/N, length.out = N)
    ub <- seq(1/N, 1, length.out = N)
    K <- N
  }else if(!is.null(nr_intervals)){
    lb <- seq(0,1-1/nr_intervals, length.out = nr_intervals)
    ub <- seq(1/nr_intervals, 1, length.out = nr_intervals)
    K <- nr_intervals
  }else if(!is.null(prediction_interval)){
    lb <- prediction_interval[1]
    ub <- prediction_interval[2]
    K <- 1
  }

  #Define index, in which interval the value lies
  index <- as.numeric(cut(values, seq(0,1,length.out = N+1), include.lowest = T, labels = 1:N))
  res <- data.frame(matrix(as.numeric(NA), nrow = length(lb), ncol = length(values)+3))
  res[,1] <- paste0("I",1:NROW(res))
  res[,2] <- lb
  res[,3] <- ub
  for(i in 1:length(values)){
    if(!is.na(index[i])){
      res[,i+3] <- MK_f[[index[i]]](ub)-MK_f[[index[i]]](lb)
    }
  }
  res <- data.frame(res)
  colnames(res) <- c("Interval", "lowerBound", "upperBound", paste0(conditioned,"=",values))

  return(list(pred_prob = res, indices = index))
}



