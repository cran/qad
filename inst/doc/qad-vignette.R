## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "70%", 
  echo = FALSE,
  warning = TRUE
)
library(ggplot2)
library(qad)


## ----eval = FALSE, echo = TRUE, class.source = "fold-hide"--------------------
# R <- 100000
# X <- sample(1:6, R, replace = T)
# Y <- sample(1:6, R, replace = T)
# 
# #Probability that Y = 1
# mean(ifelse(Y == 1, 1, 0))
# 
# #Probability that Y = 1 if we know X
# probs <- rep(0,6)
# for(i in 1:6){
#   probs[i] <- sum(ifelse(Y == 1 & X == i, 1, 0))/sum(X == i)
# }
# #Probability that Y = 1 if we know the value of X
# probs

## ----plot1, echo=F, out.width = "70%", fig.align='center', fig.width=6, fig.height=4.5, fig.cap = "Sample of size n=40 drawn from the model $Y=X^2+\\varepsilon$."----
set.seed(5)
n <- 40
X <- runif(n,-1,1)
Y <- X^2 + rnorm(n, 0, 0.05)
df <- data.frame(X,Y)
f_col <- "grey70"

f <- function(x) return(x^2)
p <- ggplot()
p <- p + geom_point(data = df, aes(x = X, y = Y))
p <- p + geom_function(fun = f, color = "blue", alpha = 0.4, linewidth = 1.2, n = 2000, xlim = c(-1,1))
p <- p + geom_path(aes(x=c(-3/4,-3/4,-1.1), y= c(-0.1, f(-3/4), f(-3/4))), color = f_col, linewidth = 1.1, 
                   arrow = arrow(length = unit(0.1, "inches"), ends = "both"))
p <- p + geom_path(aes(x=c(3/4,3/4,-1.1), y= c(-0.1,f(3/4),  f(3/4))), color = f_col, linewidth = 1.1, 
                   arrow = arrow(length = unit(0.1, "inches"), ends = "both"))
p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
p <- p + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
p <- p + ggtitle("Parabola") + xlab("X") + ylab("Y")
print(p)

## -----------------------------------------------------------------------------
myscatterplot <- function(x,y, pseudo = FALSE){
  df <- data.frame(x,y)
  if(pseudo){
    df$x <- rank(df$x)/length(df$x)
    df$y <- rank(df$y)/length(df$y)
  }
  p <- ggplot(df, aes(x=x, y=y))
  p <- p + geom_point(color = "black", alpha = 0.9, size = 1.2) #ff78b4
  p <- p + theme_bw()
  p <- p + scale_x_continuous(breaks = function(x) pretty(x,8))
  p <- p + scale_y_continuous(breaks = function(x) pretty(x,8))
  p <- p + theme(panel.grid = element_blank())
  return(p)
}
mydistributionplot <- function(fit){
  data <- fit$mass_matrix*fit$resolution
  res <- fit$resolution
  grid <- seq(0, 1, length.out = fit$resolution + 1)

  distr <- matrix(0, nrow = fit$resolution + 1, ncol = fit$resolution)
  distr[-1,] <- apply(data, 1, cumsum)
  distr <- data.frame(distr)
  names(distr) <- c(paste("Strip",1:res))

  distr0 <- distr
  distr0$x <- grid
  df0 <- reshape2::melt(distr0, variable.name = "Kernel", id.vars = c("x"))

  #df_approx
  R <- 1000
  ngrid <- seq(0, 1, length.out = R)
  distr <- data.frame(apply(distr, 2, function(x) return(approx(x = grid, y = x, xout = ngrid)$y)) )
  distr$x <- ngrid
  names(distr) <- c(paste("Strip",1:res),"x")

  df <- reshape2::melt(distr, variable.name = "Kernel", id.vars = c("x"))
  df$min <- pmin(df$value, df$x)
  df$max <- pmax(df$value, df$x)

  p <- ggplot()
  p <- p + geom_line(data = df0, aes(x = x, y = value, color = Kernel), linewidth = 1.03)
  p <- p + geom_line(data = df0, aes(x = x, y = x), linewidth = 1.05, linetype = "dashed")
  p <- p + geom_line(data = subset(df0, df0$Kernel == "Strip 1"), aes(x = x, y = value), color = 'magenta', linewidth = 1.05)

  p <- p + geom_ribbon(data = subset(df, df$Kernel == "Strip 1"), aes(x = x, ymin = min, ymax = max),fill = "magenta", alpha = 0.3)
  p <- p + scale_x_continuous(breaks = function(x) pretty(x, 8))
  p <- p + scale_y_continuous(breaks = function(x) pretty(x, 8))
  p <- p + scale_color_manual(values = rainbow(fit$resolution, start = 0.9, end = 0.4),
                              guide = guide_legend(title.position = "top", title.hjust = 0.5))
  p <- p + theme_bw() + labs(color = "Conditional distribution function for")
  p <- p + theme(legend.position = "bottom", panel.grid.minor = element_line(linetype = "dashed"),
                 panel.grid.major = element_blank())
  return(p)
}
#_____________________________________________________________________________________________

#Example 01 - less noise
x <- df$X
y <- df$Y
coef <- coef(qad(x,y, print = FALSE))


p1 <- myscatterplot(x,y) + xlab("X") + ylab("Y") + ggtitle(paste0("Sample of size n=", n))



fit0 <- qad(x,y,resolution = n, print = FALSE)
p2 <- plot(fit0, addSample = T, copula = T, panel.grid = F, density = T, color = "rainbow", rb_values = c(1,0.75,0.4))
p2 <- p2 + xlab("U:=F(X)~U(0,1)") + ylab("V:=G(Y)~U(0,1)")


p3 <- plot(fit0, addSample = T, copula = T, panel.grid = F, density = T, color = "rainbow", rb_values = c(1,0.75,0.4))
p3 <- p3 + ggtitle("Empirical copula + checkerboard grid") + xlab("U:=F(X)~U(0,1)") + ylab("V:=G(Y)~U(0,1)")
p3 <- p3 + geom_hline(yintercept = seq(0,1,1/floor(sqrt(length(x)))), linetype = "dashed", color = "red")
p3 <- p3 + geom_vline(xintercept = seq(0,1,1/floor(sqrt(length(x)))), linetype = "dashed", color = "red")


fit1 <- qad(X,Y, print = F)
p4 <- plot(fit1, addSample = T, copula = T, panel.grid = F, density = T, color = "rainbow", rb_values = c(30,0.65,0.17)) +
  ggtitle("Empirical checkerboard copula") + xlab("U:=F(X)~ U(0,1)") + ylab("V:=G(Y) ~ U(0,1)")  +
  geom_hline(yintercept = seq(0,1,1/floor(sqrt(length(X)))), linetype = "dashed", color = "red") +
  geom_vline(xintercept = seq(0,1,1/floor(sqrt(length(X)))), linetype = "dashed", color = "red")


p5 <- mydistributionplot(fit1) +
  ggtitle(label = "vertical strips",
          subtitle = paste0("q_n(X,Y) ~", round(coef[1],3))) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + xlab("t") + ylab("P(V <= t | U in I_i)") +
  theme(legend.position = "bottom", axis.title.x = element_blank())


fit2 <- qad(Y,X, print = F)
p6 <- mydistributionplot(fit2) +
  ggtitle(label = "horizontal strips",
          subtitle = paste0("q_n(Y,X) ~", round(coef[2],3))) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("t") + ylab("P(U <= t | V in I_i)") +
  theme(legend.position = "bottom", axis.title.x = element_blank())


## ----fig.width=8, fig.height=7, fig.cap="Sample of size n=40", fig.align="center"----
p1

## ----fig.width=8, fig.height=7, fig.cap="Empirical copula and normalized ranks (points); notice that the masses are uniform over the squares and that,  by construction of the empirical copula, the upper right corner of the squares are the normalized ranks", fig.align="center"----
p2

## ----fig.width=8, fig.height=7, fig.cap="Empirical copula (left panel) and checkerboard aggregation with resolution $N=6$ (right panel)", fig.align="center"----
p3
p4

## ----fig.width=8, fig.height=7, fig.cap="Distance between the conditional distriubtion functions of the checkerboard copula and the product copula, representing independence, for vertical strips (left panel) and horizontal strips (right panel)."----
p5
p6

## ----eval = F, echo = T-------------------------------------------------------
# install.packages("qad")
# library(qad)

## ----echo = T, eval = F-------------------------------------------------------
# help("qad")
# help("qad-package")

## ----example1, echo = TRUE, fig.width=8, fig.height=6, fig.align="center"-----
set.seed(1)

## Step 1: Generate sample 
n <- 100
#Underlying Model Y = sin(X) + small.error
X <- runif(n, -10, 10)
Y <- sin(X) + rnorm(n, 0, 0.1)
#Plot the sample 
plot(X,Y, pch = 16)

#Compute the dependence measure q_n (and the additional p-values obtained by testing for q=0 and a=0)
fit <- qad(X,Y, p.value = T)

## ----example2, echo = TRUE, fig.width=8, fig.height=6, fig.align="center"-----
set.seed(1)

## Step 1: Generate sample 
n <- 100
#Underlying Model Y = sin(X) + error
X <- rnorm(n, 0, 10)
Y <- rnorm(n, 0, 20)
#Plot the sample 
plot(X,Y, pch = 16)

#Compute the dependence measure q
fit <- qad(X,Y, p.value = T)

## ----echo=T-------------------------------------------------------------------
set.seed(11)
x <- runif(4)
y <- runif(4)
qad(x,y)

## ----echo=T, eval=T, results="hide"-------------------------------------------
#Generate sample
set.seed(1)
n <- 250
y <- runif(n, -10, 10)
x <- y^2 + rnorm(n, 0, 6)
df <- data.frame(x=x,y=y)

#Compute qad
fit <- qad(df, print = FALSE)

#Predict the values of Y given X=0 and Y given X=65)
pred <- predict.qad(fit, values=c(0,65), conditioned=c("x1"), pred_plot = FALSE)
#Output as data.frame
pred$prediction
#Output as plot
pp <- pred$plot + theme_classic() + theme(plot.title = element_blank(), legend.position = c(0.9,0.5))
#pp

## ----echo=F, eval=T, results="hold"-------------------------------------------
#Generate sample
pred$prediction

## ----plotforecast, fig.align="center", fig.width=8, fig.height=6, fig.cap="Sample of the data (left panel) and prediction probabilities of qad (right panel)"----
p <- ggplot(df, aes(x=x,y=y))
p <- p + geom_point()
p <- p + theme_bw()
p
pp

## ----echo = T, fig.align="center", fig.width=8, fig.height=6, fig.cap="Empirical checkerboard copula (big squares) and empirical copula (grey rectangles) together with the normalized ranks of a sample of size n=30 containing ties. Note, that the empirical copula may not longer consist of squares with equal length due to the ties."----
set.seed(1)
n <- 30
x <- sample(-10:10, n, replace = T)
y <- x^2 
fit <- qad(x,y, print = F)
plot(fit, addSample = T, copula = T, density = T, point.size = 1.1)

## ----echo = T-----------------------------------------------------------------
qad(x,y, print = T)


