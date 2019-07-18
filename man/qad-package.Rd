\name{qad-package}
\alias{qad-package}
\docType{package}
\title{
Quantification of Asymmetric Dependencies
}
\description{
A copula-based measure for quantifying asymmetry in dependence and associations.
}
\details{
\tabular{ll}{
Package: \tab qad\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2018-12-18\cr
}

}
\author{
Florian Griessenberger: \email{florian.griessenberger@sbg.ac.at},\cr
Robert R. Junker: \email{Robert.Junker@sbg.ac.at},\cr
Wolfgang Trutschnig: \email{Wolfgang.Trutschnig@sbg.ac.at}
}
\keyword{ package }

\examples{
##Create data set
# n <- 1000
# x <- rnorm(n,0,2)
# y <- x^2 + rnorm(n)
# sample <- data.frame(x,y)
# plot(sample)

##Function: empirical copula
# eval <- data.frame(x=c(0,0.1,0.2,0.45,1), y=c(0,0.1,0.5,0.23,1))
# emp_c_copula_eval(sample, eval, resolution = 10)
# mass <- emp_c_copula(sample, resolution=10)

##Function: qad()
# qad(sample, resolution = NULL, permutation = FALSE, nperm = 100, DoParallel = TRUE, ncores = NULL)
# help(qad)
# mod <- qad(sample)
# mod <- qad(sample, resolution = NULL, permutation = TRUE, nperm = 10, DoParallel = TRUE)

##Functions: summary(), coef()
# help(summary.qad)

# mod <- qad(sample)
# summary(mod)
# coef(mod)
# coef(mod, select = c('q(x2,x1)','mean.dependence'))

##Function: plot()
# help(plot.qad)

# plot(mod)
# plot(mod, addSample = TRUE, copula = FALSE, margins = TRUE, point.size = 0.7, panel.grid = FALSE)

##Function: cci()
# help(cci)

# n <- 1000
# cci(n, alternative = "one.sided")
# cci(n, alternative = "two.sided")

##Function: predict()
# help(predict.qad)
# values <- c(-2.4,1,0,2.6)
# predict.qad(mod, values = values, conditioned = 'x1')
# predict(mod, values, conditioned = "x1", nr_intervals = 10, pred_plot = TRUE, panel.grid = FALSE)

# values <- c(0.1,0.5)
# predict(mod, values, conditioned = "x2", nr_intervals = 5, copula = TRUE, pred_plot = TRUE)

##Function: pairwise.qad and heatmap.qad
# df <- iris[1:4]
# mod <- pairwise.qad(df)
# heatmap.qad(mod, select = 'dependence')
}
