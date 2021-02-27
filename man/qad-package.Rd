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
Version: \tab 1.0.0\cr
Date: \tab 2021-02-26\cr
}

}
\author{
Florian Griessenberger: \email{florian.griessenberger@sbg.ac.at},\cr
Robert R. Junker: \email{Robert.Junker@sbg.ac.at},\cr
Valentin Petztel,\cr
Wolfgang Trutschnig: \email{Wolfgang.Trutschnig@sbg.ac.at}
}
\keyword{ package }

\examples{
# #Create data set
# n <- 100
# x <- rnorm(n,0,2)
# y <- x^2 + rnorm(n)
# sample <- data.frame(x,y)
# plot(sample, pch = 16)
#
#
# #Function: qad()
# qad(sample, p.value = TRUE, p.value_asymmetry = TRUE)
# fit <- qad(sample)
# plot(fit)
# plot(fit, copula = TRUE)
# plot(fit, copula = TRUE, addSample = T)
#
# #Functions: summary(), coef()
# summary(fit)
# coef(fit)
#
# #Function: predict()
# values <- c(-2.4,1,0,2.6)
# predict.qad(fit, values = values, conditioned = 'x1')
# predict(fit, values, conditioned = "x1", pred_plot = TRUE, panel.grid = FALSE)
#
# values <- c(0.1,0.5)
# predict(fit, values, conditioned = "x2", copula = TRUE, pred_plot = TRUE)
#
# #Function: pairwise.qad and heatmap.qad
# df <- iris[1:4]
# mod <- pairwise.qad(df)
# heatmap.qad(mod, select = 'dependence')
}
