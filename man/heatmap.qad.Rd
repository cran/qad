% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{heatmap.qad}
\alias{heatmap.qad}
\title{Heatmap of dependence measures}
\usage{
heatmap.qad(
  pw_qad,
  select = c("dependence", "max.dependence", "asymmetry"),
  fontsize = 4,
  significance = FALSE,
  use_p.adjust = TRUE,
  sign.level = 0.05,
  scale = "abs",
  color = "plasma",
  white_font = 0.7,
  rb_values = c(10, 0.315, 0.15),
  title = ""
)
}
\arguments{
\item{pw_qad}{output of the function \code{pairwise.qad}().}

\item{select}{a character indicating which dependence value is plotted.
Options are c("dependence", "max.dependence", "asymmetry").}

\item{fontsize}{a numeric specifying the font size of the values.}

\item{significance}{a logical indicating whether significant values with respect to the (adjusted) qad p.values are denoted by a star.}

\item{use_p.adjust}{a logical indicating if the adjusted p.values are used (default = TRUE).}

\item{sign.level}{numeric value indicating the significance level.}

\item{scale}{character indicating whether the heatmap uses a relative or absolute scale. Options are 'rel' or 'abs' (default).}

\item{color}{Select the color palette. Options are c("plasma" (default), "viridis", "inferno", "magma", "cividis", "rainbow").}

\item{white_font}{numeric between 0 and 1 denoting the start value for white text font (default = 0.7)}

\item{rb_values}{a vector of size 3 with number of values, start value and end value in the rainbow colors space (if color = 'rainbow').}

\item{title}{The text for the title}
}
\value{
a heatmap
}
\description{
The pairwise computed dependence measures (output of the function \code{pairwise.qad()}) are illustrated by a heatmap.
}
\details{
If the output of \code{pairwise.qad}() contains p-values, significant values can be highlighted by stars by setting significance=TRUE.
}
\examples{
n <- 100
x1 <- runif(n, 0, 1)
x2 <- x1^2 + rnorm(n, 0, 0.1)
x3 <- runif(n, 0, 1)
x4 <- x3 - x2 + rnorm(n, 0, 0.1)
sample_df <- data.frame(x1,x2,x3,x4)
#Fit qad
model <- pairwise.qad(sample_df, p.value = FALSE)
heatmap.qad(model, select = "dependence", fontsize = 6)
}
