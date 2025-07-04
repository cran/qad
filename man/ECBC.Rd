% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECBC.R
\name{ECBC}
\alias{ECBC}
\alias{ECBC.eval}
\title{Calculate empirical checkerboard copula}
\usage{
ECBC(X, Y, resolution = NULL)

ECBC.eval(CB, eval.points)
}
\arguments{
\item{X}{Numeric vector of values in the first coordinate}

\item{Y}{Numeric vector of values in the second coordinate}

\item{resolution}{Integer indicating the resolution of the checkerboard aggregation, i.e.
the number of vertical/horizontal strips of the checkerboard copula. (default = NULL, sets resolution = floor(sqrt(sample size)))}

\item{CB}{A numeric mass matrix of a checkerboard copula (ECBC)}

\item{eval.points}{A numeric matrix or data.frame indicating the eval.points (x,y)}
}
\value{
\code{ECBC()} returns a matrix with the mass distribution of the empirical
(checkerboard) copula.
}
\description{
The function \code{ECBC()} computes the mass distribution of the
empirical (checkerboard) copula, given a bi-variate sample.
If resolution equals sample size, the bi-linearly extended empirical copula is returned.
Note, if there are ties in the sample an adjusted empirical copula is calculated. The function
\code{ECBC.eval()} evaluates the checkerboard copula at given points.
}
\details{
If the observations are drawn from a continuous distribution (no ties in the sample),
the function \code{ECBC()} returns the commonly used empirical checkerboard copula.
If there are ties in the sample, the empirical copula is adjusted and calculated in the following way: \cr
Let \eqn{(u_i,v_i) := (F_n(x_i),G_n(y_i))} be the pseudo-observations for \eqn{i \in \{1,\ldots,n\}} and \eqn{(u_1',v_1'),\ldots, (u_m',v_m')} the distinct pairs of pseudo-observations with m leq n. Moreover set \eqn{S_1:=\{0, u_1, \ldots, u_{m_1}\}} and \eqn{S_2:=\{0, v_1,\ldots, v_{m_2}\}} and define the quantities \eqn{t_i,r_i} and
\eqn{s_i} for \eqn{i=1,\ldots, m} by
\deqn{t_i := \sum_{j=1}^n 1_{(u_i',v_i')}(u_j,v_j)}
\deqn{r_i := \sum_{j=1}^n 1_{u_i}(u_j)}
\deqn{s_i := \sum_{j=1}^n 1_{v_i}(v_j)}
where 1 defines the indicator function.
Define the empirical subcopula \eqn{A'_n: S_1 x S_2 \to \{0,1/n, \ldots, (n-1)/n,1\}} by
\deqn{A'_n(s_1,s_2)= 1/n  \sum_{i=1}^m t_i * 1_{[0,s_1] x [0,s_2]} (u_i', v_i')=1/n \sum_{i=1}^n 1_{[0,s_1] x [0,s_2]} (u_i, v_i)}
for all \eqn{s_1 \in S_1} and \eqn{s_2 in S_2}. \cr
We extend the subcopula \eqn{A'_n} to a copula by defining the transformations \eqn{w_i:[0,1]^2 \to [u_i'-r_i/n,u_i'] x [v_i'-s_i/n,v_i']} by
\deqn{w_i(x,y)=(u_i'-r_i/n+r_i*x/n, v_i'-s_i/n + s_iy/n)}
and set the measure of the empirical copula \eqn{mu_{A_n}^B := 1/n \sum_{i=1}^m t_i mu_B^{w_i}}, where B denotes the product copula.
}
\examples{
##Generate data drawn from the product copula and compute the empirical (checkerboard) copula
n <- 100
x <- runif(n, 0, 1)
y <- runif(n, 0, 1)
mass <- ECBC(x,y, resolution = 10)
plot_density(mass)
mass <- ECBC(x,y, resolution = n)
plot_density(mass)

## Compute empirical checkerboard copula of a sample with ties and plot density
n <- 100
x <- sample(runif(n, -1, 1), n, replace=TRUE)
y <- x^2 + rnorm(n, 0, 1)
mass <- ECBC(x,y)
plot_density(mass)
}
\references{
Deheuvels, P. (1979). La fonction de dependance empirique et ses proprietas: un test non parametrique d'independance, Acad. Roy. Belg. Bull. Cl. Sci., 5th Ser. 65, 274-292.

Li, X., Mikusinski, P. and Taylor, M.D. (1998). Strong approximation of copulas, Journal of Mathematical Analysis and Applications, 255, 608-623.

Genest, C., Neshlehova J.G. and Remillard, B. (2014). On the empirical multilinear copula process for count data. Bernoulli, 20 (3), 1344-1371.

Junker, R.R., Griessenberger, F. and Trutschnig, W. (2021). Estimating scale-invariant directed dependence of bivariate distributions. Computational Statistics and Data Analysis, 153, 107058.
}
