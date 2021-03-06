% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/simulators.R
\name{volesSimulator}
\alias{volesSimulator}
\title{Simulates from the voles models of Fasiolo and Wood (2015).}
\usage{
volesSimulator(param, nsim, nObs, extraArgs = NULL, model = "full",
  nBurn = 12 * 10, nSteps = 10, T0 = 0, ...)
}
\arguments{
\item{param}{If model == "full" it must be a vector of length 9 or a matrix of size nsim by 9, if model == "standard"
it must be a vector of length 7 or a matrix of size nsim by 7.
If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
will use the i-th row of param as parameters.}

\item{nsim}{Number of simulations from the model.}

\item{nObs}{Length of each simulated time series in months.}

\item{extraArgs}{Only needed for compatibility with the "synlik" package.}

\item{model}{If model == "full" the function will simulate from the full model described in eq. 2 of Turchin and Ellner (2000),
model == "standard" the function will simulate from the full model described in eq. 6 of Turchin and Ellner (2000).}

\item{nBurn}{Number of initial steps to be discarded before saving the following nObs steps.}

\item{nSteps}{Number of finite differences steps between 2 observations (an interval which corresponds to a month).}

\item{...}{Only needed for compatibility with the "synlik" package.}
}
\value{
A list where [[""voles"]] and [["weasels"]] are two matrices of size nsim by nObs where each row is a trajectory
        simulated from the model.
}
\description{
Simulates from the voles models of Fasiolo and Wood (2015).
}
\examples{
library(volesModel)

########## Simulating from full model
set.seed(76735756)
nSteps <- 50
nObs <- 500
res <- volesSimulator(nObs = nObs, nsim = 100, nBurn = 120, model = "full",
                      param = log(c(r = 4.5, e = 1, g = 0.12, h = 0.1, a = 8, d = 0.04, s = 1.25,
                                    sigmaProc = 1.5, phi = 100)),
                      nSteps = nSteps, T0 = 0)

ii = 1
par(mfrow = c(2, 1))
plot(res[[1]][ii, ], type = 'l', main = "Voles", ylab = "Voles", xlab = "Semester")
plot(res[[2]][ii, ], type = 'l', main = "Weasels", ylab = "Weasels", xlab = "Semester")
ii = ii + 1

########## Simulating from simpler model
# Standard Model (maybe 99 steps is fine)
set.seed(452535626)
res <- volesSimulator(nObs = nObs, nsim = 100, nBurn = 120, model = "standard",
                      param = log(c(r = 3.07, e = 1.58, beta = 1.54, d = 0.6, alpha = 7.6, sigmaProc = 1.5, phi = 100)),
                      nSteps = nSteps, T0 = 0, parall = F)

ii = 1
par(mfrow = c(2, 1))
plot(res[[1]][ii, ], type = 'l', main = "Voles", ylab = "Voles", xlab = "Semester")
plot(res[[2]][ii, ], type = 'l', main = "Weasels", ylab = "Weasels", xlab = "Semester")
ii = ii + 1
}
\author{
Matteo Fasiolo <matteo.fasiolo@gmail.com>
}
\references{
Fasiolo, M and Wood, S. N. (2015). Approximate methods for dynamic ecological models.
             ArXiv XXX. To appear in the Handbook of Approximate Bayesian Computation.

Turchin, P. and S. P. Ellner (2000). Living on the edge of chaos: population
dynamics of fennoscandian voles. Ecology 81 (11), 3099–3116.
}

