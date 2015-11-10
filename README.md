This package contains the voles-weasels prey-predator model described in [Fasiolo and Wood, 2015](http://arxiv.org/abs/1511.02644), which
is a modified version of that originally proposed by Turchin and Ellner, 2000.
The documentation is fairly sparse, but most of the information is contained in the paper.

The function `volesSimulator` can be used to simulate data from the model, see `?volesSimulator` for 
informations about its arguments. The following code simulates data from the model

```R
library(volesModel)

########## Simulating from full model 
set.seed(76735756)
nSteps <- 50
nObs <- 500
res <- volesSimulator(nObs = nObs, nsim = 100, nBurn = 120, model = "full", 
                      param = log(c(r = 4.5, e = 1, g = 0.12, h = 0.1,
                                    a = 8, d = 0.04, s = 1.25, 
                                    sigmaProc = 1.5, phi = 100)), 
                      nSteps = nSteps, T0 = 0)

ii = 1
par(mfrow = c(2, 1))
plot(res[[1]][ii, ], type = 'l', main = "Voles", ylab = "Voles", xlab = "Semester")
plot(res[[2]][ii, ], type = 'l', main = "Weasels", ylab = "Weasels", xlab = "Semester")
ii = ii + 1
```

The model contains the real dataset described in ([Fasiolo and Wood, 2015](http://arxiv.org/abs/1511.02644)),
which can be loaded using `data(voles_data)`. 
To fit the model with Synthetic Likelihood (Wood, 2010), we need to load the [synlik package](https://cran.r-project.org/web/packages/synlik/index.html) and
then create an object of class `synlik`. Notice that here we use the wrapper `volesWrap`, rather than `volesSimulator`, 
for compatibility with `synlik`.

```R
# This is needed because the process is observed only in 
# the months of June and September, for 45 years
obsTiming <- sort(c(6 + seq(0, 12 * 44, by = 12), 9 + seq(0, 12 * 44, by = 12)))

# Create "synlik" object
voles_sl <- new("synlik",
                simulator = volesWrap,
                summaries = volesStats,
                param = log(c(r = 4.5, e = 1, g = 0.2, h = 0.15, a = 8,
                              d = 0.06, s = 1, sigmaProc = 1.5, phi = 100)),
                extraArgs = list("nObs" = 12 * 45,  
                                 "nBurn" = 12 * 10, 
                                 "monthsToSave" = obsTiming)
)

# Put real data into object
data(voles_data)
voles_sl@data <- round(voles_data$popDensity * 10)
voles_sl@extraArgs$obsData <- round(voles_data$popDensity * 10)
```

We can than fit the model using a Metropolis-Hastings sampler.

```R
# Load MCMC proposal
data(voleFullProp_sl)

# Run MCMC on synthetic likelhood (only 100 iterations and it takes a while)
set.seed(51554)
tim <- proc.time()
voles_sl <- smcmc(voles_sl, 
                  initPar = voleFullTrueInit_sl,
                  nsim = 1000,
                  niter = 100, 
                  burn = 0,
                  priorFun = voleFullPrior_sl, 
                  propCov = 0.1 * voleFullProp_sl)
)
tim <- proc.time() - tim

plot(voles_sl)
```


References
----------------------------
  
  * Fasiolo, M and Wood, S. N. (2015). Approximate methods for dynamic ecological models, 
  ArXiv http://arxiv.org/abs/1511.02644. To appear in the Handbook of Approximate Bayesian Computation by S. Sisson, L. Fan, and M. Beaumont.

  * Turchin, P. and S. P. Ellner (2000). Living on the edge of chaos: population
dynamics of fennoscandian voles. Ecology 81 (11), 3099–3116.

  * Wood, S. N. (2010). Statistical inference for noisy nonlinear ecological dynamic systems. Nature 466 (7310), 1102–1104.
