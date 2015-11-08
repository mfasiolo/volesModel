#################################
# Settings for the Voles model
#################################

###########
##### Settings for simulated data
###########

## Synlik
voleTrueParams <- log(c(r = 4.5, e = 1, g = 0.2, h = 0.15, a = 8, d = 0.06, s = 1, sigmaProc = 1.5, phi = 100))

voleFullPrior_sl <- function(input){ sum( input ) +
                               dnorm (exp(input[1]), 5, 1, log = TRUE)          +
                               dnorm (exp(input[2]), 1, 1, log = TRUE)          +
                               dexp  (exp(input[3]), 7, log = TRUE)             +
                               dgamma(exp(input[4]), 4, 40, log = TRUE)         +
                               dnorm (exp(input[5]), 15, 15, log = TRUE)        +
                               dnorm (exp(input[6]), 0.04, 0.04, log = TRUE)    +
                               dnorm (exp(input[7]), 1.25, 0.5, log = TRUE)     +
                               dunif (exp(input[8]), 0.5, 30, log = TRUE)       +
                               dunif (exp(input[9]), 20, 300, log = TRUE)       }

voleFullSimulInit_sl <- log(c(r = 3.5, e = 1.2, g = 0.1, h = 0.1, a = 10, d = 0.06, s = 1.1, sigmaProc = 2.5, phi = 80))

voleStdSimulInit_sl <- log( c(r = 4, e = 1.5, beta = 1.5, d = 0.6, alpha = 8, sigmaProc = 1.5, phi = 100) )


###########
##### Settings for real data
###########

voleFullTrueInit_sl <- log(c(r = 3.5, e = 1.2, g = 0.1, h = 0.1, a = 10, d = 0.06, s = 1.1, sigmaProc = 2.5, phi = 300))

voleFullTrueInit_pmcmc <- log(c(r = 3.5, e = 1.2, g = 0.1, h = 0.1, a = 10, d = 0.06, s = 1.1, sigmaProc = 2.5, phi = 300, 
                                voles.0 = 2, weasels.0 = 0.2))
