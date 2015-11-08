
####################
######## VOLES: statistics for the vole model (Turchin)
####################

volesStats  <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(orderDist(tx, obsData, np=3, diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x),lag=c(6, 6, 6, 1, 1), power=c(1, 2, 3, 1, 2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x)) # mean values of Y, # of 0's
  X0 <- cbind(X0, t(slAcf(tx, max.lag = 5)))  ## autocovariances up to lag 5 (the first element is the variance)
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x,1,diff) ), 2, diff ) )), 1, sum)/2) #Number of turning points 
  #X0 <- t( log(apply(x, 1, function(input) spec.pgram(input, plot = FALSE)$spec)) )[ , 1:20, drop = FALSE] # Spectrogram
  
  return(X0)
}  ## end of ss.model.matrix function



