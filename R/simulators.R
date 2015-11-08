
###############
#### Simulators for voles model (Turchin and Ellner (2000): LIVING ON THE EDGE OF CHAOS: POPULATION DYNAMICS OF FENNOSCANDIAN VOLES)
###############
##### Basic simulator 
#' Simulates from the voles models of Turchin and Ellner (2000).
#' 
#' @param nObs Length of each simulated time series in months.
#' @param nsim Number of simulations from the model.
#' @param param If model == "full" it must be a vector of length 9 or a matrix of size nsim by 9, if model == "standard" 
#'        it must be a vector of length 7 or a matrix of size nsim by 7.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters.
#' @param extraArgs Only needed for compatibility with the "synlik" package. 
#' @param model If model == "full" the function will simulate from the full model described in eq. 2 of Turchin and Ellner (2000),
#'              model == "standard" the function will simulate from the full model described in eq. 6 of Turchin and Ellner (2000).                   
#' @param nBurn Number of initial steps to be discarded before saving the following nObs steps.  
#' @param nSteps Number of finite differences steps between 2 observations (an interval which corresponds to a month).                
#' @param ... Only needed for compatibility with the "synlik" package. 
#' 
#' @return A list where [[""voles"]] and [["weasels"]] are two matrices of size nsim by nObs where each row is a trajectory
#'         simulated from the model.
#' @references Turchin and Ellner (2000), LIVING ON THE EDGE OF CHAOS: POPULATION DYNAMICS OF FENNOSCANDIAN VOLES.                    
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @examples
#' library(volesModel)
#'
#' ########## Simulating from full model 
#' set.seed(76735756)
#' nSteps <- 50
#' nObs <- 500
#' res <- volesSimulator(nObs = nObs, nsim = 100, nBurn = 120, model = "full", 
#'                       param = log(c(r = 4.5, e = 1, g = 0.12, h = 0.1, a = 8, d = 0.04, s = 1.25, 
#'                                     sigmaProc = 1.5, phi = 100)), 
#'                       nSteps = nSteps, T0 = 0)
#' 
#' ii = 1
#' par(mfrow = c(2, 1))
#' plot(res[[1]][ii, ], type = 'l', main = "Voles", ylab = "Voles", xlab = "Semester")
#' plot(res[[2]][ii, ], type = 'l', main = "Weasels", ylab = "Weasels", xlab = "Semester")
#' ii = ii + 1
#' 
#' ########## Simulating from simpler model
#' # Standard Model (maybe 99 steps is fine)
#' set.seed(452535626)
#' res <- volesSimulator(nObs = nObs, nsim = 100, nBurn = 120, model = "standard", 
#'                       param = log(c(r = 3.07, e = 1.58, beta = 1.54, d = 0.6, alpha = 7.6, sigmaProc = 1.5, phi = 100)), 
#'                       nSteps = nSteps, T0 = 0, parall = F)
#' 
#' ii = 1
#' par(mfrow = c(2, 1))
#' plot(res[[1]][ii, ], type = 'l', main = "Voles", ylab = "Voles", xlab = "Semester")
#' plot(res[[2]][ii, ], type = 'l', main = "Weasels", ylab = "Weasels", xlab = "Semester")
#' ii = ii + 1
#' 
#' ############# Using full model with synthetic likelihood
#' # Loading real data
#' data(voles_data)
#' 
#' 
#' 

#' 
#' # Loading real dataset
#' data(voles_data)
#' 
#' # Setting 
#' @export volesSimulator
#' 
volesSimulator <- function(param, nsim, nObs, extraArgs = NULL, model = "full", nBurn = 12 * 10, nSteps = 10, T0 = 0, ...)
{
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  if( !is.null(extraArgs) ){ 
    if( !is.null(extraArgs$model) ) model <- extraArgs$model
    if( !is.null(extraArgs$nBurn) ) nBurn <- extraArgs$nBurn
    if( !is.null(extraArgs$T0) ) T0 <- extraArgs$T0
  }
  
  if( nBurn %% 12 != 0 ) warning("nBurn % 12 != 0 but nBurn must be a multiple of 12, 
                                  otherwise the first observed month doesn't correspond to January and so on!")
  
  tmp <- switch(model,
                full = .Call("volesModel_volesFullCpp", nMon_ = nObs, nSimul_ = nsim, nBurn_ = nBurn, param_ = param, nSteps_ = nSteps, T0_ = T0,
                             randInit_ = TRUE, startVole_ = 0.0,  startWeas_ = 0.0, addObsNoise_ = TRUE, PACKAGE = "volesModel"),
                standard = .Call("volesModel_volesStdCpp", nMon_ = nObs, nSimul_ = nsim, nBurn_ = nBurn, param_ = param, nSteps_ = nSteps, T0_ = T0,
                                 randInit_ = TRUE, startVole_ = 0.0,  startWeas_ = 0.0, addObsNoise_ = TRUE, PACKAGE = "volesModel"),  
                stop("Wrong model name")
  )
  
  # Compile and load Rcpp code on the fly if needed
#   if( !exists("volesFullCpp") ) sourceCpp("voles_pomp_simulator.cpp")
#   tmp <- volesFullCpp(nMon = nObs, nSimul = nsim, nBurn = nBurn, param = param, nSteps = nSteps, T0 = T0,
#                       randInit = TRUE, startVole = 0.0,  startWeas = 0.0, addObsNoise = TRUE)
  
  return(tmp)
}


### Wrappers used by synlik 
volesWrap <- function(param, nsim, extraArgs, nObs, ...)
{
  if( !is.loaded("volesModel") ) library("volesModel")
  nObs <- extraArgs$nObs
  stopifnot( !is.null(nObs) )
  
  simul <- volesSimulator(param = param, nsim = nsim, extraArgs = extraArgs, nObs = nObs, ...)[["voles"]]
  
  # Save only the months with indexes in extraArgs$monthsToSave
  if( !is.null(extraArgs$monthsToSave) ){ 
    return( simul[ , extraArgs$monthsToSave, drop = FALSE] )
  } else {
    return( simul )
  }

}






