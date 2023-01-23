##########################################################################################
# Inverse calibration  ---------------------------------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
## Informal inversion
# - Some measure of model fit (objective function)
  # - e.g., Nash-Sutcliff etc.

  # GLUE: MCMC on informal target functions. Result: 'generalized uncertainty'

## Statistical inversion
# - Likelihood, where there is an assumption about distribution of the data
# - Bayesian as an add-on: prior uncertainty



# Recipe ------------------------------------------------------------------
# 1. Plot, check, omit.na data.
# 2. Check model output, do sensitivity analysis.
# 3. Decide on a likelihood (try to guess variance model [just like in linear regression!], covariance structure, based on 1. and 2.)
#    (Make the errors with a simple likelihood!)
# 4. MCMC sample, check burnin, check convergence, check residuals
# 5. Prediction and backward calculation of uncertainties (CI/PI etc. for time series)

## Hacks: don't calibrate all parameters,

# Stochastic models ----------------------------------------------------
#
# - Is the optimizer robust against stochastic likelihoods?
#
# - If the model is stochastic enough to produce the true variance:
#   - Pseudomarginal theorem! MCMC works with stochastic models!
#   - just use it without a statistical error/additional likelihood
#   - synthetic likelihood approach: take both mean and likelihood from the process


# Library -----------------------------------------------------------------
library(Rpreles)
library(lattice)

library(BayesianTools)
library(DEoptim)
library(numDeriv) # for derivative sensitivity analysis
# library(sensitivity)


# Load and inspect data ---------------------------------------------------------------
load('Dynamic/Data/boreal_site.rdata', verbose = T)
str(s1)
# VPD: Vapor pressure deficit
# GPPobs: from carbon fluxes

pairs(s1[c('GPPobs', 'ETobs', 'TAir', 'Precip', 'VPD', 'PAR')])
xyplot(GPPobs + ETobs ~ DOY, data = s1)

load('Dynamic/Data/par.rdata', verbose = T)
Refpar <- par


# Model ---------------------------------------------------------------
# Select parameters to include in the our analysis
selectparameter <- c(5:11, 14:18, 31, 32) # Indices for PRELES parameters (and also statistical parameters in our likelihood)
Refpar$name[selectparameter]

# function for selected parameter replacement
runPreles <- function(par){
  p <- Refpar$def
  p[selectparameter] <- par
  model <-  PRELES(DOY = s1$DOY,
                   PAR = s1$PAR,
                   TAir = s1$TAir,
                   VPD = s1$VPD,
                   Precip = s1$Precip,
                   CO2 = s1$CO2,
                   fAPAR = s1$fAPAR,
                   LOGFLAG = 0,
                   p = p[1:30],
                   control = 1)
  return(model)
}


modelvalues <- runPreles(Refpar$def[selectparameter])

plot(modelvalues$GPP, type = "l", col = "red")
points(s1$GPPobs, pch = 4)

# observed vs. predicted
barplot(modelvalues$GPP - s1$GPPobs) # sic!


plot(modelvalues$GPP, s1$GPPobs)
fit <- lm(s1$GPPobs ~ modelvalues$GPP)
summary(fit)
abline(fit, col = "red")
abline(0, 1)

# better centered
BayesianTools::GOF(observed = s1$GPPobs, predicted = modelvalues$GPP,plot = T) # only in latest version right way around


# Informal objective (GLUE) -----------------------------------------------
# Calibration on simple distance function (RSS): GLU (generalized likelihood uncertainty estimation)
# 'generalized' == informal

getCalTarget <- function(par){
  rss  <-  sum((runPreles(par)$GPP - s1$GPPobs)^2)
  return(rss)
}

## optimization with DEoptim (fast and highly parallelizable)
fit_deoptim <- DEoptim::DEoptim(getCalTarget,
                 lower = Refpar$min[selectparameter],
                 upper = Refpar$max[selectparameter])

# getCalTarget(fit_deoptim$...)

## Do the same with MCMC in order to


# BayesianTools::GOF()

# if used with a made up general objective function:
# interpretation of the uncertainty is problematic scale on the y-axis is completely arbitrary
# cite Keith Beven to defent this procedure: GLUE. (But actually there is no justification)
# Argument: 'there is no perfect model' -> 'there is no likelihood anyway'.

# When doing this, uncertainty can be arbitralily changed by dividing the likelihood by some value.

## GlUE:
## + Often not slower than optimizer
## + Shape of objective is visible (skew)
## - however, width of parameter distribution is arbitrary, not equivalent to uncertainty!






# Formal statistical calibration ------------------------------------------
## THERE IS A STATISTICAL MODEL ON TOP OF THE PROCESS MODEL (e.g. normal)

### Bayes theorem:
## p(parameter | data ) = p(data | parameter) * p(parameter)
##                           likelihood            prior
##
## Can only statistically work if the model produces probabilities.
## Often there has to be added a statistical model on top of the process based model.


### flat prior:
## Laplace: principle of insufficient reason.
## But sd should not be set flat!

## UK tradition:
## Advantage of Beta priors (scaling inside likelihood):
## posterior is prior again and can then be scaled

likelihood <- function(par){
  # rss  <-  sum((runPreles(par)$GPP - s1$GPPobs)^2 / 2 / par[13]^2)
    # actually this is the same as lognormal (compare to informal)
    # if the sd (the width) is fit with it, scaling uncertainty to the right level
  ll <- dnorm(s1$GPPobs, mean = runPreles(par)$GPP, sd = par[13], log = T)
  return(sum(ll))
}

# not fitting sd is not Bayesian but GLUE, uncertainties are not right

likelihood(Refpar$def[selectparameter]) # Always check the likelihood to take the relatively right values for known parameters


fitMLE <- DEoptim(likelihood, lower = Refpar$min[selectparameter], upper = Refpar$max[selectparameter])
fitMLE$optim$bestmem

library(BayesianTools)
setup <- createBayesianSetup(likelihood, lower = refPars$min[selectparameter], upper = refPars$max[selectparameter], names = refPars$names[selectparameter])

traceplot(MCMCout, includesProbabilities = T)

# runMCMC


# DEzs sampler: internal chains are not independent if the algorithm has not converged yet!!!
# check by multiple MCMCs to compare the convergence. The internal changed will be correlated between each other
# if time constraint, nrChains = 1 in two different sessions
# After convergence the internal chains will be independent and it is equally probable they are all in the same place, even across different MCMCs


## Checks RECIPE:
# 1. traceplots (FOR PROCESS 1: BURNIN? Traces should mix randomly)
# 2. Posterior over time (FOR PROCESS 1: BURNIN! No trend in posterior anymore )
#
# 3. Gelman diagnostics (FOR PROCESS 2: SAMPLED LONG ENOUGH after BURNIN?)
#        - Measurement whether the area of good fit is fully explored.
#        - For this purpose Neff, effective sample size!
#        - No. of iterations does not equal the no. of independent draws from posterior.
#
#        - How to get the sampler get sampling better:
#           - often this is due to parameter correlations:
              correlationPlot(MCMCout, start = 10000, thin = 20)

#           - Correlation reduction strategies
#             - 1. Increase data
#               - if there is little data and many parameters: the optimization problem is underdetermined. Correlated parameters
#             - 2. Tune sampler
#               - DEzs: tune the z-matrix

              startZ =  matrix(rnorm(13*1000, Refpar$def[selectparameter], sd = (refPars$max[selectparameter] - refPars$min[selectparameter])/100), ncol = 13, byrow = T) # just choose some sd arbitrarily divide by some values

              # run with less iterations
              MCMCout <- runMCMC(setup, settings = list(iterations = 100000, Z = startZ))
              traceplot(MCMCout, includesProbabilities = T)
# Convergence: roghly psf < 1.05


# now with stats: real parameter uncertainty
marginalPlot(MCMCout, start = 10000, prior = T, singlePanel = T, type = "v")

# for reporting
MAP(MCMCout)

## predictions
parametersample <- getSample(etc)
makePrediction <- function(par){
  gpp_predicted <- runPreles(par)$GPP
  sum(gpp_predicted/2) #?
}

gpp_priorprediction <- apply(setup$prior$sampler(nsamples), 1, makePrediction)
gpp_posteriorprediction <- apply(parametersample, 1, makePrediction)

hist(gpp_priorprediction)
hist(gpp_posteriorprediction, add = T)

## Find out which parameters are responsible for uncertainty
# fit_rf <- randomForest::randomForest(getSample(etc), prediction)
# varImpPlot(fit_rf)

###########
## Posterior predictive checks
# Check with residuals (since there are statistical assumptions!)
##########

simulateObservation <- function(par){
  modelobs <- runPreles(par)$GPP
  obs <- rnorm(modelobs, mean = modelobs, sd = par[13])
  return(obs)
}

simulateObservation(Refpar$def[selectparameter])
Data_simulated <- apply(parametersample, 1, simulateObservation)

## transform into quantile residuals - distribution should be uniform
# for informative priors the model residuals will not be exactly uniform!
res <- DHARMa::createDHARMa(Data_simulated,
                            observedResponse = s1$GPPobs,
                            fittedPredictedResponse = apply(Data_simulated, 1, median))
testZeroInflation(res)
plotResiduals(pred = s1$VPD, res$scaledResiduals)
testTemporalAutocorrelation(res, time = 1:length(s1$VPD)) # index
# acf()

# accounting for bad residuals in the likelihood (autocorrelation etc.)
# will likely result in a worse RSME (= simple normal likelihood!) because the objective function for optimization is now different!

# for something continuous positive (growth rates etc.) normally gamma distribution (or lognormal), also to improve residuals!
# (model still can not fit to the data in winter, where there are no negative observation)

### improve likelihood:
likelihoodVariance <- function(par){
  modelvalue <- runPreles(par)$GPP
  sdmodel <- par[13] + par[14] * modelvalue$GPP # without loglink (exp) when we are sure like here, that this can't be negative! (because of priors)
                                                # sd intercept + slope dependent on GPP
  ll <- dnorm(s1$GPPobs, mean = modelvalue, sd = sdmodel, log = T)
  return(sum(ll))
}

# Fit
# Do the residudal stuff again


# Predictive Intervals ----------------------------------------------------------------------
## C interval
## Prediction interval

out <-  BayesianTools::getPredictiveIntervals(sampleParameters,
                             getPredictionPreles,
                             numSamples = 1000,
                             quantiles = c(0.025, 0.975),
                             error = errorF)

BayesianTools::plotTimeSeries(observed = s1$GPPobs,
               predicted = runPreles(mapValue$parametersMAP)$GPP,
               confidenceBand = out$posteriorPredictiveCredibleInterval,
               predictionBand = out$posteriorPredictivePredictionInterval)



# Model selection ---------------------------------------------------------

## if unsure, just use the more complicated likelihood instead of fitting two MCMCs
## rather do it with priors on parameters (e.g. spike'n'slap)

# 1. if the problem can be expressed through a nested/big model, just do it.
# 2. The complicated odel selection possibilities:
#   1. Information theoretical MS: DIC, WAIC (to account for how much does the prior decreases model complexity, provided by BayesianTools), READ: GELMAN et al. 2014. Understanding predoctove onfprmation …
    DIC(MCMCout); DIC(MCMCout2)
    WIC(MCMCout); WIC(MCMCout2) # review-wise the preferred fancy measure # but must have likelihood with option sum = F
    # but: DIC stabilizes a lot slower than the median, psrf = 1.05 are for median values. More sampling for small probs!
#   2. Bayesian MS = Bayes Factor
    # For Bayes Factor: Marginal likelihood, integration/average over prior space, is needed
    ?marginalLikelihood()
    # but be aware: for uncalibrated models, wide priors, in depends wholly on prior
    # solution: fractional Bayes factor, somehow take half of the data to calibrate


