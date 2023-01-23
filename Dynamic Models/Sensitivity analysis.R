##########################################################################################
# Sensitivity analysis  ---------------------------------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
## Global SA
# vary all parameters -> hypervolume
# - Morris screening: random initial hypervolume points, random jumps

## Local SA
# vary only one parameter

## Semilocal SA
# one factor at a time: one-dimensional cuts through parameter space


## PRELES ##
# A simple semi-empirical ecosystem carbon and water balance model.
# The model predicts gross primary production and evapotranpiration (and soil water balance) based on embedded empirical relationships and basic meteorological inputs.

# Library -----------------------------------------------------------------
library(Rpreles)
library(lattice)

library(BayesianTools)
library(numDeriv) # for derivative sensitivity analysis
library(sensitivity)


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
selectparameters <- c(5:11, 14:18, 31) # Indices for PRELES parameters
Refpar$name[selectparameters]

# function for selected parameter replacement
runPreles <- function(par){
  p <- Refpar$def
  p[selectparameters] <- par
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


modelvalues <- runPreles(refpars$def[selectparameters])

plot(modelvalues$GPP, type = "l", col = "red")
points(s1$GPPobs, pch = 4)

# observed vs. predicted
plot(modelvalues$GPP, s1$GPPobs)
fit <- lm(s1$GPPobs ~ modelvalues$GPP)
barplot(modelvalues$GPP - s1$GPPobs) # sic!
summary(fit)
abline(fit, col = "red")
abline(0, 1)

# better centered
BayesianTools::GOF(observed = s1$GPPobs, predicted = modelvalues$GPP,plot = T)


# Sensitivity analysis -----------------------------------------------------------

### 1. Total local SA: calculate derivatives around reference point
getSensTarget <- function(par){
  summedgpp  <-  sum(runPreles(par)$GPP)
  return(summedgpp)
}

# in pricniple, we could just take the numerical derivatives
derivative <- numDeriv::hessian(getSensTarget, Refpar$def[selectparameters])
fields::image.plot(log(abs(derivative + 1)))

### 2. One factor at a time
setup <- BayesianTools::createBayesianSetup(getSensTarget,
                                            lower = Refpar$min[selectparameters],
                                            upper = Refpar$max[selectparameters],
                                            names = Refpar$name[selectparameters])

# plotSensitivity(setup, equalScale = F) # works from BayesianTools version > .1.6

### 3. Proper global SA
# library(sensitivity)
getSensTarget_par <- BayesianTools::generateParallelExecuter(getSensTarget)
# convenience function for parallelization expects function, like an apply
morrisresult <- sensitivity::morris(getSensTarget_par$parallelFun, # a faster 'screening' method
                                    Refpar$name[selectparameters],
                                    r = 1e3, # replicates
                                    design = list(type = "oat", levels = 5, grid.jump = 3),
                                    binf = Refpar$min[selectparameters],
                                    bsup = Refpar$max[selectparameters],
                                    scale = T) # IMPORTANT parameter! Describes whether to scale the inputs to the range [0,1].
# With scale = F it becomes independent of input uncertainty ranges
# (wide input uncertainty/high sensitivity both can lead to high sigma and mu)

# alpha varies stronger than beta in parameter space

# high number of replicates needed to become stable, dependent on the model being sensitive in certain areas
plot(morrisresult)

### lm method on outputs (maybe including quadratic effects)
# similar to Hessian matrix because including both first and second derivative
# direction of effect is easily seen from regression table

parameters <- setup$prior$sampler(1000) # will sample runif(min, max) from the parameter space; actually prior
colnames(parameters) <- Refpar$name[selectparameters]
GPP_result <- getSensTarget_par$parallelFun(parameters) # automatically created by BS
D <- data.frame(cbind(scale(parameters, scale = T), GPP = GPP_result))

summary(fit <- lm(GPP ~ (.)^2, data = D))

### some sobol methods
x1 <- data.frame(setup$prior$sampler(100)) # get two indendent samples from the same 'prior'
x2 <- data.frame(setup$prior$sampler(100))

sobolanalyses <- sobol(getSensTarget_par$parallelFun, X1 = x1, X2 = x2, order = 2) # nothing different but a regression table
sobolanalyses <- sobolSalt(getSensTarget_par$parallelFun, X1 = x1, X2 = x2)
plot(sobolanalyses, choice = 1)
ggplot2::ggplot(sobolanalyses, choice = 1)


