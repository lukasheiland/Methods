##########################################################################################
# INLA ---------------------------------------------------------------
##########################################################################################
# GLM
# Bayesian, analytical estimation of posteriors.
# Only works with normal priors.



# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(INLA)
library(DHARMa)
library(glmmTMB)

D <- createData(400, numGroups = 20, family = gaussian(), temporalAutocorrelation = 100)
D$time_fac <- as.factor(D$time) # has to be a factor in glmmTMB

fit_inla <- inla(observedResponse ~ Environment1 + f(group, model = 'iid') + f(time, model = 'ar1'), family = "gaussian", data = D, control.predictor = list(link=1))
summary(fit_inla)


fit_glmmtmb <- glmmTMB(observedResponse ~ Environment1 + (1 | group) + ar1(time_fac + 0 | ID), data = D) # in lme4 this would be faster due to an analytical trick
