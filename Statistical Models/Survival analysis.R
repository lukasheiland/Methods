##########################################################################################
# Survival analysis  ---------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
## Survival data is often right-censored and left-truncated.
# - Censoring: Information that a sample exceeds a limit is known, but not the value.
# - Truncation: Information that exceeds a limit is dropped.

# Library -----------------------------------------------------------------
library(survival)

# GLM  -------------------------------------------------------------------
## Response: Survival 1/0 at time.
## Question: is this the exact link function-offset combination?
glmer(survival ~ predictors + (1|individual), # RE only as dispersion correction
      offset = log(time), family = binomial(link = cloglog))

## Btw: library(lmec) for censored data.

# Survival --------------------------------------------------------------
## Response: Survival time. Usually with censoring.
data(veteran)
s <- Surv(veteran$time, # survival time
          veteran$status # censoring status
          )
plot(s)
sf <- survfit(s ~ 1, data = veteran) # Without predictor.
summary(sf)

## Kaplan-Meyer, estimate für survival time. Only a non-parametric curve, divided by groups.
sf_km <- survfit(s ~ celltype, data = veteran)
summary(sf_km)
plot(sf_km)

## Cox regression.
## Semi-parametric: Fit non-parametric hazard curves, compare the effect of predictors on them parametrically.
sf_cox <- coxph(s ~ celltype, data = veteran)
summary(sf_cox)
plot(sf_km)

## random effects: library(coxme)

