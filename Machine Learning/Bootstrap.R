##########################################################################################
# Bootstrap ---------------------------------------------------------------
##########################################################################################

## For calculation of CONFIDENCE INTERVALS

# Resample fractions of the data with replacement and fit models.
# The spread of the estimatesm will correspond to the estimate distribution/uncertainty.
# There is a proof, that this holds, regardless of any distributional theory.

# Use this for CIs of longer pipelines with multiple methods (e.g., including AIC selection), to get correct CI (and p-values)!

# For parametric boot strap, see there

# The boot strap ---------------------------------------------------------------
library(boot)
d <- rexp(100)
mean(d)
b <- boot(d,
          function(data, indices) mean(data[indices]), # fun
          R = 999) # Replicates
plot(b)
myconfidenceinterval <- sd(b$t)
