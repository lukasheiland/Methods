##########################################################################################
# Bootstrap ---------------------------------------------------------------
##########################################################################################

## For calculation of CONFIDENCE INTERVALS

# Resample fractions of the data __with replacement__ and fit models.
# This is based on the idea that the process of resampling from the data set on hand represents the sampling process that got the data set from the original population in the first place.
# How many bootstrap samples to get a good Monte Carlo approximation to the bootstrap result? Increase the size until you get convergence. No one number fits all problems.
# The spread of the estimates will correspond to the estimate distribution/uncertainty.
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
