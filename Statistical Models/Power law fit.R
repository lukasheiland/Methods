##########################################################################################
# Power law binning  ---------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
#
#


# Library -----------------------------------------------------------------
library(poweRlaw)


# Simulate distribution ---------------------------------------------------
n <- 1000
support <- rplcon(n, 1, 1.1) # continuous observation, i.e. support.
support <- support[support < 1000] # just drop very high values to be able to plot
hist(support, breaks = 100)

support_log <- log(support)
# count <- table(cut(support_log, breaks = 20)) # the response 'y'.
hist <- hist(support_log, breaks = 10, plot = F)
count <- hist$counts
density <- hist$density
mid <- hist$mids

# Cut off the first and the last bin!
# Especially the last point because the empirical maximum ist not the actual maximum.
i_maxbin <- length(hist$mids)
i_minbin <- 1
D <- data.frame(x = mid[-c(i_minbin, i_maxbin)], y = density[-c(i_minbin, i_maxbin)])

fit <- lm(y ~ x, data = D[is.finite(D$y),])
summary(fit)
res <- fit$coefficients[2] # this is the slope


## The above could be wrapped inside a loop with a NA vector res[i] <- … and then hist(res), mean(res)


# Regression --------------------------------------------------------------
# The problem with an lm is, that it assumes homoscededastic error! (… while we are in log-log)
# This is thought to introduce a bias.


# MLE of the distribution ----------------------------------------------------
# This is the right way.
# When there are only bins, make a two-level process with binning in the model.

