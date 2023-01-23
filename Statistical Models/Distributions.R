##########################################################################################
# Alternative parameterizations of Distributions ---------------------------------------
##########################################################################################

# Distributions, re-parameterized with mean and some dispersion parameter ------------------------------------------------------------



# Gamma -------------------------------------------------------------------
##  If x is gamma distributed, it is the sum of many exponentially-distributed variates. For example, the waiting time for many events of a Poisson process.
## The Gamma has more of a tail on the left, and less of a tail on the right;
## the far right tail of the lognormal is heavier and its left tail lighter!
##
## For relations of parameters alpha and beta (the standards in R and stan; sic!) to empirical mean E and Var
## shape == alpha; rate == inverse scale == beta
## $\alpha = \frac{E^2[X]} {\mathrm{Var}(x)}$
## $\beta = \frac{E[X]} {\mathrm{Var}(x)}$

dgamma2 <- function(x, mean, disp_inv) {
  dgamma(x, shape = disp_inv, rate = disp_inv/mean)
}

rgamma2 <- function(n, mean, disp_inv) {
  rgamma(n, shape = disp_inv, rate = disp_inv/mean)
}

rgamma2 <- function(n, mean, disp) {
  variance <- disp * (mean^2)
  shape <- (mean^2) / variance
  rate <- mean / variance
  rgamma(n, shape = shape, rate = rate)
}

dgamma2 <- function(x, mean, disp) {
  variance <- disp * (mean^2)
  shape <- (mean^2) / variance
  rate <- mean / variance
  dgamma(x, shape = shape, rate = rate)
}

dgamma3 <- function(x, mode, sd) {
  rate <-  ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
  shape <-  1 + mode * rate
  dgamma(x, shape = shape, rate = rate)
}


# NegBinomial -------------------------------------------------------------
## The precision parameter phi acts as an inverse dispersion parameter:
## Var[n] = mean + (mean^2/phi)
## phi = mean^2/(var - mean)
## Alternative parameterization in stan: https://mc-stan.org/docs/2_29/functions-reference/nbalt.html
rnbinom2 <- function(n, mean, precision) {
  rnbinom(n, mu = mean, size = precision)
}


# LogNormal ---------------------------------------------------------------
## if y ~ Normal, then exp(y) ~ LogNormal
## if x ~ LogNormal(), then y = log(x) ~ Normal()
rlnorm2 <- function(n, median, disp) {
  rlnorm(n, log(median), disp)
}

dlnorm2 <- function(x, median, disp) {
  dlnorm(x, log(median), disp)
}

rlnorm3 <- function(n, mean, disp) {
  rlnorm(n, log(mean) - 0.5 * disp^2, disp)
}

## cf.
# mean(rlnorm3(1000000, 0.5, 3))
# mean(rlnorm2(1000000, median = 0.5, 3))


# Beta ---------------------------------------------------------------

dbeta2 <- function(x, mean, samplesize) {
    dbeta(x, mean*samplesize, (1-mean)*samplesize)
}