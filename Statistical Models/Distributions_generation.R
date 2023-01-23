
## Probability distributions are related and can be derived from each other.
## See, e.g., here https://www.johndcook.com/blog/distribution_chart/,
## and here http://www.math.wm.edu/~leemis/2008amstat.pdf.


##########################################################################
# Discrete distributions ------------------------------------------------
##########################################################################


## Bernouilli --------------------------------------------------------------

bernoulli <- function(n, p) {
  rbinom(n, size = 1, p)
}

hist(bernoulli(100, 0.9))



## Binomial ----------------------------------------------------------------
binomial <- function(n, size, p) {
    replicate(n, sum(bernoulli(size, p)))
}

hist(binomial(100, 3, 0.3))
hist(rbinom(1000, 100000000000, 0.3))


## Poisson ----------------------------------------------------------------
poisson <- function(n, lambda) {
  ## Poisson experiment is the number of successful bernoullis in a sample, where the average number of successes is lambda
  ## lambda == n*p
  size_binom <- 1000
  p <- lambda/size_binom
  replicate(n, rbinom(100, size_binom, p))
}

hist(poisson(1000, 1.5))
hist(rpois(1000, 1.5))


##  Negative binomial -------------------------------------------------------
## Generative approach
## https://stats.stackexchange.com/questions/176034/negative-binomial-distribution-vs-binomial-distribution

nbinomial_one <-function(size, prob) {
  trial <- 0
  successes <- 0
  while (successes < size) {
      trial <- trial+1
      successes <- successes + bernoulli(1, prob)
    }
  return(trial - successes) # return no of failures!
}

nbinomial <- function(n, size, prob) {
  replicate(n, nbinomial_one(size, prob))
}

hist(nbinomial(1000, 12, 0.2))
hist(rnbinom(1000, 12, 0.2))



##########################################################################
# Continuous distributions ------------------------------------------------
##########################################################################

#### Plot density
#### … of a rng f
d <- function(f, ..., n = 10000) {
  plot(density(f(n = n, ...)))
}


## Uniform -----------------------------------------------------------------
uniform <- runif

d(uniform, -1, 1)

## normal ------------------------------------------------------------------

## There are also approximations from discrete distributions like the binomial
## Here is an implementation just making use of the CLT and scaling like a dirty animal


normal <- function(n, mean, sd) {
  somenormal <- replicate(n, sum(rpois(1000, 10)))
  scalednormal <- scale(somenormal) * sd + mean
  return(c(scalednormal))
}
d(normal, -5, 2.3)

normal <- function(n = 100,
                   m = 20) {
  rowSums(replicate(m, runif(n, -1, 1)))
}

hist(normal())


## exponential -------------------------------------------------------------
## If X is a uniform random variable, -λ log X is an exponential random variable with mean λ.
## More generally, applying the inverse CDF of any random variable X to a uniform random variable creates a variable with the same distribution as X.
exponential <- function(n, mean) {
  -mean*log(runif(n, 0, 1))
}

d(exponential, 20)

