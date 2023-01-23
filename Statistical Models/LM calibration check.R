
library(magrittr)

# Simulation --------------------------------------------------------------


simulateEstimation <- function(n = 200, b = 0.5, distribution = rnorm){

  x = runif(n)
  y = b * x + distribution(n)
  fit <- lm(y ~ x)

  s <- summary(fit)
  standarderror <-  s$coefficients[2,2] # confidence interval
  estimate <- s$coefficients[2,1]
  ci <- 1.96 * standarderror
  estimateinside <- abs(b - estimate) <  ci
  s <- c(b = estimate, p = s$coefficients[2,4], coverage = estimateinside)

  return(s)
}



# Check for bias --------------------------------------------------------------
b_true <- 0.5
Results <- replicate(500, simulateEstimation(b = b_true)) %>%
  t() %>%
  as.data.frame()

## Check whether the linear regression is an unbiased estimator for the slope, i.e. spreading around the true value
hist(Results$b, breaks = 50)
abline(v = b_true, col = "red")

## Explicitly calculating the bias
mean(Results$b) - b_true
sd(Results$b)


# Check for calibration of p-value  ---------------------------------------------
## Checking that the type I error rate (i.e. the rate at which p-values for the regression slope will become significant if the true slope is 0) is indeed 0.05.

b_true <- 0 # !!!
Results <- replicate(500, simulateEstimation(b = b_true)) %>%
  t() %>%
  as.data.frame()

#### p-value?
issignificant <- Results$p < 0.05
table(issignificant) # How often is the slope significant
mean(issignificant)
## logical vectors are interpreted as vectors of 0s and 1s
## -> percentage of significant models
## -> type I error rate


# Check coverage ----------------------------------------------------------
## The confidence interval is defined as the interval in which $level (95%) of the model estimates in repeated experiments land.
## Checking that nominal coverage of a CI is indeed 95%,

b_true <- 0.5 # !!!
Results <- replicate(1000, simulateEstimation(b = b_true)) %>%
  t() %>%
  as.data.frame()

iswithin <- Results$coverage
mean(iswithin) # Check!

############################################################################
# Destroying model assumptions ---------------------------------------------
############################################################################

# Simulating data with the wrong distribution -----------------------------------------------------------

### double exponential
rdexp <- function(n) {
  e <- rexp(n, rate = 0.01) * sample(c(-1,1), n, replace = T)  # randomly make a random number from an exponential distribution negative or positive
  # print(e)
  return(e)
}

b_true <- 0 #!!! for p-value check
Results <- replicate(1000, simulateEstimation(b = b_true,
                                             distribution = rdexp) #!!!
                     ) %>%
  t() %>%
  as.data.frame()

## 0. bias?
mean(Results$b) - b_true # bias?

## 1. p value?
mean(Results$p < 0.05) # strangely calibrated!

## 2. coverage
mean(Results$coverage) # strangely ok!




# AIC selection -----------------------------------------------------------
##  AIC model selection inflates type I error
## Yet to clean up/integrate!


getEstimate <- function(n = 100){
  x = matrix(runif(n * 30), ncol = 30)
  y = rnorm(100)
  dat = data.frame(x)
  dat$y = y
  fullModel <- lm(y ~ ., data = dat)
  redModel <- MASS::stepAIC(fullModel, trace = 0)
  x = summary(redModel)
  return(x$coefficients[-1,4])
}
out <- replicate(100, getEstimate())
mean(unlist(out) < 0.05)



