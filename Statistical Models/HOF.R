# Simulate  -----------------------------------------------------------
logist <- function(x) 1/(exp(x) + 1)

n <- 1000

a <- 2
b <- -4
c <- 0.7
d <- 1.6

x <- runif(n, -1, 1)
y1 <- rbinom(n, 1, logist(a + b*x + c*x^2))


# . ---------------------------------------------------------------------
## HOF algo
#
# 1. MLE fits of the formulae
# 2. Model selection with AICc
#    Within levels of equal number of parameters (III and IV, V and VI), model selection is only dependent on data fit (deviance).
# 3. Bootstrap: Model types will be collected, from bootstrapped data sets for which the original occurrences are re-sampled with replacement until the ori- ginal number of occurrences is reached.
#    The most frequent model type of these random runs is chosen, if it is different from the originally estimated model type.

m1 <-  function(pars){
  logist(pars)
}
m2 <- function(pars){
  logist(pars[1] + pars[2]*x)
}
m3 <- function(pars){
  1/((1+exp(pars[1] + pars[2]*x)) * (1+exp(pars[3])))
}
m4 <- function(pars){
  1/((1+exp(pars[1] + pars[2]*x)) * (1+exp(pars[3] + pars[2]*x)))
}
m5 <- function(pars){
  1/((1+exp(pars[1] + pars[2]*x)) * (1+exp(pars[3] + pars[4]*x)))
}

modellist <- list(m1, m2, m3, m4, m5)
npars <- c(1, 2, 3, 3, 4)

ll <- function(pars, model){
  -sum(dbinom(y1, 1, model(pars), log = T))
}

results <- mapply(function(m, n) {
  opt <-  optim(rep(1, n), ll, model = m)
  ml <- ll(opt$par, m)
  return(list(opt, ml))
}, modellist, npars, SIMPLIFY = F)

sapply(1:length(modellist), function(x) results[[x]][[2]])


# Random forest -----------------------------------------------------------
library(ranger)

D <-  data.frame(y = y1, x = x)
i <- sample.int(n, n/2)

rf <- ranger(as.factor(y) ~ x, data = D[i,], probability = T)
pred <- predict(rf, D[-i,])

llrf <- - sum(dbinom(D$y[-i], 1, as.numeric(pred$predictions[,2] + .Machine$double.eps), log = T))
llrf


# eHOF ---------------------------------------------------------------------
library('eHOF')
hof <- HOF(y1, x)
plot(hof)
eHOF::Para(hof)
