##########################################################################################
# Beta regression model and DHARMa -------------------------------------------------------
##########################################################################################

library(glmmTMB)
library(DHARMa)

# Reparameteraization of the beta distribution ------------------------------------------
## For parameterization Silvia Ferrari & Francisco Cribari-Neto (2004) Beta Regression for Modelling Rates and Proportions, Journal of Applied Statistics, 31:7, 799-815, DOI: 10.1080/0266476042000214501
## mu = E = p/(p+q)
## phi = p + q
## below, this was solved with phi != 0


dbeta2 <- function(x, mu, phi, log = F) {
  
  p <-  mu * phi
  q <-  phi - mu * phi
  
  dbeta(x, p, q, log = log)
}


rbeta2 <- function(n, mu, phi) {
  
  p <-  mu * phi
  q <-  phi - mu * phi
  
  rbeta(n, p, q)
}


# Sim ---------------------------------------------------------------------

n <- 1000

x <- runif(n, -1, 1)
f <- rep(c("apple", "banana"), length.out = n)
X <- model.matrix(~ x * f, data.frame(x = x ,f = f))
b <- c(0.5, 0.5, -0.5, 0.3)
b_phi <- c(2, 3, 10, -5) # ok fits
b_phi <- c(-5, 20, -10, -20) # bad fits for particular low precision per factor, numerical fit problems?


y_hat <- plogis(X %*% b)
phi_hat <- exp(X %*% b_phi)

y <- rbeta2(n, y_hat, phi_hat)
y <- scales::rescale(y, to = c(.Machine$double.eps, 1 - .Machine$double.eps))

plot(y ~ x, col = as.factor(f))


# Fit ---------------------------------------------------------------------
m <- glmmTMB(y ~ x * f, family = beta_family) # wrong model spec
m_disp <- glmmTMB(y ~ x * f, dispformula = ~ x * f, family = beta_family) # true model spec
m_disp_gaussian <- glmmTMB(y ~ x * f, dispformula = ~ x * f, family = gaussian)


# points(predict(m) ~ x, col = "blue")
# points(predict(m_disp) ~ x, col = "red")

r <- simulateResiduals(m)
r_disp <- simulateResiduals(m_disp) # true model spec
r_disp_gaussian <- simulateResiduals(m_disp_gaussian)


plot(r)
plot(r_disp) # true model spec
plot(r_disp_gaussian) # true model spec


boxplot(r_disp$scaledResiduals ~ r_disp$fittedModel$frame$f)


# Inspect simulations
plot(getSimulations(m_disp) ~ y)
plot(getSimulations(m_disp) ~ x, col = as.factor(f)) # == points(simulate(m_disp)[,1] ~ x, col = "red")

points(y ~ x, col = "blue")
