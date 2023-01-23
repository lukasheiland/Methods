library(gtools) # for logit()





# From linear data to I/0
n <- 100
x <- runif(n, -1, 1) # "Grad Celsius", always rescale when multiple predictors

a <- 4.2
b <- 1.4

y_hat <- a*x + b # y: some metabolism function dependent on deg Celsius
plot(y_hat ~ x)


y_logis <- plogis(y_hat) # probability of logistic dist (dist function), given quantiles y
plot(y_logis ~ x)

y <- rbinom(n, y_logis)

# imagine a process which will lead to death at
plot(y.strich ~ x)

plot(y ~ x)
points(log(y.strich/(1-y.strich)) ~ x, pch = 3, col = 2) # do the logit
# some high values somehow arent transformed back, maybe because they have just become 1.
# consider

exp(1:100)/(1+exp(1:100))

# add:
.Machine$double.eps

# glm(famnily = binomial) does even work with all prob data on [0, 1]
# but lets "polarize" anyway
# y.strich <- y.strich >= 0.5 # glm will probably do the same internally
# plot(y.strich ~ x)

# Fit a model
fit <- glm(y.strich ~ x, family = binomial)
plot(y.strich ~ x)
curve(predict(fit, newdata = data.frame(x), type = "response"), add = T)
summary(fit)

# !!!
# fit is the linear function in
# y = inverse.link(fit)

# because what is done is fit(link(y))

# intercept is 2.5 in logit space, what is survival prob at 0 deg celsius
inv.logit(2.5805)
# and at x = 1?
inv.logit(2.5805 + 2.2566)

plot(y ~ x)
points(logit(predict(fit, newdata = data.frame(x), type = "response")) ~ x, col = 3, add = T)

# ERROR TERM ADDS REGRESSION DILUTION!!!


# offenbar gibt es viele Lösungen für y = exp(y')/1-+exp(y')
# wenn ursprünglicher linearer TERM MIT INTERCEPT
# ??

exp(fit$coefficients) # the multiplicative change in the odds ratio for y=1 if the covariate associated with β increases by 1. 
a2 <- fit$coefficients[1]
b2 <- fit$coefficients[2]

fit2 <- glm(plogis(a2*x + b2) ~ x, family = binomial)

# y1 und y2 sind proportional zu x.
# y1 = ax + b AND y2 = lx + m
# x = (y1 - b)/a AND x = (y2 - m)/l
# (y1 - b)/a = (y2 - m)/l
# y1 = (a*(y2 - m)/l)-b
# y1 = (1/l)*a*y2 - (1/l)*a*m -b
# y1 und y2 sind auch untereinander proportional.




# Note: ll in anderem Skript
ll <- function(pars, model){
  -sum(dbinom(y1, 1, model(pars), log = T))
}



## NOTES zu LME
optimout <- optim(par = c(0.01),
                  fn = getLL,
                  method = "L-BFGS-B",
                  control = list(fnscale = -1))
# Hessian matrix, für CI, geht schief, wenn likelihood nicht 
# Quadratische Approximation ist nicht mehr verlässlich -> Bayes kompletter shape

library("bbmle")
getNLL <- function(rate) -getLL(rate)
out <- mle2(getNLL, start = c(rate = 0.01)) # expects minuslogl
summary(out) # Art Regressionstabelle, SD aus Hesse-Matrix
plot(profile(out)) # wie fällt likelihood ab, Annahme Normalverteilung







# now with quasibinomial
y <- a*x + b + x*e
plot(y ~ x)
y.strich <- exp(y)/(1+exp(y)) # inverse logit == logistic
y.strich <- y.strich > 0.5 # inverse logit == logistic
plot(y.strich ~ x)
summary(glm(y.strich ~ x, family = quasibinomial))

     