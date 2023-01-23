## The offset is added to the linear predictor so that ...

### This generative process
n <- 1000
lambda <- 10
x <- runif(n, 0, 2)
beta <- 3
os <- 1/3
os_rep <- rep(1/3, n)
y <- rpois(n, exp(lambda + os_rep))
y2 <- rpois(n, exp(log(lambda) + log(os_rep))) # log(a * b) = log(a) + log(b)
y3 <- rpois(n, exp(lambda + x * beta + log(os_rep)))

y4 <- rpois(n, exp(lambda + beta*x + log(1/(1+x)))) ## == log(beta*x * 1/(1+x))


### ... is recoverd by this model.
f <- glm(y ~ 1 + offset(os_rep), family = poisson(link = "log"))
f
f2 <- glm(y2 ~ 1 + offset(log(os_rep)), family = poisson(link = "log"))
f2
f3 <- glm(y3 ~ 1 + x + offset(log(os_rep)), family = poisson(link = "log"))
f3

f4 <- glm(y4 ~ 1 + x + offset(log(1/(1+x))), family = poisson(link = "log"))
f4
