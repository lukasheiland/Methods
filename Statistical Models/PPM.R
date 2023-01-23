
# Library -----------------------------------------------------------------
# library(BayesianTools)
library(spatstat)

# Simulation --------------------------------------------------------------

## Poisson point process:
## has 3 properties:

## 1) window
window <- list(
  type = "rectangle",
  xrange = c(-1, 1),
  yrange = c(-1, 1),
  units = list(singular = "m", plural = "m")
)
attr(window, "class") <- "owin"

## 1a) n, the number of points within the window.
## But in the simulation in spatstat::rpoispp(), n seems to be an emerging parameter from lambda, number of trials and rejection sampling.


## 2) point generating intensity function lambda(x,y); the most basic case lambda(x,y) = exp(a + bx)
## Rationale: lambda = n/Area <=> n = lambda*Area
## Thus Poisson process: N(Area) ~ Poisson(lambda*Area)

a <- 2
b_x <- 3
b_y <- -1.2

# lambda <- function(x, y) exp(a + b_x*x)
## And here is the point intensity function
lambda <- function(x, y) exp(a + b_x*x + b_y*y)
curve(lambda(x, y = 0), -1, 1)

## Here is a weird function to simulate other covariates
weirdfunction <- function(x,y){ 10 * x^2 + 5 * sin(10 * y) }

## And here is the point intensity function
lambda2 <- function(x, y) exp(a + b_x*x + b_y*y + weirdfunction(x, y))
curve(lambda(x, y = 0), -1, 1)


## Make a grid
# imageres <- 5
# Grid <- expand.grid(x = seq(window$xrange[1], window$xrange[2], length.out = imageres*window$xrange[2]),
#                     y = seq(window$yrange[1], window$yrange[2], length.out = imageres*window$yrange[2]))
# L <- matrix(mapply(lambda, Grid$x, Grid$y), nrow = imageres*window$xrange[2])
# image(L)


## Simulation
simPoints <- function(lambda, w = window, trials = 10000){
  x_win <- runif(trials, w$xrange[1], w$xrange[2])
  y_win <- runif(trials, w$yrange[1], w$yrange[2])

  imag <- as.im(lambda, w)
  summ <- summary(imag)
  lmax <- summ$max + 0.05 * diff(summ$range) # wtf?

  prob <- lambda(x_win, y_win)/lmax
  u <- runif(trials)
  retain <- (u <= prob)

  Points <- cbind(x = x_win[retain], y = y_win[retain])
  return(Points)
}

#

P <- simPoints(lambda)
plot(P, xlim = window$xrange, ylim = window$yrange)

## Fake a ppp object.
pointdata <- list(window = window, n = nrow(P), x = P[,"x"], y = P[,"y"])
attr(pointdata, "class") <- "ppp"
plot(pointdata)


## The real thing.
pointdata2 <- rpoispp(lambda, win = window)
plot(pointdata2)

# Fit ---------------------------------------------------------------------
?spatstat::ppm()
str(pointdata)


# fit the stationary Poisson process
# to point pattern 'nztrees'
fit <- ppm(pointdata ~ x + y)
summary(fit)

## Include a third covariate, will also work with an image
fit2 <- ppm(pointdata ~ x + y + weirdfunction)
