# Set up ------------------------------------------------------------------

## Construct a matrix of transition rates by rows
stages <-  c("rosette", "shoot", "seed")
M <- rbind(c(0.95,      0.21,    0.01),
           c(0.33,      0,      0   ),
           c(0,        10,      0.25 ))
colnames(M) <- rownames(M) <- stages
M

## Construct a vector, representing an inital population
population_init <- c(rosette = 1, shoot = 10, seed = 50)
population_init



## Matrix multiplication ------------------------------------------------------

## Matrix multiplication of the inital population vector yields the population vector after 1 year
population_2 <-  M %*% population_init
population_2
population_2 <- c(population_2) ## c() "concatenates" the 1-column-matrix into a new vector
population_2

## ... after 2 years
population_3 <- c(M %*% population_2)
population_3



# Forward simulation of a timeseries ------------------------------------------------------
n_years <- 10
n_stages <- length(population_init)

Timeseries <- matrix(NA, nrow = n_years, ncol = n_stages) ## construct an empty matrix
Timeseries[1,] <- population_init ## fill the first row of the matrix with the initial population
times <- 1:(n_years-1) ## make a vector of integer times from 1 to n_years-1

for (t in times) {
  population_before <- Timeseries[t,] ## this is automatically extracted as a vector
  Timeseries[t+1,] <- M %*% population_before
}

Timeseries

matplot(Timeseries, type = "b", pch = c("r", "F", "s"))
  ## execute '?matplot' to read what matplot() does

## The population development is dependent on the population structure.
## assuming a population where everything aboveground has been killed and only a few seeds exist
## Can the population recover to grow again exponentially after 3 years?

population_init <- c(rosette = 0, shoot = 0, seed = 5)
Timeseries[1,] <- population_init ## fill the first row of the matrix with the initial population

for (t in times) {
  population_before <- Timeseries[t,]
  Timeseries[t+1,] <- M %*% population_before
}

Timeseries
matplot(Timeseries, type = "b", pch = c("r", "F", "s"))



# Matrix analytics  -----------------------------------------------

## Compute n eigenvectors and eigenvalues of M
M
e <- eigen(M)
e
  ## The values are ordered, the greatest eigenvektor ("dominant") comes first
  ## The i indicates the imaginary part of the complex numbers that come up in the solution of polynomials, which we can ignore for our population demographics.


## Check: Do the computed values satisfy the definition Mv = kv?
v <- e$vectors[,1]
k <- e$values[1]
v
k
Mv <- c(M %*% v)
kv <- k * v
all.equal(Mv, kv)
Mv
kv


## What is the overall growth rate of the population?
Re(k) ## Re(), the real part of the dominant eigenvalue


## What is the stable proportion of the population for exponential growth?
v <- Re(v) ## Re(), the real part of the dominant eigenvector has the proportions of any population that grows exponentially
v
stableprop <- v/sum(v) ## get the proportions relative to the unit 1
stableprop

## Test whether simulation yields leads to the same proportion of states
n_years <- 100
population_init <- c(rosette = 5, shoot = 3, seed = 1)
Timeseries <- matrix(NA, nrow = n_years, ncol = n_stages) ## construct an empty matrix
Timeseries[1,] <- population_init ## fill the first row of the matrix with the initial population
times <- 1:(n_years-1) ## make a vector of integer times from 1 to n_years-1

for (t in times) {
  population_before <- Timeseries[t,]
  Timeseries[t+1,] <- M %*% population_before
}

Timeseries
matplot(Timeseries[1:15,], type = "b", pch = c("r", "F", "s"))
## express all rows proportional to their sum:
Timeseries_prop <- t(apply(Timeseries, MARGIN = 1, function(row) row/sum(row)))
Timeseries_prop
## are all population states proportioned like the dominant eigenvector?
all.equal(Timeseries_prop[100,], stableprop) ## is the 100th row of the timeseries approximately equal to the stable proportion?

