##########################################################################################
# MCMC sampler ---------------------------------------------------------------
##########################################################################################


# Library -----------------------------------------------------------------



# Sim linear model -------------------------------------------------------------------------
n <- 100

a <- 3
b <- 1
sigma <- 2


x <- (-(n - 1)/10):((n - 1)/10)
y <-  rnorm(n, a * x + b, sigma)

plot(x, y)



# -------------------------------------------------------------------------

likelihood <- function(param) {
  a <-  param[1]
  b <-  param[2]
  sd <- param[3]

  y_hat <- a*x + b
  singlelikelihoods <- dnorm(y, y_hat, sd = sd, log = T)
  sumll <- sum(singlelikelihoods)
  return(sumll)
}

predictSlope <- function(x) {
  return(likelihood(c(x, a, sigma)))
  }

slopelikelihoods <- sapply(seq(-3, 7, by=.05), predictSlope)
plot(seq(-3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")

prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior + bprior + sdprior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

# Metropolis MCMC -------------------------------------------------------------------------
# 1. starting at a random parameter value
# 2. Choosing a new parameter value close to the old value based on some probability density that is called the proposal function
# 3. Jumping to this new point with a probability p(new)/p(old), where p is the target function, and p > 1 means jumping as well


proposalfunction <- function(param){
  return(rnorm(param, mean = param, sd= c(0.1,0.5,0.3)))
}

## Metropolis MCMC algorithms
runMetropolis <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1, 3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal <- proposalfunction(chain[i,])

    probab <- exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] <- proposal
    } else {
      chain[i+1,] <- chain[i,]
    }
  }
  return(chain)
}

startvalue <- c(4,0,10)
chain <- runMetropolis(startvalue, 50000)

burnin <- 20000
acceptance <- 1 - mean(duplicated(chain[-(1:burnin),]))

# Again, working with the logarithms of the posterior might be a bit confusing at first, in particular when you look at the line where the acceptance probability is calculated (probab = exp(posterior(proposal) – posterior(chain[i,]))). To understand why we do this, note that p1/p2 = exp[log(p1)-log(p2)].

## Plot --------------------------------------------------------------------

hist(chain[-(1:burnIn),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = a, col="red" )

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = b, col="red" )

hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = sigma, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = a, col="red" )

plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = b, col="red" )

plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = sigma, col="red" )

# for comparison:
summary(lm(y~x))


# HMC (Hamilton Monte-Carlo) -------------------------------------------------------------------------
## inspired by Hamiltonian classical mechanics
# https://blogs.rstudio.com/ai/posts/2019-10-03-intro-to-hmc/
#

# U is a function that returns the potential energy given q
# grad_U returns the respective partial derivatives
# epsilon stepsize
# L number of leapfrog steps
# current_q current position


## kinetic energy is assumed to be sum(p^2/2) (mass == 1)
## Here gradient calculation grad_U is crucial

runHMC <- function (U, grad_U, epsilon, L, current_q) {
  q <- current_q
  # independent standard normal variates
  p <- rnorm(length(q), 0, 1)
  # Make a half step for momentum at the beginning
  current_p <- p
  # Alternate full steps for position and momentum
  p <- p - epsilon * grad_U(q) / 2
  for (i in 1:L) {
    # Make a full step for the position
    q <- q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i != L) p <- p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end
  p <- p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    return (q)  # accept
  } else {
    return (current_q)  # reject
  }
}


# NUTS (No U turn sampler) -------------------------------------------------------------------------

## https://blogs.rstudio.com/ai/posts/2019-10-03-intro-to-hmc/
##
##


stepLeapfrog = function(theta, r, eps, grad_f, M_diag){
  r_tilde <- r + 0.5 * eps * grad_f(theta)
  theta_tilde <- theta + eps * r_tilde / M_diag
  r_tilde <- r_tilde + 0.5 * eps * grad_f(theta_tilde)
  list(theta = theta_tilde, r = r_tilde)
}

getJointLogDensity = function(theta, r, f, M_diag){
  f(theta) - 0.5*sum(r**2 / M_diag)
}

## Find a reasonable epsilon
findEpsilon = function(theta, f, grad_f, M_diag, eps = 1, verbose = TRUE){
  r <- rnorm(length(theta), 0, sqrt(M_diag))
  proposed <- stepLeapfrog(theta, r, eps, grad_f, M_diag)
  log_ratio <- getJointLogDensity(proposed$theta, proposed$r, f, M_diag) - getJointLogDensity(theta, r, f, M_diag)
  alpha <- ifelse(exp(log_ratio) > 0.5, 1, -1)
  if(is.nan(alpha)) alpha <- -1
  count <- 1
  while(is.nan(log_ratio) || alpha * log_ratio > (-alpha)*log(2)){
    eps <- 2**alpha * eps
    proposed <- stepLeapfrog(theta, r, eps, grad_f, M_diag)
    log_ratio <- getJointLogDensity(proposed$theta, proposed$r, f, M_diag) - getJointLogDensity(theta, r, f, M_diag)
    count <- count + 1
    if(count > 100) {
      stop("Could not find reasonable epsilon in 100 iterations!")
    }
  }
  if(verbose) message("Reasonable epsilon = ", eps, " found after ", count, " steps")
  eps
}

checkNUTS = function(s, theta_plus, theta_minus, r_plus, r_minus){
  if(is.na(s)) return(0)
  condition1 <- crossprod(theta_plus - theta_minus, r_minus) >= 0
  condition2 <- crossprod(theta_plus - theta_minus, r_plus) >= 0
  s && condition1 && condition2
}

buildTree = function(theta, r, u, v, j, eps, theta0, r0, f, grad_f, M_diag, Delta_max = 1000){
  if(j == 0){
    proposed <- stepLeapfrog(theta, r, v*eps, grad_f, M_diag)
    theta <- proposed$theta
    r <- proposed$r
    log_prob <- getJointLogDensity(theta, r, f, M_diag)
    log_prob0 <- getJointLogDensity(theta0, r0, f, M_diag)
    n <- (log(u) <= log_prob)
    s <- (log(u) < Delta_max + log_prob)
    alpha <- min(1, exp(log_prob - log_prob0))
    if(is.nan(alpha)) stop()
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(theta_minus=theta, theta_plus=theta, theta=theta, r_minus=r,
                r_plus=r, s=s, n=n, alpha=alpha, n_alpha=1))
  } else{
    obj0 <- buildTree(theta, r, u, v, j-1, eps, theta0, r0, f, grad_f, M_diag)
    theta_minus <- obj0$theta_minus
    r_minus <- obj0$r_minus
    theta_plus <- obj0$theta_plus
    r_plus <- obj0$r_plus
    theta <- obj0$theta
    if(obj0$s == 1){
      if(v == -1){
        obj1 <- buildTree(obj0$theta_minus, obj0$r_minus, u, v, j-1, eps, theta0, r0, f, grad_f, M_diag)
        theta_minus <- obj1$theta_minus
        r_minus <- obj1$r_minus
      } else{
        obj1 <- buildTree(obj0$theta_plus, obj0$r_plus, u, v, j-1, eps, theta0, r0, f, grad_f, M_diag)
        theta_plus <- obj1$theta_plus
        r_plus <- obj1$r_plus
      }
      n <- obj0$n + obj1$n
      if(n != 0){
        prob <- obj1$n / n
        if(runif(1) < prob){
          theta <- obj1$theta
        }
      }
      s <- checkNUTS(obj1$s, theta_plus, theta_minus, r_plus, r_minus)
      alpha <- obj0$alpha + obj1$alpha
      n_alpha <- obj0$n_alpha + obj1$n_alpha

    } else{
      n <- obj0$n
      s <- obj0$s
      alpha <- obj0$alpha
      n_alpha <- obj0$n_alpha
    }
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(theta_minus=theta_minus, theta_plus=theta_plus, theta=theta,
                r_minus=r_minus, r_plus=r_plus, s=s, n=n, alpha=alpha, n_alpha=n_alpha))
  }
}




## https://github.com/kasparmartens/NUTS/blob/master/R/nuts.R

#' No-U-Turn sampler
#'
#' @param theta Initial value for the parameters
#' @param f log-likelihood function (up to a constant)
#' @param grad_f the gradient of the log-likelihood function
#' @param n_iter Number of MCMC iterations
#' @param M_diag Diagonal elements of the mass matrix in HMC. Defaults to ones.
#' @param M_adapt Parameter M_adapt in algorithm 6 in the NUTS paper
#' @param delta Target acceptance ratio, defaults to 0.5
#' @param max_treedepth Maximum depth of the binary trees constructed by NUTS
#' @param eps Starting guess for epsilon
#' @return Matrix with the trace of sampled parameters. Each mcmc iteration in rows and parameters in columns.
#' @export
runNUTS <- function(theta, f, grad_f, n_iter, M_diag = NULL, M_adapt = 50, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE){
  theta_trace <- matrix(0, n_iter, length(theta))
  par_list <- list(M_adapt = M_adapt)
  for(iter in 1:n_iter){
    nuts <- stepNUTS(theta, iter, f, grad_f, par_list, delta = delta, max_treedepth = max_treedepth, eps = eps, verbose = verbose)
    theta <- nuts$theta
    par_list <- nuts$pars
    theta_trace[iter, ] <- theta
  }
  theta_trace
}


## One step of the NUTS
stepNUTS <- function(theta, iter, f, grad_f, par_list, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE){
  kappa <- 0.75
  t0 <- 10
  gamma <- 0.05
  M_adapt <- par_list$M_adapt
  if(is.null(par_list$M_diag)){
    M_diag <- rep(1, length(theta))
  } else{
    M_diag <- par_list$M_diag
  }

  if(iter == 1){
    eps <- findEpsilon(theta, f, grad_f, M_diag, eps = eps, verbose = verbose)
    mu <- log(10*eps)
    H <- 0
    eps_bar <- 1
  } else{
    eps <- par_list$eps
    eps_bar <- par_list$eps_bar
    H <- par_list$H
    mu <- par_list$mu
  }

  r0 <- rnorm(length(theta), 0, sqrt(M_diag))
  u <- runif(1, 0, exp(f(theta) - 0.5 * sum(r0**2 / M_diag)))
  if(is.nan(u)){
    warning("NUTS: sampled slice u is NaN")
    u <- runif(1, 0, 1e5)
  }
  theta_minus <- theta
  theta_plus <- theta
  r_minus <- r0
  r_plus <- r0
  j=0
  n=1
  s=1
  if(iter > M_adapt){
    eps <- runif(1, 0.9*eps_bar, 1.1*eps_bar)
  }
  while(s == 1){
    # choose direction {-1, 1}
    direction <- sample(c(-1, 1), 1)
    if(direction == -1){
      temp <- buildTree(theta_minus, r_minus, u, direction, j, eps, theta, r0, f, grad_f, M_diag)
      theta_minus <- temp$theta_minus
      r_minus <- temp$r_minus
    } else{
      temp <- buildTree(theta_plus, r_plus, u, direction, j, eps, theta, r0, f, grad_f, M_diag)
      theta_plus <- temp$theta_plus
      r_plus <- temp$r_plus
    }
    if(is.nan(temp$s)) temp$s <- 0
    if(temp$s == 1){
      if(runif(1) < temp$n / n){
        theta <- temp$theta
      }
    }
    n <- n + temp$n
    s <- checkNUTS(temp$s, theta_plus, theta_minus, r_plus, r_minus)
    j <- j + 1
    if(j > max_treedepth){
      warning("NUTS: Reached max tree depth")
      break
    }
  }
  if(iter <= M_adapt){
    H <- (1 - 1/(iter + t0))*H + 1/(iter + t0) * (delta - temp$alpha / temp$n_alpha)
    log_eps <- mu - sqrt(iter)/gamma * H
    eps_bar <- exp(iter**(-kappa) * log_eps + (1 - iter**(-kappa)) * log(eps_bar))
    eps <- exp(log_eps)
  } else{
    eps <- eps_bar
  }

  return(list(theta = theta,
              pars = list(eps = eps, eps_bar = eps_bar, H = H, mu = mu, M_adapt = M_adapt, M_diag = M_diag)))
}

