##########################################################################################
# PCA ---------------------------------------------------------------
##########################################################################################

# Library -----------------------------------------------------------------


# Simulate data --------------------------------------------------------------
n_features <- 5
feature <- 1:n_features
n_obs <- 30

optima <- 1:n_features*n_obs/n_features
M <- sapply(optima, function(x) dnorm(1:n_obs, x, n_obs/n_features))
matplot(M) # Features along an ordered "gradient" of observations 1:n_obs

optima2 <- runif(n_features, 0, n_obs)
M2 <- sapply(optima2, function(x) rnorm(n_obs, dnorm(1:n_obs, x, n_obs*10) * 10, 0.0001))
matplot(M2) # Features along an ordered "gradient" of observations 1:n_obs

M <- M2

getPCA <- function(M){
  ## 0. scale data, actually only centering is necessary (by subtracting the mean from data)
  M <- scale(M)

  ## 1. compute covariance matrix
  C <- cov(M)

  ## 2. get Eigenvectors and corresponding eigenvalues
  # Eigenvector properties:
  # - def: vector that multiplied with the matrix yields the same vector scaled with a scalar (the eigenvalue) M*v = s*v
  # - scaling the vector will not change its direction
  # - all eigenvectors of a matrix orthogonal
  Evector <- eigen(C)$vectors # n_features (unit!) eigenvectors
  evalue <- eigen(C, only.values = T)$values # n_features eigenvalues, already ordered decendingly!

  ## Eigenvectors of the covariance matrix characterize the axes of variation between the data
  ## and the corresponding eigenvalues are equivalent to the sd!
  plot(M[,1], M[,2], col = 1:n_obs)
  lines(c(0, Evector[1,1]), c(0, Evector[2,1])) # 1st eigenvector
  lines(c(0, Evector[1,2]), c(0, Evector[2,2])) # 2nd eigenvector

  ## For some reason princomp() sometimes reverts the direction of some Principal component:
  # Evector[,2] <- -Evector[,2] # ???

  ## 3. New data, transformed ("rotated") along the axes of Eigenvectors
  # What this is saying is that the first row of Scores is the sum of rows of M weighted by the first eigenvector.
  Scores <- t( t(Evector) %*% t(M) )
  # plot(Scores[,1], Scores[,2],  col = 1:n_obs)

  fakeprincomp <- list(sdev = evalue, # evalue/sum(evalue) #?
                       loadings = Evector,
                       center = attr(M, "scaled:center"),
                       scale = attr(M, "scaled:scale"),
                       n.obs = nrow(M),
                       scores = Scores,
                       call = 'fuckyou')
  class(fakeprincomp) <- "princomp"
  return(fakeprincomp)
}


biplot(princomp(M))
biplot(getPCA(M))


## Compare the resulting eigenvectors.
# princomp(M)$loadings
# getPCA(M)$loadings


# Other methods --------------------------------------------------------------------
library(vegan)
nmds <- vegan::metaMDS(M)

## We know the more or less true "distance" of features, because we have ordered them along 1:n_obs
gradient <- scale(1:n_obs)
gradient_error <- rnorm(gradient, gradient, 5)

ordiplot(getPCA(M), type = "text")

# M <- M2
ordisurf(getPCA(M), gradient)
ordisurf(metaMDS(M), gradient)


