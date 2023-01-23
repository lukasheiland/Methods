library(glmmTMB)
library(mvtnorm)
library(tictoc)

## Functions ---------------------------------------------------------------

# function to create the interaction matrix
# x=number of species
create.sigma = function(x) {
  sigma1 = matrix(runif(x^2, -1,1), x, x)
  sigma2 = sigma1 %*% t(sigma1)
  sigma3 = cov2cor(sigma2)
  return(sigma3)
}

# function to arrange data structure for glmer/glmmTMB
# x=number of sites
# y=number of species
arrange.data = function(x, y, occ) {
  data = data.frame(site=rep(1:x, y),
                    env1=env.var[ ,1],
                    env2=env.var[ ,2],
                    env3=env.var[ ,3],
                    genus=factor(rep(1:y, each = x)),
                    occurrence=c(occ))
  return(data)
}




# Simulation Variables ----------------------------------------------------

n_e = 3 # number of environmental variables
n_s = 3 # number of species
n = 300*n_s # number of sites

n_se = n_s * n_e # number of regression parameters to fit (=value for the number of columns of the recov_weights matrix)
n_int = (n_s^2-n_s)/2 # number of interaction parameters to fit (= number of columns of the recov_int matrix), (n*n-n)/2 parameters

# Environment variables
env.var = matrix(runif(n_e*n, -0.5, 0.5), n, n_e)

# Species weights
species_weights = matrix(runif(n_s*n_e, -1, 1), n_e, n_s)

# occurrence (without species-species interaction)
y = env.var %*% species_weights

# Interaction matrix
sigma_interact = create.sigma(n_s)




# GLMM --------------------------------------------------------------------
iter = 100 # number of iterations

recov.int = matrix(ncol=iter, nrow=n_int) # empty matrix for recovered interaction parameters
recov.weights = matrix(ncol=iter, nrow=n_se) # empty matrix for recovered species weights (regression parameters)

#unparallelized
tic()
for (i in 1:iter){
  occurrence_logit_spp_count = apply(exp((t(apply(y, 1, function(x) rmvnorm(1, mean = x, sigma=sigma_interact))))), 1:2, function(i) rpois(1,  i))
  data = arrange.data(n, n_s, occurrence_logit_spp_count)
  fit.glmmTMB = glmmTMB(occurrence ~ 0 + genus + env1:genus + env2:genus + env3:genus + (0 + genus | site), family=poisson, data=data)
    # jede Site bekommt einen Genus random slope, für jede Site also, Quasi observation-level random effect
    # R fittet Korrelation für Faktoren mit (ursprünglich gedacht für Intercept/slope-Korrelation)
  
  fit.cov = VarCorr(fit.glmmTMB)
  fit.int = attributes(fit.cov$cond$site)$correlation
  
  diag(sdRE) <- attributes(fit.cov$cond$site)$stdev
  diag(sdRE) %*% fit.int %*% diag(sdRE)
  
  l = 0
  for (j in 1:ncol(fit.int)){
    for (k in 1:nrow(fit.int)){
      if (lower.tri(fit.int)[k,j] == TRUE){
        l = l+1
        recov.int[l,i] = fit.int[k,j]
      }
    }
  }
  fit.fixef = glmmTMB::fixef(fit.glmmTMB)$cond
  fit.weights = fit.fixef[-1:-n_s]
  for (j in 1:n_se){
    recov.weights[j,i] = fit.weights[j]
  }
  if(i==iter) toc()
}



# Plotting ----------------------------------------------------------------


# plot interaction parameters
#true.int = the real interaction values (sigma_interact) in different structure
par(mfcol=c(n_s,n_s))
true.int = matrix(nrow=n_int, ncol=1)
c=0
for (j in 1:ncol(sigma_interact)){
  for (k in 1:nrow(sigma_interact)){
    if (lower.tri(sigma_interact)[k,j] == TRUE){
      c = c+1
      true.int[c,1] = sigma_interact[k,j]
      plot(density(recov.int[c,]), main="")
      abline(v=true.int[c,])
    }
    else {
      plot.new()
    }
  }
}


# plot regression parameters
par(mfrow=c(n_e,n_s), mar=c(2,2,2,2))
for (i in 1:n_se){
  plot(density(recov.weights[i,]), main="")
  abline(v=c(t(species_weights))[i])
}
