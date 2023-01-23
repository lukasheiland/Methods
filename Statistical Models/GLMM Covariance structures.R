# Library ------------------------------------------------------------------
library(nlme)

library(glmmTMB)
library(DHARMa)

# Source ------------------------------------------------------------------
source('~/Documents/Studium/Methods/Helpers.R') # for range01()

##########################################################################################
# Covariance ---------------------------------------------------------------
##########################################################################################

# lm assumptions for residuals: iid Normal =
#   - independent[!!! (but groups treated by random effects)],
#   - identical [variance homogeneity treated by]
#   - Normal [distribution treated by GLM].

# uncorrected correlation leads to
#   - increased type I-error (p-value),
#   - wrong CI,
#  (- maybe bias when sampling unbalanced)

## Regression with correlation
# yvector ~ a * xvector + b + Errormatrix
# where correlation in Errormatrix is a function(a, distance)

# in glm: errordistribution + link[y ~ ax + Errormatrix]
# overdispersion through linear term (within link function) already having Errormatrix (only this way correlation can be included) + then another errordistribution

# A TREND does not satisfy the assumptions of a covariance autocorrelation Structure.
# fit time as predictor for detrending!
# for flexibility and speed just use a GAM for detrending!
# this will catch a lot of correlation.



# RSA ---------------------------------------------------------------------
# RSA = residual spatial autocorrelation (IN RESIDUALS, not in response!)
# response is almost always spatially autocorrelated, which is fine,
# but should not be in residuals!
# Spatial auto-correlation leads to  pseudo-replicates -> residuals all in the same direction.
# 'Four points measured in spatial vicinity might only be 1.7 points in terms of information.'
#
# Covariance Model: exponential

# Temporal autocorrelation -------------------------------------------------
# Too temporally close: pseudoreplication,
#
# Covariance Model: AR1: error = independent term + … = quadratic de


# Phylogenetic autocorrelation -------------------------------------------------
#
# Covariance Model: Brownian motion, proportional sqrt(time since phylogenetic divergence) :: PGLS

# Distance autocorrelation -------------------------------------------------
# works with distance measures in general! e.g. degree of frienbdship


# Temporal analysis ----------------------------------------------------------------
library(nlme)
library(lattice)
xyplot(follicles ~ Time | Mare, type = 'b', data = Ovary)

fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary) # just a normal lm
xyplot(residuals(fm1) ~ Time | Mare, type = 'b', data = Ovary) # there is temporal correlation, even trent
# plot(rnorm(50), type = 'b') # this would be random

## Diagnosis
# 1. Variogram: indipendence compared to distance between data points
plot(Variogram(fm1, form = ~ Time | Mare))
# they are getting more independent after 0.6
# no autocorrelation, constant 1 over distance.

# 2. formal hypothesis test woult be: dwtest, e.g.
# DHARMa::testTemporalAutocorrelation(simulateResiduals(model), time = numerictime)

# -> there is autocorrelation
# -> but there is also a trend!

## Improvement:
# 1. detrending, include as interaction to give each Mare her own trend:
fm_detrended <- gls(follicles ~ Time * Mare + sin(2*pi*Time) + cos(2*pi*Time), Ovary)
xyplot(residuals(fm_detrended) ~ Time | Mare, type = 'b', data = Ovary) # there is temporal correlation, even trent
plot(Variogram(fm_detrended, form = ~ Time | Mare))

# 2. Add correlation structure
fm_full <- gls(follicles ~ Time * Mare + sin(2*pi*Time) + cos(2*pi*Time),
               Ovary,
               correlation = corAR1(form = ~ 1 | Mare))
plot(Variogram(fm_full, form = ~ Time | Mare)) # Correlation structure THIS DOES NOT CHANGE THE VARIOGRAM, NEITHER residuals (NEITHER FORMAL TESTS)
# Fitting correlation structure does not change the residuals! (just makes the p-values 'calibrated')

# compare changed p-values
summary(fm_detrended)
summary(fm_full)


# Another one with gls ------------------------------------------------------------

str(Dialyzer)
fm1Dial.gls <- gls(rate ~ poly(pressure, 4) * QB, Dialyzer)
plot(fm1Dial.gls)
fm2Dial.gls <- update(fm1Dial.gls, weights = varPower(form = ~ pressure))
plot(fm2Dial.gls)
fm3Dial.gls <- update(fm2Dial.gls, corr = corAR1(form = ~ 1 | Subject)) # value is only a start value!
summary(fm3Dial.gls)


# A spatial example -------------------------------------------------------
# the measured thickness of coal seams at different coordinates
# soil measuring the soil quality
spdata <- read.table("https://stats.idre.ucla.edu/stat/r/faq/thick.csv", header = T, sep = ",")
plot(thick ~ soil, data = spdata)

summary(fit_lm <- lm(thick ~ soil, data = spdata))
abline(fit_lm)

# correlation in response is no problem!
# problematic is only autocorrelation in the residuals

plot(north ~ east,
     # cex = residuals(fit_lm) - min(residuals(fit_lm)),
     data = spdata,
     col = rgb(colorRamp(c(1, 'white', 2))(range01(residuals(fit_lm))),
               maxColorValue = 255)) # just scaling
  # btw: divergent scale rechts und links von 0! -> über weiß
# -> there seems to be a quadratic dependecy on north!


# Achtung: wenn Fleckengröße im Raum verschieden groß ->
# careful all structures assume homogeneity of the correlation (same correlation in space)

# Total spatial correlation
library(gstat)
vario_tanndir <- variogram(residuals(fit_lm) ~ 1, loc = ~ east + north, data = spdata)
plot(vario_tanndir)

## Moran's I
# standard formal test
Dist <- as.matrix(dist(cbind(spdata$east, spdata$north)))
Dist_inv <- 1/Dist # Eine Art power-law Modell, entspricht grob Distanzen auf Grund Brownian Motion, Dist für Moran.I kann beliebig gesetzt werden
diag(Dist_inv) <- 0

ape::Moran.I(residuals(fit_lm), Dist_inv) # tests forr weight

## remove trend
library(mgcv)
summary(fit_gam <- gam(thick ~ soil + te(east, north), data = spdata)) # discussion worthy: Raum absorbiert Prädiktoren. 'Prädiktoren sollten besser sein als Raum.'
plot(fit_gam)
vario2 <- variogram(residuals(fit_gam) ~ 1, loc = ~ east + north, data = spdata)
plot(vario2)

summary(lm(thick ~ soil + north + I(north^2), data = spdata))

# + Moran's I

# detrending is enough (no RSA left)! But here is spatial model:

## just for syntax: autocor
# standard isotrop spatial model: exponential decay with space

correlation = corExp(form = ~ east + north)
fit_auto <- gls(thick ~ soil , correlation = corExp(form = ~ east + north), data = spdata)
summary(fit_auto)
plot(variogram(residuals(fit_auto) ~ 1, loc = ~ east + north, data = spdata)) # bringt ja nix!
ape::Moran.I(residuals(fit_auto), Dist_inv)



# Another spatial example -------------------------------------------------
# ordinary least square = lm vs. gls
S <- read.table(header = T, file = '/Users/heiland/Dropbox/Studium/Exercise/Statistikeinführung/2019 Advanced course/Models/snouterdata.txt')
str(S)


fit <- gls(snouter1.1 ~ rain + djungle , correlation = corExp(form = ~ X + Y), data = S)
summary(fit)
image(xtabs(residuals(fit) ~ X+Y, data = S)) # immer die Residuen vor der correlation. Nur sind pvalues verschieden.

####
### Spatial data

# 1. Ganzes Modell fitten
# 2. Resiuden checken auf räumliche Korrelation
# 3. Trends rausrechnen
# 4. Wieder Residuen checken
# 5. Wenn notwendig Korrelationsstruktur einbauen

summary(nicemodel <- glmmTMB(snouter2.1 ~ rain + djungle, family = binomial, data = S))
plot(res <- simulateResiduals(nicemodel))
testSpatialAutocorrelation(res, x = S$X, y = S$Y) #

S$ID <- 1:nrow(S)
S$point <- numFactor(S$X, S$Y) # glmmTMB probably checks for the same factor level!

summary(nicermodel <- glmmTMB(snouter2.1 ~ rain + djungle + exp(point + 0 | ID), family = binomial, data = S))
plot(res <- simulateResiduals(nicemodel))
testSpatialAutocorrelation(res, x = S$X, y = S$Y)

## Phylogenetic data
gls(trait1 ~ trait2, correlation = corBrownian(phy = phylotree)) # trade off





# Sim ---------------------------------------------------------------------
library(nlme)
set.seed(0)

n_obs <- 400
n_x <- 2

X <- matrix(runif(n_obs*n_x, -1, 1), ncol = 2)

Coord <- data.frame(lon = runif(n_obs, -10, 10), lat = runif(n_obs, 0, 20))
plot(Coord)

#### Correlation structure
range <- 1.5 ## range parameter determines the exponential spatial deoendency
# curve(exp(-x/range), 0, 10)
Distance <- as.matrix(dist(Coord, method = 'euclidean'))
Sigma <- exp(-Distance/range)

## Ensure to be positive definite, but should be anyway
# diag(Sigma) <- 0
# Sigma <- sfsmisc::posdefify(Sigma) # Finds a close positive definite matrix
# Sigma <- cov2cor(Sigma)

## Equivalent to
# cs1Exp <- corExp(0.5, form = ~ lon + lat)
# cs1Exp <- Initialize(cs1Exp, Coord)
# Sigma <- corMatrix(cs1Exp)

K <- MASS::mvrnorm(n = 1, mu = rep(0, n_obs), Sigma = Sigma)

beta = c(2, -3)
y_hat <- X %*% beta

## Here, only K will be added to
## Adding another non-spatial error term would still be fitted within K in nlme::gls. Anywhere?
# sigma <- 2
# e <- rnorm(y_hat, 0, sd = sigma)
y <- y_hat + K
D <- data.frame(id = 1:n_obs, y, x = X, Coord)

plot(Coord, cex = range01(getData()$y)*2) # bigger circles clump together

#### Fit
summary(fit_gls0 <- gls(y ~ x.1 + x.2, data = D))
summary(fit_gls <- gls(y ~ x.1 + x.2, correlation = corExp(form = ~ lon + lat), D))


# fit_gls$sigma
# coef(fit_gls$modelStruct$corStruct, unconstrained = F)

# remember: this does not change the residuals compared to the model without correlation structure, it only makes the p-values calibrated!

#### glmmTMB
D$point <- numFactor(D$lon, D$lat) # glmmTMB probably checks for the same factor level!
summary(fit_tmb <- glmmTMB(y ~ x.1 + x.2 + exp(point + 0 | id), data = D))
VarCorr(fit_tmb)

fit_tmb$fit$par
summary(fit_tmb)
# plot(res <- simulateResiduals(fit_tmb0))
testSpatialAutocorrelation(res, x = D$coord.1, y = D$coord.2)

## BEWARE:
# When resampling, e.g. cross validation, in data with correlation structure,
# block the folds according to structure! (Otherwise they won't be random!)
# E.g. in spatial data: folds stratified spatially. Phylogenetically: folds = clades.
