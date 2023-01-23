##########################################################################################
# GAMs - generalized additive models -------------------------------------------------
##########################################################################################

library(mgcv)

# GAM fit -----------------------------------------------------------------
# with smooth terns

fit <- gam(Ozone ~ s(Wind) + s(Temp), data = airquality) # exactly like the formula for a GLM except that smooth terms, s, te, ti and t2, can be added to the right hand side to specify that the linear predictor depends on smooth functions of predictors (or linear functionals of these).
par(mfrow = c(1,2))
plot(fit)


# Different interaction smooths -----------------------------------------------------------------
data(trees)
ct5 <- gam(Volume ~ te(Height, Girth, k=5), family=Gamma(link=log), data=trees)
ct5
vis.gam(ct5)

ct5 <- gam(Volume ~ s(Height, Girth, k=5), family=Gamma(link=log), data=trees)
ct5
vis.gam(ct5)

plot(ct5, too.far=0.15)

# s(Y) + s(X)
# s(Y, X) # multivariate smooth with interactions
# te(Y, X) # Tensor product smooths are for modeling responses to interactions of multiple inputs with different units (scales!).
# ti(Y, X) # The purpose of ti() is to separate the interactions from individual univariate effects.
## In general, the tensor product constructor implementing a ti() or te() term will reparameterize the marginal smooths so that the coefficients are interpretable as values of the smooth at evenly spaced covariate values (see section 5.6.2. Wood, 2017, Generalized Additive Models: An introduction with R)


# GAM catching variance with nuisance variables ---------------------------
# wenn Variablen nerven (nuisance variables) aber nicht interessieren: spline als spline fitten und der Effekt ist weg:
# e.g., confounders, detrending

# GAM is basically a GLM with splines:
fit2 <- gam(Ozone ~ Wind + s(Temp), data = airquality) # e.g. when interested in Wind:

summary(fit2)
par(mfrow = c(1,1))


# Mixed effects GAMs ------------------------------------------------------
## for GAM with lme4 syntax
# library(gamm4)

