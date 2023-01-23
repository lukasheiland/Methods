##########################################################################################
# Multi Model Inference  ---------------------------------------------------------------
##########################################################################################


# Library -----------------------------------------------------------------
library(MuMIn)


# MuMIn -------------------------------------------------------------------------
options(na.action = "na.fail") #  change the default "na.omit" to prevent models
#  from being fitted to different datasets in
#  case of missing values.


# Fits --------------------------------------------------------------------

fm1 <- lm(y ~ ., data = Cement)
ms1 <- dredge(fm1) # Generate a model selection table of models with combinations (subsets) of fixed effect terms in the global model, with optional rules for model inclusion.

# Visualize the model selection table:
plot(ms1, labAsExpr = TRUE)
ms1
model.avg(ms1, subset = delta < 4)
## There is no statistical justification to estimate the Std. Error like it is done here:
## Here it is just the SD of estimates.

confset.95p <- get.models(ms1, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)



# Revert MuMIn option -----------------------------------------------------
options(na.action = "na.omit") #  change the default "na.omit" to prevent models

