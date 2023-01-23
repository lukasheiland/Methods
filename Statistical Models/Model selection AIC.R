##########################################################################################
# Model selection ---------------------------------------------------------------
##########################################################################################

library(MASS)

selectedmodel <- MASS::stepAIC(completemodel)
summary(selectedmodel)
plot(selectedmodel)

# for extended AIC seletion: library(MuMIn)