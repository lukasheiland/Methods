##########################################################################################
# Random Forest ------------------------------------------------------------
##########################################################################################


# Library -----------------------------------------------------------------
library(randomForest)
library(forestFloor)
# devtools::install_github('araastat/reprtree')
library(reprtree)



# Continuous random forest -----------------------------------------------------------
fit_rf <-  randomForest(Ozone ~ .,
                        data = airquality[complete.cases(airquality),],
                        mtry = 3,
                        ntree = 500,
                        importance = T,
                        keep.inbag = T)

fit_rf
plot(fit_rf)

# Predictions
prediction <- predict(fit_rf, airquality)
plot(Ozone ~ Temp, data = airquality)
lines(airquality$Temp[order(airquality$Temp)], prediction[order(airquality$Temp)], col = 2)

## Variable importance
# often problematic: interpretation as effect size
randomForest::varImpPlot(fit_rf)

# Tree Visualization
reprtree:::plot.getTree(fit_rf) # some consensus tree
# forestFLoor(fit_rf)


# Problem: Kolinearität. RF kann, im Ggs zu lm mit richtigen Annahmen, keine Kausalität trennen.
# RF teilt die Varianz halb/halb auf.



# Classification ----------------------------------------------------------
classification <- randomForest(Species ~ .,
                               data = iris,
                               mtry = 3,
                               ntree = 500,
                               importance = T,
                               keep.inbag = T)
reprtree:::plot.getTree(classification)
plot(classification)
## see advanced statistics script for classification visualization methods!

