##########################################################################################
# Cross validation ---------------------------------------------------------------
##########################################################################################
# Daten in Folds einteilen,
# Modell auf alle Folds außer einem fitten, auf ausgelassenen fold predicten und Fehler messen.
# -> Selektion oder Quatifikation der Prediktionsqualität


# Nested cross validation ----------------------------------------------------------
## inner: optimize model for prediction
## outer: optimize hyperparameters (like thresholds)
## there must be two nested loops, and another exclusively reserved testing set for each run of the outer loop.
# { training outer [                  training  inner                     ] [ testing inner ]}    [ testing outer ]

# Data ---------------------------------------------------------------------
AQ <- airquality[complete.cases(airquality),]
n_obs <- nrow(AQ)


# Ten-fold cross validation -----------------------------------------------
n_folds <- 10
squarederror_l <- 0
squarederror_q <- 0

for (i in 1:n_folds) {
  fold <- sample.int(n_folds, size = n_obs, replace = T) #!
  traindata <- AQ[fold != i,]
  testdata <- AQ[fold == i,]

  fit_l <- lm(Ozone ~ Temp, data = traindata)
  fit_q <- lm(Ozone ~ poly(Temp^2, degree = 2) + Solar.R, data = traindata)

  squarederror_l <- squarederror_l + sum((predict(fit_l, newdata = testdata) - testdata$Ozone)^2)
  squarederror_q <- squarederror_q + sum((predict(fit_q, newdata = testdata) - testdata$Ozone)^2)
}

RMSE_l <- sqrt(squarederror_l/n_obs) # Root-mean-square deviation/error. sqrt((mean(res^2)))
RMSE_q <- sqrt(squarederror_q/n_obs) # Root-mean-square deviation/error. sqrt((mean(res^2)))

RMSE_l > RMSE_q


# Leave one out cross validation -----------------------------------------------------------

n_folds <- n_obs #!
squarederror_l <- 0
squarederror_q <- 0

for (i in 1:n_folds) {
  fold <- sample.int(n_folds, size = n_obs, replace = F) #!
  traindata <- AQ[fold != i,]
  testdata <- AQ[fold == i,]

  fit_l <- lm(Ozone ~ Temp, data = traindata)
  fit_q <- lm(Ozone ~ poly(Temp^2, degree = 2) + Solar.R, data = traindata)

  squarederror_l <- squarederror_l + sum((predict(fit_l, newdata = testdata) - testdata$Ozone)^2)
  squarederror_q <- squarederror_q + sum((predict(fit_q, newdata = testdata) - testdata$Ozone)^2)
}

RMSE_l <- sqrt(squarederror_l/n_obs) # Root-mean-square deviation/error. sqrt((mean(res^2)))
RMSE_q <- sqrt(squarederror_q/n_obs) # Root-mean-square deviation/error. sqrt((mean(res^2)))

RMSE_l > RMSE_q

plot(Ozone ~ Temp, data = AQ)
points(predict(fit_l), col = 2)
points(predict(fit_q), col = 3)

# ROC curves ----------------------------------------------------------
# library(pROC)
# ?pROC
