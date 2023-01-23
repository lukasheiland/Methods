##########################################################################################
# Regression trees ------------------------------------------------------------
##########################################################################################

# Library -----------------------------------------------------------------
library(tree)
library(rpart)


# Data --------------------------------------------------------------------
data <- airquality[complete.cases(airquality$Ozone) &
                    complete.cases(airquality$Solar.R), ]


# Do the thing ------------------------------------------------------------
rt <- tree(Ozone ~ .,
          data = data,
          control = tree.control(mincut = 30, nobs = nrow(data)))
plot(rt)
text(rt)

pred <- predict(rt, data)
plot(data$Temp, data$Ozone)
lines(data$Temp[order(data$Temp)], pred[order(data$Temp)], col = "red")

sqrt(mean((data$Ozone - pred) ^ 2))

## Exercise 1
## - change the mincut control options. What happens?
## - compare rmse

rt <- tree(Ozone ~ .,
          data = data,
          control = tree.control(mincut = 1L, nobs = nrow(data)))
plot(rt)
text(rt)

pred <- predict(rt, data)
plot(data$Temp, data$Ozone)
lines(data$Temp[order(data$Temp)], pred[order(data$Temp)], col = "red")


sqrt(mean((data$Ozone - pred) ^ 2))



rt <- tree(Ozone ~ .,
           data = data,
           control = tree.control(mincut = 50L, nobs = nrow(data)))
plot(rt)
text(rt)

pred <- predict(rt, data)
plot(data$Temp, data$Ozone)
lines(data$Temp[order(data$Temp)], pred[order(data$Temp)], col = "red")
sqrt(mean((data$Ozone - pred) ^ 2))

