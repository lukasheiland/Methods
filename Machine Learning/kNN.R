##########################################################################################
# kNN ------------------------------------------------------------
##########################################################################################

# Library -----------------------------------------------------------------
library(kknn)


# Data for labelling a new point ------------------------------------------
X = scale(iris[, 1:4])
Y = iris[, 5]
plot(X[-100, 1], X[-100, 3], col = Y)
points(X[100, 1],
       X[100, 3],
       col = "blue",
       pch = 18,
       cex = 1.3)
# Find the k nearest neighbors and for instance choose the label by voting
# Disadvantage, we have to calculate a distance matrix for all points
dist_X = dist(X)



# Split:
## Scaling is important with distances!
data = iris
data[, 1:4] = apply(data[, 1:4], 2, scale)

indices = sample.int(nrow(data), 0.7 * nrow(data))
train = data[indices, ]
test = data[-indices, ]


knn <- kknn(Species ~ ., train = train, test = test)
summary(knn)
table(test$Species, fitted(knn))

oldpar = par()
par(mfrow = c(1, 2))
plot(test$Sepal.Length,
     test$Petal.Length,
     col =  predict(knn),
     main = "predicted")
plot(test$Sepal.Length,
     test$Petal.Length,
     col =  test$Species,
     main = "observed")
par(oldpar)

## Exercise 1 - airquality
## - split airquality in train and test
## - fit knn
## - predict test values
## - visualize result (see previous exercise for trees)

data = airquality[complete.cases(airquality$Ozone) &
                    complete.cases(airquality$Solar.R), ]
data[, 2:6] = apply(data[, 2:6], 2, scale)

indices = sample.int(nrow(data), 0.7 * nrow(data))
train = data[indices, ]
test = data[-indices, ]

knn = kknn(Ozone ~ ., train = train, test = test)
pred = predict(knn)
plot(test$Temp, test$Ozone)
lines(test$Temp[order(test$Temp)], pred[order(test$Temp)], col = "red")






## Exercise - 2 airquality
## - split airquality in train and test
## - fit knn
## - predict test values
## - visualize result (see previous exercise for trees)

data = airquality[complete.cases(airquality$Ozone) &
                    complete.cases(airquality$Solar.R), ]
data[, 2:6] = apply(data[, 2:6], 2, scale)

indices = sample.int(nrow(data), 0.7 * nrow(data))
train = data[indices, ]
test = data[-indices, ]

sm = svm(Ozone ~ ., data = train)
pred = predict(sm, newdata = test)
plot(test$Temp, test$Ozone)
lines(test$Temp[order(test$Temp)], pred[order(test$Temp)], col = "red")

