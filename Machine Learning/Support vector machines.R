##########################################################################################
# Support vector machines SVM ---------------------------------------------------------------
##########################################################################################
# Hyperplanes dividing Data.
# (Finding a hyperplane which maximes the margin, the distance between the plane and the points close to it.)

# Kernel-trick: SVM works fine as long as the problem is linear seperable.
# If not we have to map the problem into a space in which it is linear seperable. That is called the kernel trick


# Library -----------------------------------------------------------------
library(e1071)


# Data --------------------------------------------------------------------
Data = iris
Data[, 1:4] = apply(Data[, 1:4], 2, scale)

indices = sample.int(nrow(Data), 0.7 * nrow(Data))
train = Data[indices,]
test = Data[-indices,]


sm = svm(Species ~ ., data = train, kernel = "linear")
pred = predict(sm, newdata = test)

oldpar = par()
par(mfrow = c(1, 2))
plot(test$Sepal.Length,
     test$Petal.Length,
     col =  pred,
     main = "predicted")
plot(test$Sepal.Length,
     test$Petal.Length,
     col =  test$Species,
     main = "observed")
par(oldpar)

mean(pred == test$Species)
