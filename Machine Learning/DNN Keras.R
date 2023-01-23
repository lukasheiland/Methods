##########################################################################################
# Neural networks ---------------------------------------------------------------
##########################################################################################

# NN ---------------------------------------------------------------------
# Input layer: 'Predictors'
# Hidden layer: Aktivierungsschwellen, Verbindungen, werden optimiert
# Output layer: Response/ Target

# DNN: viele deep layers, die weit vom Input weg sind

library(keras)
library(tensorflow)

tf$enable_eager_execution()
use_session_with_seed(42,disable_parallel_cpu = FALSE)



# Classification ----------------------------------------------------------

plot(iris)

# Prepare data
summary(iris)
X <- scale(iris[,1:4])
Y <- iris[, 5]

# 1. Build model
library(keras)

# this is an object, and will be modified in the environment
model <- keras_model_sequential() # init empty model
model %>%
  layer_dense(units = 50L, input_shape = list(ncol(X)), activation = 'relu') %>% # 1. hidden layer
  layer_dense(units = 50L, activation = 'relu') %>% # no new input_shape necessary
  layer_dense(units = 3L, activation = 'softmax') # no. of output neurons (= no. of classes, [for regressions it's 1])
summary(model)                                           #

# 2. Compile model - define loss and optimizer
model %>%
  compile(loss = loss_categorical_crossentropy, optimizer = optimizer_adamax(0.001)) # crossentropy

# explicit dummy coding
Y_dummy <- keras::to_categorical(as.integer(Y)-1L, num_classes = 3L) # helper function for dummy coding for contrasts

# 3. Fit model
model_history <- model %>%
  fit(x = X, y = Y_dummy, epochs = 200L, shuffle = T)

plot(model_history)

model %>%
  evaluate(X, Y)

predictions = predict(model, X) # probabilities for each class
preds = apply(predictions, 1, which.max)


oldpar = par()
par(mfrow = c(1,2))
plot(iris$Sepal.Length, iris$Petal.Length, col = iris$Species, main = "Observed")
plot(iris$Sepal.Length, iris$Petal.Length, col = preds, main = "Predicted")
par(oldpar)

# Predict
predict(model, X)


oldpar = par()
par(mfrow = c(1,2))
plot(iris$Sepal.Length, iris$Petal.Length, col = iris$Species, main = "Observed")
plot(iris$Sepal.Length, iris$Petal.Length, col = preds, main = "Predicted")
par(oldpar)




# Linear regression ------------------------------------------------------------

# Prepare data
data <- na.omit(airquality)
X <- apply(as.matrix(data[,2:ncol(data)]), 2, scale)
Y <- data$Ozone

# 1. Build model
model2 <- keras_model_sequential() # init empty model
model2 %>%
  layer_dense(units = 50L, input_shape = list(ncol(X)), activation = 'relu') %>% # 1. hidden layer
  layer_dense(units = 50L) %>% # no new input_shape necessary
  layer_dense(units = 50L) %>%
  layer_dense(units = 1L, activation = 'linear') # no. of output neurons for regressions  = 1
summary(model2)                                  # 'lnear' = default means simply sums of weights, ('no activation function')

# 2. Compile model - define loss and optimizer
model2 %>%
  compile(loss = loss_mean_squared_error, optimizer = optimizer_adamax(0.1))
# 

# 3. Fit model
model_history <- model2 %>% fit(x = X, y = Y, epochs = 500L, shuffle = T)

plot(model_history)
prediction <- predict(model2, X)

plot(data$Ozone, prediction, main = "Predicted vs Observed")
abline(0, 1, col = "red")


# Convolutional NN ---------------------------------------------------------------------
#

