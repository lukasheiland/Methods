##########################################################################################
# Regularization with Keras ------------------------------------------------------------
##########################################################################################

# Library -----------------------------------------------------------------
library(keras)
library(tensorflow)

library(Metrics)

tf$enable_eager_execution()


# Data --------------------------------------------------------------------
data = airquality[complete.cases(airquality),]
X = scale(data[,-1])
Y = data$Ozone

load("Data/Titanic.RData")

test$survived <-  NA
data = rbind(train, test)
data[,-2] <-  missRanger::missRanger(data[,-2])

sub_data <-  data[, c(1, 2, 4, 5, 9)]

sub_data <-  cbind(sub_data, to_categorical(as.integer(sub_data$sex)-1L, 2L))
sub_data <-  sub_data[,-3]

train <-  sub_data[1:nrow(train),]
head(train)
X  <-  scale(train[,-2])
Y  <-  train[,2]
Y  <-  to_categorical(as.integer(Y), 2L)

indices <-  sample.int(nrow(X), 0.5*nrow(X))

X_train  <- X[indices,]
Y_train  <- Y[indices,]

X_test <- X[-indices,]
Y_test <- Y[-indices,]



# L1, L2 on linear model  ----------------------------------

#### 1. unconstrained
model <-  keras_model_sequential()
model %>%
  layer_dense(units = 200L, activation = "relu", input_shape = list(ncol(X_train))) %>%
  layer_dense(units = 200L, activation = "relu") %>%
  layer_dense(units = 2L, activation = "softmax")

summary(model)

model %>%
  compile(loss = loss_categorical_crossentropy, optimizer = optimizer_adam(0.01))
history <-
  model %>%
  fit(x = X_train, y = Y_train,
      validation_data = list(X_test, Y_test), epochs = 200L,
      batch_size = 32L, shuffle = TRUE)


plot(history) # without regularization


#### 2. L1 (Lasso)
lambda <-  0.01 # fiddle with lambda!
model <-  keras_model_sequential()
model %>%
  layer_dense(units = 200L, activation = "relu", input_shape = list(ncol(X_train)),
              kernel_regularizer = regularizer_l1(l = lambda)) %>%
  layer_dense(units = 200L, activation = "relu",
              kernel_regularizer = regularizer_l1(l = lambda)) %>%
  layer_dense(units = 2L, activation = "softmax",
              kernel_regularizer = regularizer_l1(l = lambda))
summary(model)

model %>%
  compile(loss = loss_categorical_crossentropy, optimizer = optimizer_adam(0.01))

history <- model %>%
  fit(x = X_train, y = Y_train,
      validation_data = list(X_test, Y_test), epochs = 200L,
      batch_size = 32L, shuffle = TRUE)

plot(history)

pred <- model %>%
  predict(X_test)

pred_classes <-  ifelse(pred[,2] < 0.5, 0, 1)

Metrics::auc(Y_test[,2], pred[,1])

image(t(model$weights[[1]]$numpy()))


#### 3. L2 (Ridge)

lambda <-  0.01 # fiddle with lambda!
model <-  keras_model_sequential()
model %>%
  layer_dense(units = 20L, activation = "relu", input_shape = list(ncol(X_train)),
              kernel_regularizer = regularizer_l2(l = lambda)) %>%
  layer_dense(units = 20L, activation = "relu",
              kernel_regularizer = regularizer_l2(l = lambda)) %>%
  layer_dense(units = 2L, activation = "softmax",
              kernel_regularizer = regularizer_l2(l = lambda))
summary(model)

model %>%
  compile(loss = loss_categorical_crossentropy, optimizer = optimizer_adam(0.01))

history <- model %>%
  fit(x = X_train, y = Y_train,
      validation_data = list(X_test, Y_test), epochs = 200L,
      batch_size = 32L, shuffle = TRUE)

plot(history)

pred <- model %>%
  predict(X_test)

pred_classes <-  ifelse(pred[,2] < 0.5, 0, 1)

Metrics::auc(Y_test[,2], pred[,1])



# L1, L2, elasticnet in TensorFlow ----------------------------------------

# - try to understand the following tensorflow core code
# - implement l1 and l2 in tensorflow core
# - implement elastic net in tensorflow core (lambda * ( (1-alpha)/2*l2 + alpha*l1 ))
data = airquality[complete.cases(airquality$Ozone) & complete.cases(airquality$Solar.R),]
X = scale(data[,-1])
Y = data$Ozone


W = tf$Variable(
  tf$constant(runif(ncol(X), -1, 1), shape = list(ncol(X), 1L), "float64")
)
B = tf$Variable(tf$constant(runif(1,-1,1), "float64"))

epochs = 200L
optimizer = tf$keras$optimizers$Adamax(1)

get_batch = function(batch_size = 32L){
  indices = sample.int(nrow(X), size = batch_size)
  return(list(bX = tf$constant(X[indices,], "float64"), bY = tf$constant(Y[indices], "float64", list(batch_size, 1L))))
}

steps = floor(nrow(X)/32) * epochs

zero = tf$constant(0.0, "float64")
one = tf$constant(1.0, "float64")
two = tf$constant(2.0, "float64")

l1_tf = function(W, B, lambda = tf$constant(1.0, "float64")) lambda * tf$reduce_mean(tf$abs(W) + tf$abs(B))
l2_tf = function(W, B, lambda = tf$constant(1.0, "float64")) lambda * tf$reduce_mean(tf$square(W) + tf$square(B))
elastic_tf = function(W, B, alpha = tf$constant(0.5, "float64"), lambda = tf$constant(1.0, "float64")) {
  lambda * ((one - alpha)/two * l2_tf(W, B, one) + alpha* l1_tf(W, B, one))
}

for(i in 1:steps){
  batch = get_batch()
  bX = batch$bX
  bY = batch$bY

  with(tf$GradientTape() %as% tape,{
    pred = tf$matmul(bX, W) + B
    loss = tf$reduce_mean(tf$keras$losses$mean_squared_error(bY, pred)) + elastic_tf(W, B)
  })

  gradients = tape$gradient(loss, c(W, B))

  optimizer$apply_gradients(purrr::transpose(list(gradients, c(W, B))))

  if(i %% floor(nrow(X)/32)*20 == 0) cat("Loss: ", loss$numpy(), "\n")

}

