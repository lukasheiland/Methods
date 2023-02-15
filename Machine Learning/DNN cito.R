library(torch)
library(cito)

n <- 1000
D <- data.frame(x = runif(n, -1, 1))
D$y_hat <- D$x * 0.6 + 2
D$y <- rnorm(D$x, D$y_hat, 1)
D$train <- as.logical(rbinom(D$x, 1, 0.2))

# Build and train  Network
fit <- dnn(y ~ x + 1, data = D[D$train,],
           hidden = c(10L, 10L),
           # alpha = 0,
           # lambda = 0.2,
           loss = "mse", activation = "relu",
           epochs = 40)
fit$weights

# continue training for another 32 epochs
# fit <- continue_training(fit)

# Use model on validation set
predictions <- predict(fit, D[!D$train,])

# Scatterplot
# plot(D[!D$train, "y"], predictions)
plot(D[!D$train, "y_hat"], predictions)
# plot(D[!D$train, "x"], predictions)
# MAE
mean(abs(predictions - D[!D$train, "y"])) ## overfitting

