##########################################################################################
# Neural networks ---------------------------------------------------------------
##########################################################################################

library(tensorflow)
exists("tf")
tf$enable_eager_execution()


# TensorFlow data structures
a <- tf$constant(5)
b <- tf$constant(10)
result <- tf$add(a, b)
class(result)

result$numpy() # get an R object back


r_ma <- matrix(runif(100), 10, 10)
ma <- tf$constant(r_ma, dtype = 'float32')

## Data types - good practise with R-TF
tf$reshape(ma, shape = list(100L, 1L))

## tensorflow arguments require also exact/explicit data types:


# Exercise 1 - TensorFlow
## exploring tf, tf$math and tf$linalg
## translate R into TF operations
set.seed(42)
x = rnorm(100)
max(x)
tf$math$reduce_max(x)$numpy()

min(x)
tf$math$reduce_min(x)$numpy()

mean(x)
tf$math$reduce_mean(x)$numpy()
tf$reduce_mean(x)$numpy()

which.max(x)
tf$argmax(x)

which.min(x)
tf$argmin(x)

order(x)
tf$argsort(x)


## Bonus
m = matrix(runif(9), 3, 3)
solve(m)


diag(m)


diag(diag(m))


eigen(m)


det(m)


# lm in Tensorflow --------------------------------------------------------
x <- tf$constant(matrix(runif(1000,-1,1), ncol = 10, nrow = 100), tf$float32)
w <-  tf$constant(matrix(runif(10, 0, 2), ncol = 1), tf$float32)

tf$numpy()

tfp = reticulate::import("tensorflow_probability")
y <- tf$matmul(x, w) + tfp$distributions$Normal(0., scale = 0.5)$sample(100)

w_query <-  tf$Variable(matrix(0, ncol = 1, nrow = 10), dtype = tf$float32)

for(steps in 1:100){
  with(tf$GradientTape(persistent = T) %as% tape,
       {
         y_hat = tf$matmul(x, w_query)
         res = tf$sqrt(tf$square(y - y_hat))
         loss = tf$reduce_mean(res)
       })

  gradients <- tape$gradient(loss, w_query)

  optimizer <- tf$train$AdamOptimizer(learning_rate = 0.01)
  up <- list(list(gradients, w_query))
  optimizer$apply_gradients(up, tf$train$get_or_create_global_step())
}
cbind(w$numpy(), w_query$numpy())
w_query

