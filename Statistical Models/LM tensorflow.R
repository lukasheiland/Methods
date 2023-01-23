library(tensorflow)

tf$enable_eager_execution()


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
