library(torch)
library(luz)
library(torchvision)
library(torchdatasets)

## Example from https://torch.mlverse.org/start/guess_the_correlation/


# Img functions ------------------------------------------------------------

add_channel_dim <- function(img) {
  img$unsqueeze(1)
}

crop_axes <- function(img) {
  transform_crop(img, top = 0, left = 21, height = 131, width = 130)
}


# Data --------------------------------------------------------------------
## divide the same downloaded data into train, validation and test set

train_indices <- 1:10000
val_indices <- 10001:15000
test_indices <- 15001:20000

## download to root
root <- file.path(tempdir(), "correlation")
train_ds <- guess_the_correlation_dataset(
  # where to unpack
  root = root,
  # additional preprocessing 
  transform = function(img) crop_axes(img) %>% add_channel_dim(),
  # don't take all data, but just the indices we pass in
  indexes = train_indices,
  download = TRUE
)

valid_ds <- guess_the_correlation_dataset(
  root = root,
  transform = function(img) crop_axes(img) %>% add_channel_dim(),
  indexes = val_indices,
  download = FALSE
)

test_ds <- guess_the_correlation_dataset(
  root = root,
  transform = function(img) crop_axes(img) %>% add_channel_dim(),
  indexes = test_indices,
  download = FALSE
)

## How many items in each set?
length(train_ds) 
length(valid_ds)
length(test_ds)

train_ds[1] ## one observation. An image with a label "id". Torch tensor list format


# Divide into batches. ----------------------------------------------------
##  Train with batches, not with all data. Memory reasonst but also convergence
train_dl <- dataloader(train_ds, batch_size = 64, shuffle = TRUE) 
length(train_dl)
batch <- dataloader_make_iter(train_dl) %>% dataloader_next()
dim(batch$x)
dim(batch$y)


# Plot --------------------------------------------------------------------
par(mfrow = c(8,8), mar = rep(0, 4))
images <- as.array(batch$x$squeeze(2))

images[1,,] %>%
  #purrr::array_tree(1) %>%
  purrr::map(as.raster) %>%
  purrr::iwalk(~{plot(.x)})

images %>%
  purrr::array_tree(1) %>%
  purrr::map(as.raster) %>%
  purrr::iwalk(~{plot(.x)})

## Truth
batch$y %>% as.numeric() %>% round(digits = 2)
par(mfrow = c(1,1))


# New batches -------------------------------------------------------------------------
valid_dl <- dataloader(valid_ds, batch_size = 64)
length(valid_dl)

test_dl <- dataloader(test_ds, batch_size = 64)
length(test_dl)


# Create the model -------------------------------------------------------------------------

torch_manual_seed(777)

net <- nn_module(
  
  "corr-cnn",
  
  initialize = function() {
    
    self$conv1 <- nn_conv2d(in_channels = 1, out_channels = 32, kernel_size = 3)
    self$conv2 <- nn_conv2d(in_channels = 32, out_channels = 64, kernel_size = 3)
    self$conv3 <- nn_conv2d(in_channels = 64, out_channels = 128, kernel_size = 3)
    
    self$fc1 <- nn_linear(in_features = 14 * 14 * 128, out_features = 128)
    self$fc2 <- nn_linear(in_features = 128, out_features = 1)
    
  },
  
  forward = function(x) {
    
    x %>% 
      self$conv1() %>%
      nnf_relu() %>%
      nnf_avg_pool2d(2) %>%
      
      self$conv2() %>%
      nnf_relu() %>%
      nnf_avg_pool2d(2) %>%
      
      self$conv3() %>%
      nnf_relu() %>%
      nnf_avg_pool2d(2) %>%
      
      torch_flatten(start_dim = 2) %>%
      self$fc1() %>%
      nnf_relu() %>%
      
      self$fc2()
  }
)



model <- net()
model(batch$x)

## Train the network 
fitted <- net %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true$unsqueeze(2)),
    optimizer = optim_adam
  ) %>%
  fit(train_dl, epochs = 2, valid_data = test_dl)


preds <- predict(fitted, test_dl)


preds <- preds$to(device = "cpu")$squeeze() %>% as.numeric()
test_dl <- dataloader(test_ds, batch_size = 5000)
targets <- (test_dl %>% dataloader_make_iter() %>% dataloader_next())$y %>% as.numeric()
df <- data.frame(preds = preds, targets = targets)
library(ggplot2)
ggplot(df, aes(x = targets, y = preds)) +
  geom_point(size = 0.1) +
  theme_classic() +
  xlab("true correlations") +
  ylab("model predictions")


