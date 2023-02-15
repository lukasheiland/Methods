#!/usr/bin/env python3
import pandas as pd
import numpy as np
# import scipy.stats as stats
# from sklearn.preprocessing import StandardScaler as scaler

import torch
from torch import nn

import matplotlib.pyplot as plt


## Simulate a linear process
n = 3000
X = np.random.uniform(-1, 1, [n, 2])
tt = np.random.binomial(1, 0.2, n)
## y_hat = X * [3, 2] ## this would be elementwise
## y_hat = np.dot(X, [1, 0]) ## matmul differs from dot in two important ways: Multiplication by scalars is not allowed, Stacks of matrices are broadcast together as if the matrices were elements
y_hat = np.matmul(X, [3, 2])
y = np.random.normal(y_hat, 1)
df = pd.DataFrame({'y': y, 'y_hat': y_hat, 'tt': np.array(tt, dtype = bool)})
df = pd.concat([df, pd.DataFrame(X)], axis = 1)
n_x = X.shape[1]

## Set up
# All hidden layers usually use the same activation function.
# However, the output layer will typically use a different activation function from the hidden layers. The choice depends on the goal or type of prediction made by the model.

# stochastic gradient decent: minimizing the total loss of the network
# Forward Propagation: info (firings) move from features through intermediate activation functions to the output layer
# Back propagation: loss is minimized to inform weights [and parameters of the activation functions] (feedback with derivatives of the loss!)
# the purpose of an activation function is to add non-linearity (otherwise it would just be weights, i.e. (compositions of, which are) linear functions)
# non-linear activation functions allow backpropagation because with non-linear functions the derivatives change related to the input -> feedback (linears would be constant)
# activation functions include:
#   logistic: output of the logistic function is not symmetric around zero. So the output of all the neurons will be of the same sign -> unstable training
#   tanh (hyperbolic tangent): output of the tanh activation function is Zero centered; hence we can easily map the output values as strongly negative, neutral, or strongly positive; but vanishing gradients at the extremes!
#   ReLU (rectified linear unit): differentiable (backpropagating!) version of only activating when positive (then linear); accelerates the convergence of gradient descent towards the global minimum of the loss function due to its linear, non-saturating property.
#   Leaky ReLU: instead of flat, small slope in the negative range, avoiding the problems that some neurons might never be activated (backpropagation for negatives!)
#   parameteric relu, where the slope of the negative part is learned
#   ELU: ReLu with a log curve in the negative range, (exploding gradient [of the loss function] problem: unstable minimization due to error gradients accumulating during an update and result in very large gradients. These in turn result in large updates to the network weights, and in turn, an unstable network)
#   swish: better ReLu, smooth continuous function
#   sofmax: sigmoid outputting [0, 1]; mostly used for last layer in multiclass classification


df_train = df[df['tt']]
df_test = df[~df['tt']]


## Get cpu or gpu device for training.
# device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"


## Define model by subclassing nn.Module
class NeuralNetwork(nn.Module):
	def __init__(self, input_dim, hidden_dim, output_dim):
		
		super().__init__()
		
		# self.flatten = nn.Flatten() ## use the Flatten class
		
		## define layers, three layer network (input layer not counted)
		self.linear_relu_stack = nn.Sequential(
			nn.Linear(input_dim, hidden_dim),
			nn.ReLU(),
			nn.Linear(hidden_dim, hidden_dim),
			nn.ReLU(),
			nn.Linear(hidden_dim, output_dim)
		)
		
	def forward(self, x): ## specify how the input x is forwarded through the layers.
		# x = self.flatten(x)
		x = self.linear_relu_stack(x) ## Sequential just expects x
		return x

input_dim = n_x ## Dimensions for the network: Number of features!
hidden_dim = 10
output_dim = 1

net = NeuralNetwork(input_dim, hidden_dim, output_dim).to("cpu")
print(net)

def weights_init(m):
    if isinstance(m, nn.Linear):
       nn.init.xavier_uniform_(m.weight)

net.apply(weights_init) ## Applies ``fn`` recursively to every submodule

## Loss functions and purposes
# Regression: MSE(Mean Squared Error), MAE(Mean Absolute Error) Hubber loss.
# Classification: Binary cross-entropy (log-loss), Categorical cross-entropy.
# AutoEncoder: KL Divergence.
# Object detection: Focal loss.
# Word embeddings

loss_fn = nn.MSELoss() ## suitable for regression task


learning_rate = 0.1 ## "step size" for optimizer
optimizer = torch.optim.SGD(net.parameters(), lr = learning_rate)
num_epochs = 50
loss_values = []

# creating tensor from targets_df 
X_train = torch.tensor(X[df["tt"]]).float()
# X = scaler.fit_transform(X)
y = torch.tensor(df_train[["y"]].values).float()

## Inner loop
def run_epochs(num_epochs):
  for epoch in range(num_epochs):
    ## zero the parameter gradients (gradients are accumulated during the backward pass, so if you don't zero them at the beginning of each epoch, they will be accumulated across epochs)
    optimizer.zero_grad()
    
    ## forward pass
    predictions = net(X_train)
    loss = loss_fn(predictions, y)
    loss_values.append(loss.item())
    
    ## backward pass & optimization
    loss.backward()
    optimizer.step()

## RUN
run_epochs(num_epochs)

## Plot loss
plt.clf()
plt.plot(range(num_epochs), loss_values, color='red')
plt.ylabel("Loss")
plt.show()

## Plot pred vs. true
X_test = torch.tensor(X[~df["tt"]]).float()
y_test = df.loc[~df["tt"], "y"]
predictions_train = net(X_train).detach().numpy()
predictions_test = net(X_test).detach().numpy()

plt.clf()
plt.scatter(y_test, predictions_test, color='blue')
plt.scatter(y, predictions_train, color='red')
plt.xlabel("True Values")
plt.show()


# print(list(net.parameters()))

## test MAE
mae_test = abs(predictions_test.squeeze() - y_test.to_numpy()).mean()






