#!/usr/bin/env python3
import pandas as pd
import numpy as np
# import scipy.stats as stats

import torch
from torch import nn
from torch.utils.data import DataLoader


## Simulate a linear process
n = 1000
x = np.random.uniform(-1, 1, n)
tt = np.random.binomial(1, 0.2, n)
y_hat = x * 3
y = np.random.normal(y_hat, 2)
df = pd.DataFrame({'x': x, 'y': y, 'y_hat': y_hat, 'tt': np.array(tt, dtype = bool)})


## Set up


df[df['tt']]
df[0 == df['tt']]

# Get cpu or gpu device for training.
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
# print(f"Using {device} device")

## Define model
class NeuralNetwork(nn.Module):
	def __init__(self):
		super().__init__()
		self.flatten = nn.Flatten()
		self.linear_relu_stack = nn.Sequential(
			nn.Linear(28*28, 512),
			nn.ReLU(),
			nn.Linear(512, 512),
			nn.ReLU(),
			nn.Linear(512, 10)
		)
		
	def forward(self, x):
		x = self.flatten(x)
		logits = self.linear_relu_stack(x)
		return logits
	
model = NeuralNetwork().to(device)
print(model)

## Loss function
loss_fn = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)
