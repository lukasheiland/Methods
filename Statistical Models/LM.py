import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

n = 100
slopes = []

for i in range(100):
    
    x = np.random.uniform(-1, 1, n)
    y_hat = x * 3
    y = np.random.normal(y_hat, 2)
    
    df = pd.DataFrame({'x': x, 'y': y, 'y_hat': y_hat})
    result = stats.linregress(x, y)
    slopes = slopes + [result.slope]
    print(result.slope)


# print(df['y_hat'])
print(np.mean(slopes))
print(np.std(slopes))
print(type(df['x']))

## Plotting
# Classes that are similar to arrays ('array-like') such as pandas data objects and numpy.matrix may not work as intended. Common convention is to convert these to numpy.array objects prior to plotting.
# print(type(y))
data = {'x': x, 'y': y}

fig, ax = plt.subplots(2) # returns a tuple with 1. Figure, 2. ONE subplot, for multiple: plt.subplots(2, 2)

c = np.ndarray((x.size,))
c.fill(1)

ax[0].scatter('x', 'y', c = c, alpha = np.random.uniform(0.4, 1, n), data=data)

plt.show()

ax[1].hist(slopes, histtype = 'step')
ax[1].axvline(3)
plt.show()

plt.clf() ## clear figure
plt.scatter(x, y)
plt.show()
