import numpy as np
import pandas as pd
import scipy.stats as stats

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