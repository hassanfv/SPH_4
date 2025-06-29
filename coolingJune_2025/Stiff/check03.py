# Forward Euler !

import numpy as np
import matplotlib.pyplot as plt

def f1(x, t):
    return -10. * x


#print('Result = ', f1(1e-3))

h = 0.15

t = 0.0
x = 1.0 # initial value of u at t = 0 is 1.0 ! 

res = [[t, x]]

while t < 2:

    t += h
    x = x + h * f1(x, t)

    res.append([t, x])

res = np.array(res)

xgrid = np.linspace(0, 1, 100)
yy = np.exp(-10.*xgrid)
print(res[:, 1])

plt.plot(xgrid, yy, color = 'k', linewidth = 2)
plt.scatter(res[:, 0], res[:, 1], color = 'r', s = 50)
plt.show()
