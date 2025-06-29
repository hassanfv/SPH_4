import numpy as np
import matplotlib.pyplot as plt

def u(x):
    return 2.0 * np.exp(-x) - np.exp(-1000. * x)

def v(x):
    return -np.exp(-x) + np.exp(-1000. * x)


print("Check = ", u(0.1))


x = np.arange(0, 10, 0.0001)

#plt.plot(x, u(x), marker = 'o')
plt.plot(x, v(x), marker = 's')

plt.xlim(0.0, 0.006)


plt.show()