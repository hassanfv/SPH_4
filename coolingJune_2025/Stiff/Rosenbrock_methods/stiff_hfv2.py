import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
import time


import numpy as np
from scipy.integrate import solve_ivp

def f(t, y):
    y0, y1, y2 = y
    return [
        -0.013*y0 - 1000*y0*y2,
        -2500*y1*y2,
        -0.013*y0 - 1000*y0*y2 - 2500*y1*y2
    ]

y0 = [1.0, 1.0, 0.0]
TA = time.time()
sol = solve_ivp(f, (0, 100), y0, method='Radau', rtol=1e-4, atol=1e-4,
                dense_output=True)          # Radau = implicit Runge–Kutta
print('Elapsed time in solve_ivp = ', time.time() - TA)

t = np.linspace(0, 100, 600)
y = sol.sol(t)          # y[0]=y0, y[1]=y1, y[2]=y2


# Plot the results
plt.plot(y[1], y[2], label = 'python solve_ivp')


df = pd.read_csv("outX.csv") # From C++ stiff_hfv.cpp code !

x = df["x"]
y = df["y"]

#plt.figure(figsize = (16, 12))

plt.plot(x, y, c = 'k', label = 'C++ stiff_hfv code')
plt.scatter(x, y, c = 'r', s = 30, label = 'C++ stiff_hfv code')


plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig('python_ca_cpp.png', dpi = 200, bbox_inches = 'tight')

plt.show()




