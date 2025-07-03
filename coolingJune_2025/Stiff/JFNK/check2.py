

import numpy as np
from scipy.optimize import newton_krylov

# Define the nonlinear system F(x) = 0
def F(x):
    return np.array([
        x[0]**2 + x[1]**2 - 1,
        x[0] - x[1]**2
    ])

# Initial guess
x0 = np.array([1.0, 1.0])

# Solve using JFNK (GMRES is the default Krylov method)
x_solution = newton_krylov(F, x0)

print("Solution (using scipy.optimize.newton_krylov):", x_solution)
