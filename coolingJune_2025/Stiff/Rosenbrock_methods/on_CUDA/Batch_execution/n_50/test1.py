import numpy as np
from scipy.integrate import solve_ivp

n = 50

def derivs(t, y):
    dydt = np.zeros_like(y)
    for i in range(n):
        a = 0.1 + 0.01 * i
        b = 0.2 + 0.005 * i
        ip1 = (i + 1) % n
        ip2 = (i + 2) % n
        dydt[i] = -a * y[i] + b * y[ip1] * np.sin(y[ip2])
    return dydt

def jacobian(t, y):
    J = np.zeros((n, n))
    for i in range(n):
        a = 0.1 + 0.01 * i
        b = 0.2 + 0.005 * i
        ip1 = (i + 1) % n
        ip2 = (i + 2) % n
        J[i, i] = -a
        J[i, ip1] = b * np.sin(y[ip2])
        J[i, ip2] = b * y[ip1] * np.cos(y[ip2])
    return J

# Initial condition
y0 = np.ones(n)
t_span = (0.0, 1.0)

sol = solve_ivp(derivs, t_span, y0, method='Radau', jac=jacobian, rtol=1e-6, atol=1e-8)

print("Final state:", sol.y[:, -1])

