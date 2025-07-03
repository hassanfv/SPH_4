import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the Robertson ODE system
def robertson(t, y, k1, k2, k3):
    y1, y2, y3 = y
    dy1dt = -k1 * y1 + k3 * y2 * y3
    dy2dt = k1 * y1 - k2 * y2**2 - k3 * y2 * y3
    dy3dt = k2 * y2**2
    return [dy1dt, dy2dt, dy3dt]

# Initial conditions and parameters
y0 = [1.0, 0.0, 0.0]
t_span = (0, 1e5)
t_eval = np.logspace(-2, 5, 1000)  # log-spaced time points for plotting
params = (0.04, 3e7, 1e4)  # k1, k2, k3

# Solve the ODE
sol = solve_ivp(robertson, t_span, y0, args=params, method='BDF', t_eval=t_eval)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='y1')
plt.plot(sol.t, sol.y[1], label='y2')
plt.plot(sol.t, sol.y[2], label='y3')
plt.xscale('log')
plt.xlabel('Time')
plt.ylabel('Values')
plt.title('Robertson Problem')
plt.legend()
plt.grid(True)
plt.show()

