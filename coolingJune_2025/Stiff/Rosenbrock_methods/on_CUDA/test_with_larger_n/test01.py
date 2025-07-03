import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the ODE system
def derivs(t, y):
    dydt = np.zeros_like(y)
    dydt[0] = -10.0 * y[0] + y[1] * y[2]
    dydt[1] = -0.1 * y[1] + y[0] * y[2]
    dydt[2] = -100.0 * y[2] + y[3]
    dydt[3] = -0.5 * y[3] + y[4] * y[0]
    dydt[4] = -y[4] + y[5] * y[6]
    dydt[5] = -0.01 * y[5] + y[6]
    dydt[6] = -20.0 * y[6] + y[7]
    dydt[7] = -0.2 * y[7] + y[8] * y[1]
    dydt[8] = -y[8] + y[9] * y[0]
    dydt[9] = -5.0 * y[9] + y[0] * y[2]
    return dydt

# (Optional) Provide the Jacobian for better performance with stiff solvers
def jacobian(t, y):
    n = 10
    J = np.zeros((n, n))

    J[0, 0] = -10.0
    J[0, 1] = y[2]
    J[0, 2] = y[1]

    J[1, 0] = y[2]
    J[1, 1] = -0.1
    J[1, 2] = y[0]

    J[2, 2] = -100.0
    J[2, 3] = 1.0

    J[3, 0] = y[4]
    J[3, 3] = -0.5
    J[3, 4] = y[0]

    J[4, 4] = -1.0
    J[4, 5] = y[6]
    J[4, 6] = y[5]

    J[5, 5] = -0.01
    J[5, 6] = 1.0

    J[6, 6] = -20.0
    J[6, 7] = 1.0

    J[7, 1] = y[8]
    J[7, 7] = -0.2
    J[7, 8] = y[1]

    J[8, 0] = y[9]
    J[8, 8] = -1.0
    J[8, 9] = y[0]

    J[9, 0] = y[2]
    J[9, 2] = y[0]
    J[9, 9] = -5.0

    return J

# Initial condition
y0 = np.array([1.0]*10)
t_span = (0.0, 5.0)
t_eval = np.linspace(t_span[0], t_span[1], 500)

# Solve using solve_ivp with BDF (good for stiff systems)
sol = solve_ivp(derivs, t_span, y0, method='BDF', jac=jacobian, t_eval=t_eval, rtol=1e-6, atol=1e-9)

# Print final values
print("Final state y(t=50):")
for i, yi in enumerate(sol.y[:, -1]):
    print(f"y[{i}] = {yi:.6f}")

# Plot selected components
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 6))
for i in range(10):
    plt.plot(sol.t, sol.y[i], label=f'y[{i}]')
plt.xlabel('Time')
plt.ylabel('y values')
plt.title('ODE system solution (n=10)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


