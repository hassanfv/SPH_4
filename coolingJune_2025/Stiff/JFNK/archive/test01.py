import numpy as np
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator

# Right-hand side function from CUDA code
def f(y):
    dydx = np.zeros_like(y)
    dydx[0] = -0.013 * y[0] - 1000.0 * y[0] * y[2]
    dydx[1] = -2500.0 * y[1] * y[2]
    dydx[2] = dydx[0] + dydx[1]  # since dydx[2] is sum of both
    return dydx

# Residual for backward Euler
def residual(y_new, y_old, h):
    return y_new - y_old - h * f(y_new)

# Jacobian-vector product approximation (finite difference)
def jvp(v, y, y_old, h, eps=1e-8):
    return (residual(y + eps * v, y_old, h) - residual(y, y_old, h)) / eps

# Jacobian-Free Newton-Krylov solver
def solve_jfnk_system(y0, h, steps):
    ys = [y0.copy()]
    y_current = y0.copy()
    
    for step in range(steps):
        y_old = y_current.copy()
        y_new = y_old.copy()

        for newton_iter in range(10):
            F = residual(y_new, y_old, h)

            # Define matrix-vector product for GMRES
            def mv(v):
                return jvp(v, y_new, y_old, h)

            n = len(y_new)
            J = LinearOperator((n, n), matvec=mv, dtype=np.float64)
            dx, info = gmres(J, -F, atol=1e-4, maxiter=20)

            y_new += dx

            if np.linalg.norm(dx) < 1e-10:
                break

        y_current = y_new
        ys.append(y_current.copy())

    return np.array(ys)

# Run the solver
y0 = np.array([1.0, 1.0, 0.0])
h = 0.001   # small step size due to stiffness
steps = 10000
results = solve_jfnk_system(y0, h, steps)

# Optional: plot
import matplotlib.pyplot as plt

ts = np.linspace(0, h * steps, steps + 1)
plt.plot(ts, results[:, 0], label='y0')
plt.plot(ts, results[:, 1], label='y1')
plt.plot(ts, results[:, 2], label='y2')
plt.xlabel("Time")
plt.ylabel("Solution")
plt.legend()
plt.title("JFNK Solution to Stiff ODE System")
plt.grid(True)
plt.show()

