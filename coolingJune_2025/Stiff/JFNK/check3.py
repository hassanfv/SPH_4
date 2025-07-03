import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator

# Nonlinear system
def F(x):
    return np.array([
        x[0]**2 + x[1]**2 - 1,
        x[0] - x[1]**2
    ])

# Jacobian-vector product approximation (finite difference)
def Jv(x, v, eps=1e-7):
    return (F(x + eps * v) - F(x)) / eps

# JFNK solver
def jfnk(F, x0, tol=1e-8, max_iter=20):
    x = x0.copy()
    for k in range(max_iter):
        Fx = F(x)
        if np.linalg.norm(Fx) < tol:
            break
        
        # Define Jacobian-vector product as a LinearOperator
        def Jv_prod(v):
            return Jv(x, v)
        
        # Wrap Jv_prod in a LinearOperator
        J_operator = LinearOperator(
            shape=(len(x), len(x)),
            matvec=Jv_prod
        )
        
        # Solve J(x) δ = -F(x) using GMRES (without 'tol' parameter)
        delta, _ = gmres(
            J_operator,
            -Fx,
            maxiter=100  # Removed 'tol' to avoid the error
        )
        
        x += delta
        print(f"Iteration {k}: x = {x}, ||F(x)|| = {np.linalg.norm(Fx)}")
    
    return x

# Initial guess
x0 = np.array([1.0, 1.0])

# Solve
x_solution = jfnk(F, x0)
print("\nFinal solution (manual JFNK):", x_solution)
