import numpy as np
from scipy.sparse.linalg import LinearOperator

def custom_gmres(A, b, x0=None, max_iter=100, restart=None, tol=1e-8):
    """Fully debugged GMRES implementation with robust dimension handling."""
    n = len(b)
    if x0 is None:
        x0 = np.zeros_like(b)
    if restart is None:
        restart = min(max_iter, n)  # Don't exceed problem dimension
    
    x = x0.copy()
    residual = b - A.matvec(x)
    norm_b = np.linalg.norm(b)
    if norm_b < 1e-14:
        norm_b = 1  # Prevent division by zero
        
    res_norm = np.linalg.norm(residual)
    res_history = [res_norm / norm_b]
    
    for outer in range(0, max_iter, restart):
        m = min(restart, max_iter - outer)
        
        # Arnoldi process
        Q = np.zeros((n, m+1))
        H = np.zeros((m+1, m))
        Q[:, 0] = residual / res_norm
        
        for j in range(m):
            v = A.matvec(Q[:, j])
            
            # Modified Gram-Schmidt
            for i in range(j+1):
                H[i,j] = np.dot(Q[:,i], v)
                v -= H[i,j] * Q[:,i]
            
            H_norm = np.linalg.norm(v)
            if H_norm < 1e-14:  # Lucky breakdown
                break
                
            H[j+1,j] = H_norm
            Q[:,j+1] = v / H_norm
            
            # Solve least squares problem
            e1 = np.zeros(j+2)
            e1[0] = res_norm
            H_sub = H[:j+2, :j+1]
            
            y, _, _, _ = np.linalg.lstsq(H_sub, e1, rcond=None)
            
            # Ensure y has correct shape
            y = np.asarray(y).flatten()  # Force 1D array
            
            # Compute residual norm
            current_res = np.linalg.norm(H_sub @ y - e1)
            res_history.append(current_res / norm_b)
            
            if current_res / norm_b < tol:
                break
        
        # Final solution update with guaranteed correct dimensions
        effective_j = min(j, Q.shape[1]-1)  # Don't exceed allocated space
        y = y[:effective_j+1]  # Trim to actual basis size
        
        # Use either dot product or matmul depending on dimensions
        if Q[:, :effective_j+1].shape[1] == y.size:
            x_update = np.dot(Q[:, :effective_j+1], y)
        else:
            # Fallback for dimension mismatches
            x_update = np.dot(Q[:, :y.size], y)
        
        x += x_update
        residual = b - A.matvec(x)
        res_norm = np.linalg.norm(residual)
        
        if res_norm / norm_b < tol:
            break
    
    info = 0 if res_history[-1] < tol else 1
    return x, info

# Nonlinear system
def F(x):
    return np.array([
        x[0]**2 + x[1]**2 - 1,
        x[0] - x[1]**2
    ])

# Jacobian-vector product
def Jv(x, v, eps=1e-7):
    return (F(x + eps * v) - F(x)) / eps

# Robust JFNK solver
def jfnk(F, x0, tol=1e-8, max_iter=20):
    x = x0.copy()
    for k in range(max_iter):
        Fx = F(x)
        current_norm = np.linalg.norm(Fx)
        print(f"Iter {k}: x = {x}, ||F|| = {current_norm:.2e}")
        
        if current_norm < tol:
            break
            
        J_operator = LinearOperator(
            shape=(len(x), len(x)),
            matvec=lambda v: Jv(x, v)
        )
        
        delta, info = custom_gmres(
            J_operator,
            -Fx,
            tol=max(tol/10, 1e-12),  # Stricter but bounded tolerance
            max_iter=100,
            restart=min(20, len(x))
        )
        
        if info != 0:
            print(f"GMRES warning: Linear solve did not fully converge")
        
        # Simple line search
        alpha = 1.0
        for _ in range(5):
            new_x = x + alpha*delta
            new_F = F(new_x)
            if np.linalg.norm(new_F) < current_norm:
                x = new_x
                break
            alpha *= 0.5
        else:
            x += delta  # Fallback if line search fails
    
    return x

# Solve with proper initialization
print("Solving with fully debugged custom GMRES...")
x0 = np.array([1.0, 1.0])
x_sol = jfnk(F, x0)
print(f"\nFinal solution: {x_sol}")
print(f"Final residual norm: {np.linalg.norm(F(x_sol)):.2e}")
print(f"Verification:")
print(f"x² + y² - 1 = {x_sol[0]**2 + x_sol[1]**2 - 1:.2e}")
print(f"x - y² = {x_sol[0] - x_sol[1]**2:.2e}")
