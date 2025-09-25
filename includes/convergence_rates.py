import numpy as np

def calculate_rates(name, T):
    """Calculates and prints spectral radius and convergence rate for a given matrix T."""
    # It's possible for eigenvalues to be complex, so we take the absolute value
    eigenvalues = np.linalg.eigvals(T)
    spectral_radius = np.max(np.abs(eigenvalues))
    
    # The method converges only if the spectral radius is less than 1
    if spectral_radius < 1.0:
        convergence_rate = -np.log10(spectral_radius)
        print(f"{name}:")
        print(f"  - Spectral Radius (rho): {spectral_radius:.6f}")
        print(f"  - Convergence Rate (R):  {convergence_rate:.6f}")
    else:
        print(f"{name}:")
        print(f"  - Spectral Radius (rho): {spectral_radius:.6f} (>= 1, DOES NOT CONVERGE)")
        print(f"  - Convergence Rate (R):  N/A")
    print("-" * 30)


# --- Main Script ---

# Read the size and the matrix A from fem.c op
with open("kinfo.txt", "r") as f:#get the dimensions
    n = int(f.read())

A = np.loadtxt("Kmat.txt").reshape((n, n))#coeff matrix

# Decompose the matrix A into D (diagonal), L (strict lower), and U (strict upper)
D = np.diag(np.diag(A))
L = -np.tril(A, k=-1) # tril gets the lower triangle, - sign because A = D - L - U
U = -np.triu(A, k=1) # triu gets the upper triangle


D_inv = np.linalg.inv(D)

#Iter matrices
# 1. Jacobi Method
T_J = D_inv @ (L + U)
calculate_rates("Jacobi", T_J)

# 2. Gauss-Seidel Method
T_GS = np.linalg.inv(D - L) @ U
calculate_rates("Gauss-Seidel", T_GS)

# 3. SOR Method for different omega values
omegas = [1.20, 1.50, 1.80]
for w in omegas:
    # T_SOR = (D - wL)^-1 * ((1-w)D + wU)
    M = np.linalg.inv(D - w * L)
    N = (1 - w) * D + w * U
    T_SOR = M @ N
    calculate_rates(f"SOR (omega={w:.2f})", T_SOR)