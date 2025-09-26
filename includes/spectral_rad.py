import numpy as np

# Load the iteration matrix from the file created by your C program
T = np.loadtxt(r'D:\4th year\7th sem\Hpss\LAB_HPSC\matrix.txt')

print("Matrix T:\n", T)
#with open(r'C:\Users\Chiradip Biswas\Desktop\kgp\7thSem\high performance computing\lab_class2\matrix.txt', 'r') as f:
#    n = int(f.read())
#    T = np.loadtxt(f)

# Calculate eigenvalues
eigenvalues = np.linalg.eigvals(T)

# Spectral radius is the maximum absolute eigenvalue
spectral_radius = np.max(np.abs(eigenvalues))

# Asymptotic rate of convergence
convergence_rate = -np.log10(spectral_radius)

# Perform Singular Value Decomposition of matrix T
U, S, VT = np.linalg.svd(T, full_matrices=True)
#1st r diagonal elements of s will give eig values of AAT=> lambda_A^2. for the linear ly indep r cols of A

print("eigen values: ",eigenvalues)
print("U:\n", U)
print("Singular values:\n", S)
print("V^T:\n", VT)

print(f"Spectral Radius (rho): {spectral_radius}")
print(f"Convergence Rate (R): {convergence_rate}")