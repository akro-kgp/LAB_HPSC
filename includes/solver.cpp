#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm> // For std::fill
#include <iomanip>   // For std::fixed, std::setprecision

// Use the standard namespace
using namespace std;

// --- Type Definitions for clarity and compatibility ---
// The space in 'vector< vector<double> >' is kept for older compilers.
typedef vector< vector<double> > matrix;
typedef vector<double> vec;


// --- Function Prototypes ---
void read_matrix(const string& filename, matrix& A);
void read_vector(const string& filename, vec& b);
void write_matrix(const string& filename, const matrix& T);
double get_norm_diff(const vec& x_old, const vec& x_new);
int solve_jacobi(const matrix& A, const vec& b, vec& x, double tol, int max_iter);
int solve_gs(const matrix& A, const vec& b, vec& x, double tol, int max_iter);
int solve_sor(const matrix& A, const vec& b, vec& x, double omega, double tol, int max_iter);
void create_jacobi_iteration_matrix(const matrix& A, matrix& T_J);


int main() {
    int n; // Size of the system
    ifstream f_info("kinfo.txt");
    if (!f_info.is_open()) {
        cerr << "Could not open kinfo.txt" << endl;
        return 1;
    }
    f_info >> n;
    f_info.close();

    cout << "System size n = " << n << endl;

    // ---- Allocate Memory using C++ vectors ----
    matrix A(n, vec(n));
    vec b(n);
    vec x(n);

    // ---- Read Matrix and Vector from files ----
    // The previous code generated these files.
    read_matrix("Kmat.txt", A);
    read_vector("Fvec.txt", b);

    // ---- Set Solver Parameters ----
    double tolerance = 1e-6;
    int max_iterations = 20000;
    int iterations;

    cout << "\n--- Solving the System ---" << endl;

    // ---- Jacobi Method ----
    fill(x.begin(), x.end(), 0.0); // Reset initial guess to zeros
    iterations = solve_jacobi(A, b, x, tolerance, max_iterations);
    if (iterations >= max_iterations) {
         cout << "Jacobi failed to converge in " << max_iterations << " iterations." << endl;
    } else {
         cout << "Jacobi converged in " << iterations << " iterations." << endl;
    }


    // ---- Gauss-Seidel Method ----
    fill(x.begin(), x.end(), 0.0); // Reset initial guess
    iterations = solve_gs(A, b, x, tolerance, max_iterations);
    if (iterations >= max_iterations) {
         cout << "Gauss-Seidel failed to converge in " << max_iterations << " iterations." << endl;
    } else {
         cout << "Gauss-Seidel converged in " << iterations << " iterations." << endl;
    }

    // ---- SOR Method ----
    double omegas[] = {1.2, 1.5, 1.8};
    int num_omegas = sizeof(omegas)/sizeof(omegas[0]);
    for (int i = 0; i < num_omegas; i++) {
        fill(x.begin(), x.end(), 0.0); // Reset initial guess
        iterations = solve_sor(A, b, x, omegas[i], tolerance, max_iterations);
        if (iterations >= max_iterations) {
             cout << "SOR (omega=" << omegas[i] << ") failed to converge in " << max_iterations << " iterations." << endl;
        } else {
             cout << "SOR (omega=" << omegas[i] << ") converged in " << iterations << " iterations." << endl;
        }
    }

    // ---- Prepare for Spectral Radius Calculation ----
    cout << "\n--- Generating Iteration Matrix for Spectral Analysis ---" << endl;
    matrix T_J(n, vec(n));

    create_jacobi_iteration_matrix(A, T_J);
    write_matrix("T_Jacobi.txt", T_J);
    cout << "Jacobi iteration matrix written to T_Jacobi.txt" << endl;

    // Memory is automatically freed by vector destructors.

    return 0;
}

// ---- Function Implementations ----

void read_matrix(const string& filename, matrix& A) {
    ifstream file(filename.c_str());
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> A[i][j];
        }
    }
    file.close();
}

void read_vector(const string& filename, vec& b) {
    ifstream file(filename.c_str());
    int n = b.size();
    for (int i = 0; i < n; i++) {
        file >> b[i];
    }
    file.close();
}

void write_matrix(const string& filename, const matrix& T) {
    ofstream file(filename.c_str());
    int n = T.size();
    file << n << endl; // Write size for easy reading in other tools
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file << fixed << setprecision(15) << T[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

double get_norm_diff(const vec& x_old, const vec& x_new) {
    double norm = 0.0;
    for (size_t i = 0; i < x_old.size(); i++) {
        norm += (x_new[i] - x_old[i]) * (x_new[i] - x_old[i]);
    }
    return sqrt(norm);
}

int solve_jacobi(const matrix& A, const vec& b, vec& x, double tol, int max_iter) {
    int n = A.size();
    vec x_new(n);
    int k = 0;
    while (k < max_iter) {
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sigma += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) / A[i][i];
        }

        if (get_norm_diff(x, x_new) < tol) {
            x = x_new;
            return k + 1;
        }

        x = x_new; // Update x for the next iteration
        k++;
    }
    return k; // Return max_iter if not converged
}

int solve_gs(const matrix& A, const vec& b, vec& x, double tol, int max_iter) {
    int n = A.size();
    vec x_old(n);
    int k = 0;
    while (k < max_iter) {
        x_old = x; // Store previous iteration's values
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            // Use newly computed x[j] for j < i
            for (int j = 0; j < i; j++) {
                sigma += A[i][j] * x[j];
            }
            // Use old x_old[j] for j > i
            for (int j = i + 1; j < n; j++) {
                sigma += A[i][j] * x_old[j];
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }

        if (get_norm_diff(x_old, x) < tol) {
            return k + 1;
        }
        k++;
    }
    return k;
}

int solve_sor(const matrix& A, const vec& b, vec& x, double omega, double tol, int max_iter) {
    int n = A.size();
    vec x_old(n);
    int k = 0;
    while (k < max_iter) {
        x_old = x;
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if(i != j) {
                    // This uses the most recently updated values of x
                    sigma += A[i][j] * x[j];
                }
            }
            x[i] = (1 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma);
        }

        if (get_norm_diff(x_old, x) < tol) {
            return k + 1;
        }
        k++;
    }
    return k;
}

void create_jacobi_iteration_matrix(const matrix& A, matrix& T_J) {
    int n = A.size();
    // T_J = D^-1 * (L + U) = I - D^-1 * A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                T_J[i][j] = 0.0;
            } else {
                T_J[i][j] = -A[i][j] / A[i][i];
            }
        }
    }
}

