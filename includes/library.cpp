// #include <iostream>
// #include <vector>
// #include <stdexcept> // for std::runtime_error
// #include <cmath>     // for sqrt, abs, log10
// #include <iomanip>   // for std::fixed, std::setprecision
// #include <utility>   // for std::pair

// #include <fstream>
// #include <string>
// #include <cmath>
// #include <algorithm> // For std::fill
#include<bits/stdc++.h>

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
double vec_norm(const vec& v);
vec solve_least_squares(const matrix& H, const vec& beta_e1);
void arnoldi_iteration(const matrix& A, int m, const vec& v1, matrix& V, matrix& H, int& actual_m);


using namespace std;

// Define alias
typedef vector< vector<double> > matrix; //this one is supported, not the using matrix=...

// This function takes the compressed diagonal storage and returns the full matrix.
matrix reconstruct_from_diag(const matrix& diag, const vec& ioff, int n) {
    matrix A(n, vec(n, 0.0)); // Initialize n x n matrix with zeros

    if (diag.empty() || ioff.empty() || diag[0].size() != ioff.size()) {
        throw runtime_error("DIAG and IOFF dimensions are incompatible.");
    }
    if (diag.size() != n) {
         throw runtime_error("Number of rows in DIAG must match the matrix size n.");
    }

    int num_diagonals = ioff.size();

    // Loop through each column of DIAG, which represents a diagonal in A
    for (int j = 0; j < num_diagonals; ++j) {
        int offset = ioff[j];
        // Loop through each element in that column
        for (int i = 0; i < n; ++i) {
            int target_row = i;
            int target_col = i + offset;

            // Place the element if it's within the matrix bounds
            if (target_row >= 0 && target_row < n && target_col >= 0 && target_col < n) {
                A[target_row][target_col] = diag[i][j];
            }
        }
    }
    return A;
}

// --- NEW: Helper function to get DIAG/IOFF input from the user ---
// This function prompts the user for the sparse matrix data.
matrix get_A_from_diag_input() {
    int n;
    cout << "Enter the size (n x n) of the full matrix A: ";
    cin >> n;

    int num_diags;
    cout << "Enter the number of diagonals (size of IOFF vector): ";
    cin >> num_diags;

    vec ioff(num_diags);
    cout << "Enter the elements of the IOFF vector (e.g., -1 0 1 2): ";
    for(int i = 0; i < num_diags; ++i) {
        cin >> ioff[i];
    }

    matrix diag(n, vec(num_diags));
    cout << "Enter the " << n << "x" << num_diags << " DIAG matrix row by row." << endl;
    cout << "NOTE: For positions with '*', you can enter 0 or any number, it will be ignored." << endl;
    for (int i = 0; i < n; ++i) {
        cout << "Row " << i + 1 << ": ";
        for (int j = 0; j < num_diags; ++j) {
            cin >> diag[i][j];
        }
    }

    cout << "\nReconstructing full matrix A..." << endl;
    matrix A = reconstruct_from_diag(diag, ioff, n);
    return A;
}
matrix make_mat(int r1, int c1) {
    matrix a(r1, vector<double>(c1, 0));
    cout << "Enter the elements:\n";
    for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c1; j++) {
            cin >> a[i][j];
        }
    }
    return a;
}
vec make_vec(int r){
    vec a;
    cout<<"---- Input b elements ----\n";
    for(int i=0;i<r;i++)
    {
        int g;
        cin>>g;
        a.push_back(g);
    }
    return a;
}
void print_mat(const matrix& a) {   // pass by const reference (not copy)
    for (int i =0 ; i<a.size();i++) {
        for (int j =0 ; j<a[0].size();j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

matrix mat_mul(const matrix& A, const matrix& B) {
    if (A.empty() || B.empty()) {
        throw runtime_error("Empty matrix provided.");
    }
    if (A[0].size() != B.size()) {
        throw runtime_error("Matrix dimensions do not allow multiplication.");
    }

    size_t m = A.size();
    size_t n = A[0].size();
    size_t p = B[0].size();

    matrix C(m, vector<double>(p, 0));

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < p; j++) {
            for (size_t k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}



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
             for(int i=0;i<x_new.size();i++) cout<<x_new[i]<<" ";
    cout<<endl;
            return k + 1;
        }

        x = x_new; // Update x for the next iteration
        k++;
    }
    for(int i=0;i<x_new.size();i++) cout<<x_new[i]<<" ";
    cout<<endl;
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
             for(int i=0;i<x.size();i++) cout<<x[i]<<" ";
    cout<<endl;
            return k + 1;
        }
        k++;
    }
    for(int i=0;i<x.size();i++) cout<<x[i]<<" ";
    cout<<endl;
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
            for(int i=0;i<x.size();i++) cout<<x[i]<<" ";
    cout<<endl;
            return k + 1;
        }
        k++;
    }
    for(int i=0;i<x.size();i++) cout<<x[i]<<" ";
    cout<<endl;
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

matrix transpose(const matrix& A) {
    if (A.empty()) return matrix();
    size_t m = A.size();
    size_t n = A[0].size();
    matrix T(n, vec(m));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            T[j][i] = A[i][j];
        }
    }
    return T;
}

void text_from_mat(const matrix A){
    ofstream fout("matrix.txt");
    for (int i=0;i<A.size();i++) {
        for (int j=0;j<A[0].size();j++) {
            fout << A[i][j] << " ";
        }
        fout << "\n";
    }
    fout.close();
}
void vec_print(const vec b){
    for(int i=0;i<b.size();i++){
        cout<<b[i]<<' ';
    }
    cout<<endl;
}

double dot_product(const vec& a, const vec& b) {
    if (a.size() != b.size()) {
        throw runtime_error("Vector sizes must match for dot product.");
    }
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

vec vec_add(const vec& a, const vec& b) {
    if (a.size() != b.size()) {
        throw runtime_error("Vector sizes must match for addition.");
    }
    vec result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

vec vec_sub(const vec& a, const vec& b) {
    if (a.size() != b.size()) {
        throw runtime_error("Vector sizes must match for subtraction.");
    }
    vec result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

vec scalar_vec_mul(double s, const vec& v) {
    vec result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = s * v[i];
    }
    return result;
}



vec mat_vec_mul(const matrix& A, const vec& x) {
    if (A.empty() || x.empty() || A[0].size() != x.size()) {
        throw runtime_error("Matrix and vector dimensions are not compatible for multiplication.");
    }
    size_t m = A.size();
    vec result(m, 0.0);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < x.size(); ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

int solve_steepest_descent(const matrix& A, const vec& b, vec& x, double tol, int max_iter) {
    if (A.empty()) return 0;
    
    vec r = vec_sub(b, mat_vec_mul(A, x)); // r = b - A*x
    
    for (int k = 0; k < max_iter; ++k) {
        
        double r_dot_r = dot_product(r, r);

        
        vec Ar = mat_vec_mul(A, r);
        double r_dot_Ar = dot_product(r, Ar);

        double alpha = r_dot_r / r_dot_Ar;
        
        x = vec_add(x, scalar_vec_mul(alpha, r));
        r = vec_sub(r, scalar_vec_mul(alpha, Ar)); 

        if (r_dot_Ar == 0 || sqrt(r_dot_r) < tol) {
            cout<<"solution:-----\n";
            vec_print(x);
            return k; 
        }

        

    }
     cout<<"solution:-----\n";
            vec_print(x);
    return max_iter;
}
int solve_minimum_residual(const matrix& A, const vec& b, vec& x, double tol, int max_iter) {
    if (A.empty()) return 0;

    vec r = vec_sub(b, mat_vec_mul(A, x)); // r = b - A*x

    for (int k = 0; k < max_iter; ++k) {
        double r_norm = sqrt(dot_product(r, r));

        if (r_norm < tol) {
             cout<<"solution:-----\n";
            vec_print(x);
            return k + 1; // Converged
        }

        vec Ar = mat_vec_mul(A, r);
        double Ar_dot_r = dot_product(Ar, r);
        double Ar_dot_Ar = dot_product(Ar, Ar);

        if (Ar_dot_Ar == 0) {
             cout<<"solution:-----\n";
            vec_print(x);
            return k + 1; // Stagnation or exact solution
        }

        double alpha = Ar_dot_r / Ar_dot_Ar;

        x = vec_add(x, scalar_vec_mul(alpha, r));
        // Recompute residual for better numerical stability
        r = vec_sub(b, mat_vec_mul(A, x)); 
    }
     cout<<"solution:-----\n";
            vec_print(x);
    return max_iter; // Failed to converge
}
int solve_residual_norm_steepest_descent(const matrix& A, const vec& b, vec& x, double tol, int max_iter) {
    if (A.empty()) return 0;

    matrix At = transpose(A);

    for (int k = 0; k < max_iter; ++k) {
        vec r = vec_sub(b, mat_vec_mul(A, x));

        if (sqrt(dot_product(r, r)) < tol) {
            cout<<"solution:-----\n";
            vec_print(x);
            return k + 1; // Converged
        }

        vec V = mat_vec_mul(At, r);
        vec AV = mat_vec_mul(A, V);

        double V_norm_sq = dot_product(V, V);
        double AV_norm_sq = dot_product(AV, AV);

        if (AV_norm_sq == 0) {
            return k + 1; // Stagnation or exact solution
        }

        double alpha = V_norm_sq / AV_norm_sq;
        
        x = vec_add(x, scalar_vec_mul(alpha, V));
    }
    cout<<"solution:-----\n";
            vec_print(x);
    return max_iter; // Failed to converge
}

// --- Arnoldi Iteration ---
// Generates an orthonormal basis V for the Krylov subspace K_m(A, v1)
// and an upper Hessenberg matrix H such that A*V_m = V_{m+1}*H_m.
void arnoldi_iteration(const matrix& A, int m, const vec& v1, matrix& V, matrix& H, int& actual_m) {
    int n = A.size();
    V.assign(n, vec(m + 1, 0.0));
    H.assign(m + 1, vec(m, 0.0));
    actual_m = m;

    // Set first vector of the basis
    for(int i=0; i<n; ++i) V[i][0] = v1[i];

    for (int j = 0; j < m; ++j) {
        vec current_v(n);
        for(int i=0; i<n; ++i) current_v[i] = V[i][j];

        vec w = mat_vec_mul(A, current_v);

        // Gram-Schmidt process
        for (int i = 0; i <= j; ++i) {
            vec v_i(n);
            for(int k=0; k<n; ++k) v_i[k] = V[k][i];

            H[i][j] = dot_product(w, v_i);
            w = vec_sub(w, scalar_vec_mul(H[i][j], v_i));
        }

        H[j + 1][j] = vec_norm(w);

        // Check for breakdown (lucky breakdown)
        if (H[j + 1][j] < 1e-12) {
            actual_m = j + 1;
            break;
        }
        
        // Normalize to get the next vector in the basis
        vec v_next = scalar_vec_mul(1.0 / H[j + 1][j], w);
        for(int i=0; i<n; ++i) V[i][j+1] = v_next[i];
    }
}


// --- Helper Function Implementations ---

double vec_norm(const vec& v) {
    return sqrt(dot_product(v, v));
}

vec solve_least_squares(const matrix& H, const vec& beta_e1) {
    int n = H.size();
    matrix LU = H;
    vec y = beta_e1;
    vec x(n, 0.0);

    // LU decomposition (Doolittle's method)
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += LU[i][k] * LU[k][j];
            }
            LU[i][j] -= sum;
        }
        for (int j = i + 1; j < n; j++) {
             if (LU[i][i] == 0) throw std::runtime_error("Matrix is singular in LS solve.");
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += LU[j][k] * LU[k][i];
            }
            LU[j][i] = (LU[j][i] - sum) / LU[i][i];
        }
    }
    
    // Forward substitution (Ly = b)
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += LU[i][j] * y[j];
        }
        y[i] = (beta_e1[i] - sum);
    }

    // Backward substitution (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += LU[i][j] * x[j];
        }
        if (LU[i][i] == 0) throw std::runtime_error("Matrix is singular in LS solve.");
        x[i] = (y[i] - sum) / LU[i][i];
    }
    
    return x;
}


// --- Full Orthogonalization Method (FOM) using Arnoldi Iteration ---
// Solves Ax = b for a general square matrix A.
int solve_fom(const matrix& A, const vec& b, vec& x, int m_restart, double tol, int max_iter) {
    int n = A.size();
    
    for (int iter = 0; iter < max_iter; ++iter) {
        vec r0 = vec_sub(b, mat_vec_mul(A, x));
        double beta = vec_norm(r0);

        if (beta < tol) {
            cout<<"solution:-----\n";
            vec_print(x);
            return iter + 1; // Converged
        }

        vec v1 = scalar_vec_mul(1.0 / beta, r0);
        
        matrix V;
        matrix H;
        int actual_m = 0; // Will be updated if breakdown occurs

        // 1. Build the orthonormal basis V and Hessenberg matrix H using Arnoldi
        arnoldi_iteration(A, m_restart, v1, V, H, actual_m);
        
        // 2. Solve the smaller system Hy = beta*e1
        matrix H_m(actual_m, vec(actual_m, 0.0));
        for(int i=0; i<actual_m; ++i) {
            for(int k=0; k<actual_m; ++k) {
                H_m[i][k] = H[i][k];
            }
        }
        
        vec beta_e1(actual_m, 0.0);
        beta_e1[0] = beta;

        vec y = solve_least_squares(H_m, beta_e1);

        // 3. Update solution x = x + V_m * y
        for (int j = 0; j < actual_m; ++j) {
            vec v_j_update(n);
             for(int i=0; i<n; ++i) v_j_update[i] = V[i][j];
            x = vec_add(x, scalar_vec_mul(y[j], v_j_update));
        }

        // Check for convergence again after update
        r0 = vec_sub(b, mat_vec_mul(A, x));
        if (vec_norm(r0) < tol) {
            cout<<"solution:-----\n";
            vec_print(x);
            return iter + 1;
        }
    }
    cout<<"solution:-----\n";
            vec_print(x);
    return max_iter;
}


int solve_conjugate_gradient(const matrix& A, const vec& b, vec& x, double tol, int max_iter) {
    vec r = vec_sub(b, mat_vec_mul(A, x));
    vec p = r;
    double rs_old = dot_product(r, r);

    if (sqrt(rs_old) < tol) {
        return 0; // Already converged
    }

    for (int i = 0; i < max_iter; ++i) {
        vec Ap = mat_vec_mul(A, p);
        double alpha = rs_old / dot_product(p, Ap);

        x = vec_add(x, scalar_vec_mul(alpha, p));
        r = vec_sub(r, scalar_vec_mul(alpha, Ap));

        double rs_new = dot_product(r, r);
        if (sqrt(rs_new) < tol) {

            cout<<"solution:-----\n";
            vec_print(x);
            return i + 1;
        }

        p = vec_add(r, scalar_vec_mul(rs_new / rs_old, p));
        rs_old = rs_new;
    }
    cout<<"solution:-----\n";
    vec_print(x);    
    return max_iter;
}


//TDMA
// Thomas algorithm (TDMA) solver for a tridiagonal matrix A.
// Assumes A is n x n and tridiagonal (only A[i][i-1], A[i][i], A[i][i+1] are used).
// Input:  A (full matrix but only tridiagonal entries are used), b (RHS)
// Output: x (solution vector, resized to n)
// Return: 0 on success, throws runtime_error on invalid input or zero pivot.
int solve_tdma(const matrix& A, const vec& b, vec& x) {
    int n = A.size();
    if (n == 0) throw runtime_error("Empty matrix A in TDMA.");
    if ((int)A[0].size() != n) throw runtime_error("Matrix A must be square.");
    if ((int)b.size() != n) throw runtime_error("Vector b size must match A.");

    // Handle trivial case n == 1
    if (n == 1) {
        if (fabs(A[0][0]) < 1e-18) throw runtime_error("Zero pivot detected in TDMA (n==1).");
        x.assign(1, b[0] / A[0][0]);
        cout << "solution:-----\n"; vec_print(x);
        return 0;
    }

    x.assign(n, 0.0);
    vec cprime(n, 0.0); // modified super-diagonal coefficients
    vec dprime(n, 0.0); // modified RHS

    // First row
    double a00 = A[0][0];
    if (fabs(a00) < 1e-18) throw runtime_error("Zero pivot detected in TDMA at row 0.");
    cprime[0] = A[0][1] / a00;           // c_0' = c0 / b0
    dprime[0] = b[0] / a00;              // d_0' = d0 / b0

    // Forward elimination
    for (int i = 1; i < n; ++i) {
        double ai_im1 = A[i][i-1];      // a_i (sub-diagonal)
        double ai_i   = A[i][i];        // b_i (diagonal)
        double denom  = ai_i - ai_im1 * cprime[i-1];
        if (fabs(denom) < 1e-18) {
            throw runtime_error("Zero pivot detected in TDMA during forward elimination at row " + to_string(i));
        }
        // Only set cprime for rows that have a super-diagonal (i < n-1)
        if (i < n - 1) {
            cprime[i] = A[i][i+1] / denom;  // c_i' = c_i / (b_i - a_i * c_{i-1}')
        }
        dprime[i] = (b[i] - ai_im1 * dprime[i-1]) / denom; // d_i' = (d_i - a_i * d_{i-1}') / (b_i - a_i * c_{i-1}')
    }

    // Back substitution
    x[n-1] = dprime[n-1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dprime[i] - cprime[i] * x[i+1];
    }

    cout << "solution:-----\n";
    vec_print(x);
    return 0;
}
// Restarted GMRES (uses arnoldi_iteration already in your library)
// A : matrix (n x n)
// b : RHS (n)
// x : initial guess (n) and on exit contains the approximate solution
// restart : m (dimension of Krylov subspace before restart)
// tol : tolerance for residual norm
// max_iter : maximum number of outer iterations * restart (total mat-vec ops approx)
// Returns: number of Arnoldi steps performed (or max_iter if not converged)
int solve_gmres(const matrix& A, const vec& b, vec& x, int restart, double tol, int max_iter) {
    int n = A.size();
    if (n == 0) return 0;
    if ((int)A[0].size() != n) throw runtime_error("Matrix A must be square in GMRES.");
    if ((int)b.size() != n || (int)x.size() != n) throw runtime_error("Dimension mismatch in GMRES.");

    int total_steps = 0;
    int max_outer = std::max(1, (max_iter + restart - 1) / restart);

    for (int outer = 0; outer < max_outer; ++outer) {
        // 1) compute residual r0 = b - A*x, beta = ||r0||
        vec r0 = vec_sub(b, mat_vec_mul(A, x));
        double beta = vec_norm(r0);
        if (beta < tol) {
            cout << "solution:-----\n";
            vec_print(x);
            return total_steps;
        }

        // Build v1 = r0 / beta
        vec v1 = scalar_vec_mul(1.0 / beta, r0);

        // 2) Arnoldi to build V (n x (m+1)) and H ((m+1) x m)
        matrix V;
        matrix H;
        int actual_m = 0;
        arnoldi_iteration(A, restart, v1, V, H, actual_m);
        // arnoldi_iteration returns V as n x (actual_m+1) and H as (actual_m+1) x actual_m
        total_steps += actual_m;

        // 3) Form beta*e1 (size actual_m+1)
        vec beta_e1(actual_m + 1, 0.0);
        beta_e1[0] = beta;

        // 4) Solve least squares min || beta_e1 - H * y || by normal equations:
        //    (H^T H) y = H^T (beta_e1)
        int m = actual_m;
        if (m == 0) { // No Krylov vectors generated (rare), treat as stagnation
            break;
        }

        // Compute HtH (m x m) and rhs (m)
        matrix HtH(m, vec(m, 0.0));
        vec rhs(m, 0.0);

        // H is (m+1) x m: H[i][j] (i=0..m, j=0..m-1)
        // HtH[p][q] = sum_{i=0..m} H[i][p] * H[i][q]
        for (int p = 0; p < m; ++p) {
            for (int q = 0; q < m; ++q) {
                double sum = 0.0;
                for (int i = 0; i <= m; ++i) sum += H[i][p] * H[i][q];
                HtH[p][q] = sum;
            }
        }
        // rhs[p] = sum_{i=0..m} H[i][p] * (beta_e1[i])
        for (int p = 0; p < m; ++p) {
            double sum = 0.0;
            for (int i = 0; i <= m; ++i) sum += H[i][p] * beta_e1[i];
            rhs[p] = sum;
        }

        // 5) Solve HtH * y = rhs (dense small system) by Gaussian elimination w/ partial pivoting
        auto solve_dense = [&](matrix M, vec bvec)->vec {
            int N = M.size();
            if ((int)M[0].size() != N) throw runtime_error("solve_dense: need square matrix");
            vec xsol(N, 0.0);

            // forward elimination with partial pivoting
            for (int k = 0; k < N; ++k) {
                // pivot
                int piv = k;
                double maxv = fabs(M[k][k]);
                for (int i = k+1; i < N; ++i) {
                    double av = fabs(M[i][k]);
                    if (av > maxv) { maxv = av; piv = i; }
                }
                if (maxv < 1e-18) throw runtime_error("Singular or ill-conditioned small system in GMRES solve.");
                if (piv != k) {
                    swap(M[k], M[piv]);
                    swap(bvec[k], bvec[piv]);
                }
                // eliminate
                for (int i = k+1; i < N; ++i) {
                    double fac = M[i][k] / M[k][k];
                    for (int j = k; j < N; ++j) M[i][j] -= fac * M[k][j];
                    bvec[i] -= fac * bvec[k];
                }
            }
            // back substitution
            for (int i = N-1; i >= 0; --i) {
                double s = bvec[i];
                for (int j = i+1; j < N; ++j) s -= M[i][j] * xsol[j];
                xsol[i] = s / M[i][i];
            }
            return xsol;
        };

        vec y;
        try {
            y = solve_dense(HtH, rhs);
        } catch (const exception& e) {
            // If solving normal equations fails, fall back to trivial break (stagnation)
            cerr << "GMRES: small system solve failed: " << e.what() << "\n";
            break;
        }

        // 6) Update x = x + V_m * y  (V has n x (m+1) columns, use first m cols)
        for (int j = 0; j < m; ++j) {
            // extract v_j
            vec vj(n);
            for (int i = 0; i < n; ++i) vj[i] = V[i][j];
            x = vec_add(x, scalar_vec_mul(y[j], vj));
        }

        // 7) Check convergence: compute new residual
        r0 = vec_sub(b, mat_vec_mul(A, x));
        beta = vec_norm(r0);
        if (beta < tol) {
            cout << "solution:-----\n";
            vec_print(x);
            return total_steps;
        }

        // continue with next restart
    }

    // finished (max iter reached)
    cout << "solution:-----\n";
    vec_print(x);
    return total_steps;
}