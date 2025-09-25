#include<iostream>
#include<bits/stdc++.h>
#include "includes/library.cpp"

using namespace std;


int main(){
    // int n; // Size of the system
    // ifstream f_info("includes/kinfo.txt");
    // if (!f_info.is_open()) {
    //     cerr << "Could not open kinfo.txt" << endl;
    //     return 1;
    // }
    // f_info >> n;
    // f_info.close();

    // cout << "System size n = " << n << endl;

    // ---- Allocate Memory using C++ vectors ----
    

    // ---- Read Matrix and Vector from files ----
    // The previous code generated these files.
    //IOFF START
    matrix A;
    vec b;
    vec x;
    int n;

    // --- Get matrix A from DIAGONAL storage format ---
    cout << "--- Matrix Reconstruction from Diagonal Storage ---" << endl;
    A = get_A_from_diag_input();
    n = A.size(); // Get the size from the reconstructed matrix

    cout << "\nReconstructed Matrix A:" << endl;
    print_mat(A);

    // --- Get the b vector from the user ---
    b.resize(n);
    x.resize(n);
    cout<<"\n---- Input " << n << " b vector elements ----\n";
    for(int i=0; i<n; ++i) {
    cin >> b[i];
    }
    int r1=n,c1=n;
    //IOFF END

    //GENERAL START
    // int r1,c1;
    // cout<<"enter the A matrix dimensions\n";
    // cin>>r1>>c1;
    // matrix A(r1, vec(c1));
    // vec b(r1);
    // vec x(c1);

    // A=make_mat(r1,c1);
    // // // read_matrix("includes/Kmat.txt", A);
    // // // read_vector("includes/Fvec.txt", b);
    // cout<<"enter the b matrix dimension\n";
    // int r2;
    // cin>>r2;
    // b=make_vec(r2);
    //GENERAL END
    // // ---- Set Solver Parameters ----
    double tolerance = 1e-9;
    int max_iterations = 20000;
    int iterations;

    // cout << "\n--- Solving the System ---" << endl;

    // ---- Jacobi Method ----
    fill(x.begin(), x.end(), 0.0); // Reset initial guess to zeros
    iterations = solve_jacobi(A, b, x, tolerance, max_iterations);
    if (iterations >= max_iterations) {
         cout << "Jacobi failed to converge in " << max_iterations << " iterations." << endl;
    } else {
         cout << "Jacobi converged in " << iterations << " iterations." << endl;
    }


    // // ---- Gauss-Seidel Method ----
    // fill(x.begin(), x.end(), 0.0); // Reset initial guess
    // iterations = solve_gs(A, b, x, tolerance, max_iterations);
    // if (iterations >= max_iterations) {
    //      cout << "Gauss-Seidel failed to converge in " << max_iterations << " iterations." << endl;
    // } else {
    //      cout << "Gauss-Seidel converged in " << iterations << " iterations." << endl;
    // }

    // ---- SOR Method ----
    // vector<double> omegas;
    // double l=1;
    // for(int i=0;i<10;i++)
    // {
    //     l+=0.1;
    //     omegas.push_back(l);
    // }
    // int num_omegas = omegas.size();
    // for (int i = 0; i < num_omegas; i++) {
    //     fill(x.begin(), x.end(), 0.0); // Reset initial guess
    //     iterations = solve_sor(A, b, x, omegas[i], tolerance, max_iterations);
    //     if (iterations >= max_iterations) {
    //          cout << "SOR (omega=" << omegas[i] << ") failed to converge in " << max_iterations << " iterations." << endl;
    //     } else {
    //          cout << "SOR (omega=" << omegas[i] << ") converged in " << iterations << " iterations." << endl;
    //     }
    // }

    // // ---- Prepare for Spectral Radius Calculation ----
    // cout << "\n--- Generating Iteration Matrix for Spectral Analysis ---" << endl;
    // matrix T_J(r1, vec(c1));

    // create_jacobi_iteration_matrix(A, T_J);
    // write_matrix("T_Jacobi.txt", T_J);
    // cout << "Jacobi iteration matrix written to T_Jacobi.txt" << endl;


    //make vec->text->python svd
    // int r,c;
    // cin>>r>>c;
    // matrix A=make_mat(r,c);
    // text_from_mat(A);

    //for steepest descent
    // fill(x.begin(), x.end(), 0.0); // Reset initial guess to zeros
    // iterations = solve_steepest_descent(A, b, x, tolerance, max_iterations);
    // if (iterations >= max_iterations) {
    //      cout << "Steepest failed to converge in " << max_iterations << " iterations." << endl;
    // } else {
    //      cout << "Steepest converged in " << iterations << " iterations." << endl;
    // }
    // fill(x.begin(), x.end(), 0.0); // Reset initial guess to zeros
    // iterations = solve_minimum_residual(A, b, x, tolerance, max_iterations);
    // if (iterations >= max_iterations) {
    //      cout << "solve_minimum_residual failed to converge in " << max_iterations << " iterations." << endl;
    // } else {
    //      cout << "solve_minimum_residualconverged in " << iterations << " iterations." << endl;
    // }
    // fill(x.begin(), x.end(), 0.0); // Reset initial guess to zeros
    // iterations = solve_residual_norm_steepest_descent(A, b, x, tolerance, max_iterations);
    // if (iterations >= max_iterations) {
    //      cout << "solve_residual_norm_steepest_descent failed to converge in " << max_iterations << " iterations." << endl;
    // } else {
    //      cout << "solve_residual_norm_steepest_descent converged in " << iterations << " iterations." << endl;
    // }

    // for(int m=1;m<=r1;m++)//Try till all dimensions of the parent vector space 
    // {   fill(x.begin(), x.end(), 0.0); // Reset initial guess to zeros
    //     int iters = solve_fom(A, b, x,m, tolerance, max_iterations);

    //     if (iters < max_iterations) {
    //         cout <<"m= "<<m<<" Converged in " << iters << " iteration(s)." << endl;
           
            
    //     } else {
    //         cout <<"m= "<<m<< "FOM did not converge within " << max_iterations << " iterations." <<endl;
    //     }
    // }
    
    //reset the intial guess x0 each time before running the method
    // fill(x.begin(), x.end(), 0.0);
    // int iters=solve_conjugate_gradient(A, b, x, tolerance, max_iterations);
    // if (iters < max_iterations) {
    //     cout << " Converged in " << iters << " iteration(s)." << endl;
        
        
    // } else {
    //     cout<< "FOM did not converge within " << max_iterations << " iterations." <<endl;
    // }
    return 0;
}
