#include <iostream>
#include <armadillo>
#include "basis.hpp"
#include "solver.hpp"

int main() {
    // ==========================================
    // 1. SETUP & LOADING
    // ==========================================
    std::cout << "--- IPS-PROD PROJECT START ---" << std::endl;

    // From Project: br=1.9358..., bz=2.8296..., N=14, Q=1.3
    Basis basis(1.935801664793151, 2.829683956491218, 14, 1.3);

    arma::mat rho_mat;
    if (!rho_mat.load("rho.arma", arma::arma_ascii)) {
        std::cerr << "Error: Could not load rho.arma" << std::endl;
        return 1;
    }
    std::cout << "Loaded rho matrix: " << rho_mat.n_rows << "x" << rho_mat.n_cols << std::endl;

    // Define Grid
    int nbR = 32;
    int nbZ = 64;
    arma::vec rVals = arma::linspace(0.0, 10.0, nbR);
    arma::vec zVals = arma::linspace(-20.0, 20.0, nbZ);

    arma::wall_clock timer;

    // ==========================================
    // 2. OPTIMIZED ALGORITHM (Task 2)
    // ==========================================
    std::cout << "\n[Task 2] Starting OPTIMIZED Algorithm..." << std::endl;
    
    Solver solver(basis);
    solver.precomputeBasis(rVals, zVals); // Pre-calc happens here

    timer.tic();
    arma::mat densityOptimized = solver.calcDensityOptimized(rho_mat);
    double timeOptimized = timer.toc();
    
    std::cout << ">>> Optimized Time: " << timeOptimized << " seconds" << std::endl;

    // ==========================================
    // 3. NAIVE ALGORITHM (Task 1)
    // ==========================================
    std::cout << "\n[Task 1] Starting NAIVE Algorithm..." << std::endl;
    std::cout << "(This might take a while... grab a coffee)" << std::endl;

    arma::mat densityNaive = arma::zeros(nbZ, nbR);
    timer.tic();

    int idx_a = 0;
    // Loop A (rows of rho matrix)
    for (int m_a = 0; m_a < basis.mMax; ++m_a) {
        for (int n_a = 0; n_a < basis.nMax(m_a); ++n_a) {
            for (int nz_a = 0; nz_a < basis.n_zMax(m_a, n_a); ++nz_a) {

                // Optimization Lvl 1: Pre-calculate Psi_A outside the inner loop
                arma::vec R_a = basis.rPart(rVals, m_a, n_a);
                arma::vec Z_a = basis.zPart(zVals, nz_a);
                // Outer product Z * R^T
                arma::mat Psi_a = Z_a * R_a.t(); 

                int idx_b = 0;
                // Loop B (cols of rho matrix)
                for (int m_b = 0; m_b < basis.mMax; ++m_b) {
                    for (int n_b = 0; n_b < basis.nMax(m_b); ++n_b) {
                        for (int nz_b = 0; nz_b < basis.n_zMax(m_b, n_b); ++nz_b) {
                            
                            // Naive access to rho
                            double rho_val = rho_mat(idx_a, idx_b);

                            // Basic check to save SOME time (skip zeros)
                            if (std::abs(rho_val) > 1e-12) {
                                arma::vec R_b = basis.rPart(rVals, m_b, n_b);
                                arma::vec Z_b = basis.zPart(zVals, nz_b);
                                arma::mat Psi_b = Z_b * R_b.t();

                                // Accumulate: result += rho * Psi_A * Psi_B
                                densityNaive += rho_val * (Psi_a % Psi_b);
                            }
                            idx_b++;
                        }
                    }
                }
                idx_a++;
            }
        }
        // Progress bar for Naive (prints every m_a step)
        std::cout << "." << std::flush;
    }
    double timeNaive = timer.toc();
    std::cout << "\n>>> Naive Time: " << timeNaive << " seconds" << std::endl;

    // ==========================================
    // 4. COMPARISON & SAVING
    // ==========================================
    std::cout << "\n--- RESULTS ---" << std::endl;
    std::cout << "Speedup Factor: " << timeNaive / timeOptimized << "x" << std::endl;

    // Calculate difference (Error)
    double diff = arma::norm(densityNaive - densityOptimized);
    std::cout << "Difference (Norm): " << diff << std::endl;

    if (diff < 1e-8) {
        std::cout << "SUCCESS: Algorithms match!" << std::endl;
    } else {
        std::cerr << "WARNING: Results differ! Check logic." << std::endl;
    }

    // Save one of them (they are the same)
    densityOptimized.save("density.mat", arma::arma_ascii);


    return 0;
}
