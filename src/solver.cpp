/**
 * @file solver.cpp
 * @brief Implementation of the Solver class for nuclear density calculation.
 * * This file contains the logic for both pre-calculating the basis functions
 * and running the optimized algorithm using partial summation.
 */

#include "solver.hpp"
#include <iostream>

Solver::Solver(Basis& basis) : _basis(basis) {}

void Solver::precomputeBasis(const arma::vec& rVals, const arma::vec& zVals) {
    _rVals = rVals;
    _zVals = zVals;

    std::cout << "Pre-calculating Basis Functions..." << std::endl;

    // ==========================================
    // 1. Precompute Vertical (Z) Functions
    // ==========================================
    // The Z-part depends only on n_z. We find the maximum n_z used across all states
    // so we can size our matrix correctly.
    int max_nz_overall = _basis.n_zMax.max();
    
    // Allocate matrix: Rows = z-points, Cols = n_z
    _storedZ.set_size(zVals.n_elem, max_nz_overall);

    for (int nz = 0; nz < max_nz_overall; ++nz) {
        // We calculate the vector once and store it. 
        // Later, we just look it up in O(1) time.
        _storedZ.col(nz) = _basis.zPart(zVals, nz);
    }

    // ==========================================
    // 2. Precompute Radial (R) Functions
    // ==========================================
    // The R-part depends on m and n.
    // We allocate a Field (2D array of vectors) of size (mMax, max_n)
    int max_n_overall = _basis.nMax.max();
    _storedR.set_size(_basis.mMax, max_n_overall);

    for (int m = 0; m < _basis.mMax; ++m) {
        for (int n = 0; n < _basis.nMax(m); ++n) {
            // Store R_{mn}(r)
            _storedR(m, n) = _basis.rPart(rVals, m, n);
        }
    }
}

arma::mat Solver::calcDensityOptimized(const arma::mat& rho) {
    int nbR = _rVals.n_elem;
    int nbZ = _zVals.n_elem;
    
    // Result matrix: Rows correspond to Z, Columns correspond to R
    arma::mat density = arma::zeros(nbZ, nbR);

    // ==========================================
    // Index Mapping Strategy
    // ==========================================
    // The rho matrix is flat (2D), but physics is 6D.
    // To access rho(a, b) efficiently without complex math inside the loop,
    // we pre-calculate the linear starting index for each 'm' block.
    // The structure is: Block(m=0) -> Block(m=1) -> ...
    
    std::vector<int> m_start_indices(_basis.mMax + 1, 0);
    int current_idx = 0;
    for (int m = 0; m < _basis.mMax; ++m) {
        m_start_indices[m] = current_idx; // Start of block m
        for (int n = 0; n < _basis.nMax(m); ++n) {
            current_idx += _basis.n_zMax(m, n); // Add number of nz states
        }
    }
    m_start_indices[_basis.mMax] = current_idx; // End of matrix

    std::cout << "Running Optimized Calculation..." << std::endl;

    // ==========================================
    // The Optimized Loop Structure
    // ==========================================
    // Assumption: Density matrix is m-diagonal (rho_{ab} = 0 if m_a != m_b).
    // We loop over 'm' once, effectively diagonalizing the problem.
    for (int m = 0; m < _basis.mMax; ++m) {
        
        int start_idx = m_start_indices[m];
        
        // Loop over radial quantum numbers n_a and n_b
        for (int n_a = 0; n_a < _basis.nMax(m); ++n_a) {
            for (int n_b = 0; n_b < _basis.nMax(m); ++n_b) {
                
                // --- STEP 1: Partial Sum over Vertical (Z) States ---
                // We factorize the sum: Sum_all = Sum_R * ( Sum_Z )
                // This vector 'partialZ' will hold the result of the inner Z sum.
                arma::vec partialZ = arma::zeros(nbZ);
                
                // Calculate the exact linear offset for n_a and n_b within this m-block
                int idx_a_base = start_idx;
                for(int k=0; k<n_a; ++k) idx_a_base += _basis.n_zMax(m, k);
                
                int idx_b_base = start_idx;
                for(int k=0; k<n_b; ++k) idx_b_base += _basis.n_zMax(m, k);

                bool block_is_zero = true;
                
                // Inner loops over vertical quantum numbers (nz)
                for (int nz_a = 0; nz_a < _basis.n_zMax(m, n_a); ++nz_a) {
                    for (int nz_b = 0; nz_b < _basis.n_zMax(m, n_b); ++nz_b) {
                        
                        // Access the density matrix element using our pre-calculated offsets
                        double val = rho(idx_a_base + nz_a, idx_b_base + nz_b);
                        
                        // Optimization: Skip trivial zeros (sparse matrix optimization)
                        if (std::abs(val) > 1e-10) {
                            block_is_zero = false;
                            // Accumulate the weighted product of Z functions.
                            // We use element-wise multiplication (%) which is vectorized.
                            partialZ += val * (_storedZ.col(nz_a) % _storedZ.col(nz_b));
                        }
                    }
                }

                // --- STEP 2: Combine with Radial (R) Part ---
                // If the partial sum was all zeros, we skip the expensive outer product.
                if (!block_is_zero) {
                    // product R = R_a(r) * R_b(r) (Element-wise)
                    arma::vec prodR = _storedR(m, n_a) % _storedR(m, n_b);
                    
                    // Add to total density using Outer Product:
                    // partialZ (col vec) * prodR^T (row vec) -> Matrix(nbZ, nbR)
                    // This updates the entire 2D density grid in one operation.
                    density += partialZ * prodR.t();
                }
            }
        }
    }
    
    return density;
}