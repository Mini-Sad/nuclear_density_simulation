#include "../include/solver.hpp"
#include <iostream>

Solver::Solver(Basis& basis) : _basis(basis) {}

void Solver::precomputeBasis(const arma::vec& rVals, const arma::vec& zVals) {
    _rVals = rVals;
    _zVals = zVals;

    std::cout << "Pre-calculating Basis Functions..." << std::endl;

    // --- 1. Precompute Z functions ---
    // Find the absolute maximum n_z across all m, n
    int max_nz_overall = _basis.n_zMax.max();

    // Resize matrix: Rows = number of z-points, Cols = max_nz
    _storedZ.set_size(zVals.n_elem, max_nz_overall);

    for (int nz = 0; nz < max_nz_overall; ++nz) {
        _storedZ.col(nz) = _basis.zPart(zVals, nz);
    }

    // --- 2. Precompute R functions ---
    // Field size: (mMax, max_n_overall)
    int max_n_overall = _basis.nMax.max();
    _storedR.set_size(_basis.mMax, max_n_overall);

    for (int m = 0; m < _basis.mMax; ++m) {
        for (int n = 0; n < _basis.nMax(m); ++n) {
            _storedR(m, n) = _basis.rPart(rVals, m, n);
        }
    }
}

arma::mat Solver::calcDensityOptimized(const arma::mat& rho) {
    int nbR = _rVals.n_elem;
    int nbZ = _zVals.n_elem;

    // Initialize Result Density
    arma::mat density = arma::zeros(nbZ, nbR);


    struct QuantumState { int m, n, nz; };
    std::vector<QuantumState> states;

    // Reconstruct the order used in rho matrix (Naive order)
    for (int m = 0; m < _basis.mMax; ++m)
        for (int n = 0; n < _basis.nMax(m); ++n)
            for (int nz = 0; nz < _basis.n_zMax(m, n); ++nz)
                states.push_back({m, n, nz});


    
    std::vector<int> m_start_indices(_basis.mMax + 1, 0);
    int current_idx = 0;
    for (int m = 0; m < _basis.mMax; ++m) {
        m_start_indices[m] = current_idx;
        for (int n = 0; n < _basis.nMax(m); ++n) {
            current_idx += _basis.n_zMax(m, n);
        }
    }
    m_start_indices[_basis.mMax] = current_idx; // End of matrix

    // --- THE OPTIMIZED LOOPS ---
    std::cout << "Running Optimized Calculation..." << std::endl;

    for (int m = 0; m < _basis.mMax; ++m) {



        int start_idx = m_start_indices[m];

        // We need to iterate over pairs of (n_a, n_b)
        for (int n_a = 0; n_a < _basis.nMax(m); ++n_a) {
            for (int n_b = 0; n_b < _basis.nMax(m); ++n_b) {

                // Calculate the "Vertical Part" (Partial Sum over Z)
                // This sum depends only on z.
                arma::vec partialZ = arma::zeros(nbZ);




                int idx_a_base = start_idx;
                for(int k=0; k<n_a; ++k) idx_a_base += _basis.n_zMax(m, k);

                int idx_b_base = start_idx;
                for(int k=0; k<n_b; ++k) idx_b_base += _basis.n_zMax(m, k);

                // Now sum over nz_a and nz_b
                bool block_is_zero = true;

                for (int nz_a = 0; nz_a < _basis.n_zMax(m, n_a); ++nz_a) {
                    for (int nz_b = 0; nz_b < _basis.n_zMax(m, n_b); ++nz_b) {

                        double val = rho(idx_a_base + nz_a, idx_b_base + nz_b);

                        if (std::abs(val) > 1e-10) {
                            block_is_zero = false;
                            
			    
			    partialZ += val * (_storedZ.col(nz_a) % _storedZ.col(nz_b));
                        }
                    }
                }

                // 2. Multiply by Radial Part (if the partial sum wasn't zero)
                if (!block_is_zero) {
                    // product R = R_a * R_b (Element-wise)
                    arma::vec prodR = _storedR(m, n_a) % _storedR(m, n_b);

                    // Add to total density
                    // Outer product: partialZ (col) * prodR (row) -> Matrix
                    density += partialZ * prodR.t();
                }
            }
        }
    }

    return density;
}
