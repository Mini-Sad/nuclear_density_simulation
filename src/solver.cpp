#include "../include/solver.h"
#include <iostream>
#include <cmath>

Solver::Solver(Basis& basis) : _basis(basis) {}

void Solver::precomputeBasis(const arma::vec& rVals, const arma::vec& zVals) {
    _rVals = rVals;
    _zVals = zVals;

    std::cout << "Pre-calculating Basis Functions..." << std::endl;

    // --- 1. Precompute Z functions ---
    int max_nz_overall = _basis.n_zMax.max();
    _storedZ.set_size(zVals.n_elem, max_nz_overall);

    for (int nz = 0; nz < max_nz_overall; ++nz) {
        _storedZ.col(nz) = _basis.zPart(zVals, nz);
    }

    // --- 2. Precompute R functions ---
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

    // Map pour trouver rapidement l'index global de départ pour chaque m
    std::vector<int> m_start_indices(_basis.mMax + 1, 0);
    int current_idx = 0;
    for (int m = 0; m < _basis.mMax; ++m) {
        m_start_indices[m] = current_idx;
        for (int n = 0; n < _basis.nMax(m); ++n) {
            current_idx += _basis.n_zMax(m, n);
        }
    }
    m_start_indices[_basis.mMax] = current_idx; // End of matrix

    std::cout << "Running Optimized Calculation (Full Cross-Terms)..." << std::endl;

    // --- BOUCLE PRINCIPALE (Corrigée pour inclure m_a != m_b) ---
    // On doit sommer sur m_a ET m_b pour correspondre à l'algo naïf (coupe phi=0)
    
    for (int m_a = 0; m_a < _basis.mMax; ++m_a) {
        int start_idx_a = m_start_indices[m_a];

        for (int m_b = 0; m_b < _basis.mMax; ++m_b) {
            int start_idx_b = m_start_indices[m_b];

            // Boucle sur n_a et n_b
            for (int n_a = 0; n_a < _basis.nMax(m_a); ++n_a) {
                
                // Calcul de l'offset précis pour le bloc (m_a, n_a)
                int idx_a_base = start_idx_a;
                for(int k=0; k<n_a; ++k) idx_a_base += _basis.n_zMax(m_a, k);

                for (int n_b = 0; n_b < _basis.nMax(m_b); ++n_b) {
                    
                    // Calcul de l'offset précis pour le bloc (m_b, n_b)
                    int idx_b_base = start_idx_b;
                    for(int k=0; k<n_b; ++k) idx_b_base += _basis.n_zMax(m_b, k);

                    // 1. Partie Verticale (Z)
                    arma::vec partialZ = arma::zeros(nbZ);
                    bool block_is_zero = true;

                    // Somme sur nz_a et nz_b
                    for (int nz_a = 0; nz_a < _basis.n_zMax(m_a, n_a); ++nz_a) {
                        for (int nz_b = 0; nz_b < _basis.n_zMax(m_b, n_b); ++nz_b) {

                            // On récupère rho(global_idx_a, global_idx_b)
                            double val = rho(idx_a_base + nz_a, idx_b_base + nz_b);

                            if (std::abs(val) > 1e-12) {
                                block_is_zero = false;
                                // Accumulation Z : val * Z_a * Z_b
                                partialZ += val * (_storedZ.col(nz_a) % _storedZ.col(nz_b));
                            }
                        }
                    }

                    // 2. Partie Radiale (R) et ajout au total
                    if (!block_is_zero) {
                        // Produit R = R_a * R_b (element-wise)
                        // Note: Ici on suppose phi=0, donc les termes exponentiels imaginaires s'annulent ou valent 1
                        arma::vec prodR = _storedR(m_a, n_a) % _storedR(m_b, n_b);

                        // Outer product: partialZ (col) * prodR (row) -> Matrix
                        density += partialZ * prodR.t();
                    }
                }
            }
        }
    }

    return density;
}