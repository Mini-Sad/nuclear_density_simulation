#include "basis.hpp"
#include <cmath> // for tgamma, pow, sqrt, exp

Basis::Basis(double br, double bz, int N, double Q)
    : _br(br), _bz(bz), _N(N), _Q(Q) {
    
    // 1. Initialize Factorials
    initFactorials();

    // 2. Determine Truncation Limits (Test #01)
    
    // Calculate mMax: sup { i : nu(i) >= 1 }
    // We increment i until it becomes invalid
    int i = 0;
    while (valNu(i) >= 1.0) {
        i++;
    }
    // i is now the first invalid index (e.g., 15)
    // The Supremum (highest valid index) is i - 1 (e.g., 14)
    // The project convention for mMax seems to be this Supremum value.
    mMax = i - 1;

    // Calculate nMax[m]
    nMax.set_size(mMax);
    for (int m = 0; m < mMax; ++m) {
        // Formula: nMax = 1/2 * (mMax - m - 1) + 1
        nMax(m) = (int)(0.5 * (mMax - m - 1)) + 1;
    }

    // Calculate n_zMax[m, n]
    int max_n_overall = nMax.max();
    n_zMax.set_size(mMax, max_n_overall); 
    n_zMax.zeros(); // <--- CRITICAL FIX: Initialize empty spots to 0

    for (int m = 0; m < mMax; ++m) {
        for (int n = 0; n < nMax(m); ++n) {
            // Formula: n_zMax = nu(m + 2n + 1)
            double nu_val = valNu(m + 2 * n + 1);
            n_zMax(m, n) = (int)nu_val;
        }
    }
}

double Basis::valNu(int i) const {
    // Note: 2.0 / 3.0 ensures floating point division
    return (_N + 2.0) * std::pow(_Q, 2.0 / 3.0) + 0.5 - i * _Q;
}

void Basis::initFactorials() {
    // Pre-calculate factorials up to a generous limit
    int limit = 60; 
    _fact.set_size(limit);
    for(int i=0; i<limit; ++i) {
        _fact(i) = std::tgamma(i + 1);
    }
}

// ---------------------------------------------------------
// Z-Function
// Z(z, n_z) = Norm * Gaussian * H_{n_z}(z/b_z)
// ---------------------------------------------------------
arma::vec Basis::zPart(const arma::vec& z, int n_z) {
    arma::vec zeta = z / _bz;

    _poly.calcHermite(n_z + 1, zeta);
    arma::vec H_val = _poly.hermite(n_z);

    double denom = _bz * std::pow(2.0, n_z) * std::sqrt(M_PI) * _fact(n_z);
    double norm = 1.0 / std::sqrt(denom);

    arma::vec gaussian = arma::exp(-0.5 * arma::square(zeta));

    return norm * gaussian % H_val;
}

// ---------------------------------------------------------
// R-Function
// R(r, m, n) = Norm * Gaussian * Power * L_n^m(r^2/b^2)
// ---------------------------------------------------------
arma::vec Basis::rPart(const arma::vec& r, int m, int n) {
    arma::vec r_scaled = r / _br;
    arma::vec eta = arma::square(r_scaled); 

    _poly.calcLaguerre(m + 1, n + 1, eta);
    arma::vec L_val = _poly.laguerre(m, n);

    double part1 = 1.0 / (_br * std::sqrt(M_PI));
    double part2 = std::sqrt(_fact(n) / _fact(n + m));
    double norm = part1 * part2;

    arma::vec gaussian = arma::exp(-0.5 * eta);

    arma::vec powerTerm = arma::pow(r_scaled, (double)m);

    return norm * gaussian % powerTerm % L_val;
}
