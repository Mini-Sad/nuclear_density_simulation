#include "../include/poly.hpp"

Poly::Poly() {
    // Constructor (empty for now)
}

// ---------------------------------------------------------
// Hermite Polynomials
// Recurrence: H_n(z) = 2*z * H_{n-1}(z) - 2*(n-1) * H_{n-2}(z)
// ---------------------------------------------------------
void Poly::calcHermite(int nMax, const arma::vec& z) {
    int numPoints = z.n_elem;
    
    // Resize matrix: Rows = points, Cols = order n
    // This makes retrieving a specific 'n' contiguous in memory
    _hermiteVals.set_size(numPoints, nMax);

    // Order n = 0: H_0(z) = 1
    if (nMax > 0) {
        _hermiteVals.col(0).ones();
    }

    // Order n = 1: H_1(z) = 2z
    if (nMax > 1) {
        _hermiteVals.col(1) = 2.0 * z;
    }

    // Recurrence for n >= 2
    for (int n = 2; n < nMax; ++n) {
        // H_n = 2zH_{n-1} - 2(n-1)H_{n-2}
        // Note the use of '%' for element-wise multiplication
        _hermiteVals.col(n) = (2.0 * z % _hermiteVals.col(n - 1)) 
                              - (2.0 * (n - 1) * _hermiteVals.col(n - 2));
    }
}

arma::vec Poly::hermite(int n) {
    return _hermiteVals.col(n);
}

// ---------------------------------------------------------
// Generalized Laguerre Polynomials
// Recurrence: 
// L^m_0 = 1
// L^m_1 = 1 + m - z
// L^m_n = (2 + (m-1-z)/n) * L^m_{n-1} - (1 + (m-1)/n) * L^m_{n-2}
// ---------------------------------------------------------
void Poly::calcLaguerre(int mMax, int nMax, const arma::vec& z) {
    int numPoints = z.n_elem;

    // Resize Cube: (Rows=points, Cols=m, Slices=n)
    _laguerreVals.set_size(numPoints, mMax, nMax);

    // We iterate over every 'm' first, as the recurrence is on 'n' for a fixed 'm'
    for (int m = 0; m < mMax; ++m) {
        
        // Order n = 0: L^m_0(z) = 1
        if (nMax > 0) {
            _laguerreVals.slice(0).col(m).ones();
        }

        // Order n = 1: L^m_1(z) = 1 + m - z
        if (nMax > 1) {
            _laguerreVals.slice(1).col(m) = (1.0 + m) - z;
        }

        // Recurrence for n >= 2
        for (int n = 2; n < nMax; ++n) {
            // Term A: (2 + (m - 1 - z) / n)
            // We compute this term vectorised over z
            arma::vec termA = 2.0 + ((m - 1.0 - z) / (double)n);
            
            // Term B: (1 + (m - 1) / n)
            double termB = 1.0 + ((double)(m - 1) / (double)n);

            // Apply formula
            // L_n = TermA * L_{n-1} - TermB * L_{n-2}
            arma::vec L_prev1 = _laguerreVals.slice(n - 1).col(m);
            arma::vec L_prev2 = _laguerreVals.slice(n - 2).col(m);

            _laguerreVals.slice(n).col(m) = (termA % L_prev1) - (termB * L_prev2);
        }
    }
}

arma::vec Poly::laguerre(int m, int n) {
    // Access the slice for order n, then the column for index m
    return _laguerreVals.slice(n).col(m);
}
