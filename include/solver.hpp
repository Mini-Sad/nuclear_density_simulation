#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <armadillo>
#include "basis.h"

class Solver {
private:
    Basis& _basis; 

    // Pre-calculated storage
    // R_funcs(m, n) -> vector of size nbR
    arma::field<arma::vec> _storedR;

    // Z_funcs column nz -> vector of size nbZ
    arma::mat _storedZ;

    // Grid Mapping (needed for precalculation)
    arma::vec _rVals;
    arma::vec _zVals;

public:
    Solver(Basis& basis);

    //  Pre-calculate all basis functions on the given grid
    void precomputeBasis(const arma::vec& rVals, const arma::vec& zVals);

    //  Run the Optimized Algorithm
    // Returns the density matrix (Rows=Z, Cols=R)
    arma::mat calcDensityOptimized(const arma::mat& rho);
};

#endif
