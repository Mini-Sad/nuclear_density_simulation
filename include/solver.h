#ifndef SOLVER_H
#define SOLVER_H

#include <armadillo>
#include "basis.h" 

/**
 * @class Solver
 * @brief Handles the calculation of the nuclear local density.
 * * This class implements an optimized algorithm to calculate the density
 * \f$ \rho(\mathbf{r}) \f$ using the partial sum factorization method.
 * It manages the pre-calculation of basis functions to avoid redundant evaluations.
 */
class Solver {
private:
    /// Reference to the physics engine (Basis) containing parameters and polynomials.
    Basis& _basis; 

    // --- Internal Cache ---

    /**
     * @brief Pre-calculated Radial functions R(r).
     * * Stored as a field of vectors.
     * Access: _storedR(m, n) returns the vector for the whole grid rVals.
     */
    arma::field<arma::vec> _storedR;
    
    /**
     * @brief Pre-calculated Vertical functions Z(z).
     * * Stored as a matrix where each column corresponds to a quantum number n_z.
     * Access: _storedZ.col(nz) returns the vector for the whole grid zVals.
     */
    arma::mat _storedZ; 

    /// The radial grid points used for pre-calculation.
    arma::vec _rVals;
    /// The vertical grid points used for pre-calculation.
    arma::vec _zVals;

public:
    /**
     * @brief Constructor.
     * @param basis Reference to an initialized Basis object.
     */
    Solver(Basis& basis);

    /**
     * @brief Pre-calculates all basis functions on the given grid.
     * * This method iterates over all valid quantum numbers defined in the Basis
     * and stores the result of rPart() and zPart() in internal memory.
     * This must be called before calcDensityOptimized().
     * * @param rVals Vector of radial coordinates (r).
     * @param zVals Vector of vertical coordinates (z).
     */
    void precomputeBasis(const arma::vec& rVals, const arma::vec& zVals);

    /**
     * @brief Calculates the nuclear local density using the Optimized Algorithm.
     * * Uses the "Partial Sum" strategy:
     * \f[
     * \rho(r, z) = \sum_{m} \sum_{n_a, n_b} R_{a}(r)R_{b}(r) \left( \sum_{n_{za}, n_{zb}} \rho_{ab} Z_{a}(z)Z_{b}(z) \right)
     * \f]
     * This reduces complexity by factoring out the Z-summation.
     * * @param rho The input density matrix (must be loaded from file).
     * @return A matrix representing the density \f$ \rho(z, r) \f$ (Rows=Z, Cols=R).
     */
    arma::mat calcDensityOptimized(const arma::mat& rho);
};

#endif // SOLVER_HPP