#ifndef BASIS_HPP
#define BASIS_HPP

#include <armadillo>
#include "poly.hpp"

class Basis {
private:
    // Physical parameters
    double _br;
    double _bz;
    int _N;
    double _Q;

    // Polynomial calculator
    Poly _poly;

    // Pre-calculated factorials (stored as double for calculations)
    arma::vec _fact; 

    // Helper to calculate truncation function nu(i)
    double valNu(int i) const;

    // Helper to pre-calculate factorials up to a safe limit
    void initFactorials();

public:
    // Public Truncation Parameters (as required by Test #01)
    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;

    /**
     * @brief Constructor
     * @param br Radial deformation parameter
     * @param bz Axial deformation parameter
     * @param N  Truncation parameter N
     * @param Q  Truncation parameter Q
     */
    Basis(double br, double bz, int N, double Q);

    /**
     * @brief Calculate the Z-part of the basis function: Z(z, n_z)
     * @param z   Vector of z-coordinates
     * @param n_z Quantum number n_z
     * @return    Vector of function values
     */
    arma::vec zPart(const arma::vec& z, int n_z);

    /**
     * @brief Calculate the R-part of the basis function: R(r, m, n)
     * @param r   Vector of r-coordinates
     * @param m   Quantum number m
     * @param n   Quantum number n
     * @return    Vector of function values
     */
    arma::vec rPart(const arma::vec& r, int m, int n);
};

#endif // BASIS_HPP
