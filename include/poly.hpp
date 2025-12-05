#ifndef POLY_HPP
#define POLY_HPP

#include <armadillo>

class Poly {
private:
    // Storage for Hermite polynomials
    // Structure: Rows = coordinate points, Cols = order n
    // This allows efficient column extraction: _hermiteVals.col(n)
    arma::mat _hermiteVals;

    // Storage for Generalized Laguerre polynomials
    // Structure: A 3D Cube
    // Dimension 0 (Rows): coordinate points (z/eta)
    // Dimension 1 (Cols): index m
    // Dimension 2 (Slices): order n
    // Access: _laguerreVals.slice(n).col(m)
    arma::cube _laguerreVals;

public:
    // Default constructor
    Poly();

    /**
     * @brief Calculate Hermite polynomials H_n(z)
     * @param nMax The maximum order (exclusive, calculates 0 to nMax-1)
     * @param z The vector of points to evaluate
     */
    void calcHermite(int nMax, const arma::vec& z);

    /**
     * @brief Retrieve calculated Hermite polynomial
     * @param n The order
     * @return Vector of values H_n(z)
     */
    arma::vec hermite(int n);

    /**
     * @brief Calculate Generalized Laguerre polynomials L^m_n(z)
     * @param mMax The maximum m index (exclusive)
     * @param nMax The maximum n order (exclusive)
     * @param z The vector of points to evaluate
     */
    void calcLaguerre(int mMax, int nMax, const arma::vec& z);

    /**
     * @brief Retrieve calculated Laguerre polynomial
     * @param m The upper index
     * @param n The lower index (order)
     * @return Vector of values L^m_n(z)
     */
    arma::vec laguerre(int m, int n);
};

#endif // POLY_HPP
