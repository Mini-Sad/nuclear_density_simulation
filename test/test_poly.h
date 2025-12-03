#include <cxxtest/TestSuite.h>
#include "../include/poly.h"
#include "../include/basis.h"

class TestPoly:public CxxTest::TestSuite
{
    // Mandatory test #00 - Hermite and Laguerre polynomials
    public:
    void testPoly(void){
        Poly poly;

        arma::vec zVals, calcVals, targetVals;

        zVals = {-3.1, -2.3, -1.0, -0.3, 0.1, 4.3, 9.2, 13.7};
        poly.calcHermite(6, zVals); // compute Hermite polynomials for n in {0 ... 5}

        calcVals   = poly.hermite(4); // n = 4
        targetVals = {  1.02835360e+03,  2.05825600e+02, -2.00000000e+01,  7.80960000e+00,
                        1.15216000e+01,  4.59456160e+03,  1.10572154e+05,  5.54643458e+05};

        TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);

        calcVals   = poly.hermite(5); // n = 5
        targetVals = { -4.76676832e+03, -3.88909760e+02,  8.00000000e+00, -3.17577600e+01,
                        1.18403200e+01,  3.48375818e+04,  1.98557479e+06,  1.50339793e+07};
        TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);

        zVals = {0.1, 0.3, 1.2, 1.8, 2.0, 2.5, 7.1, 11.1};
        poly.calcLaguerre(6, 4, zVals); // compute generalized Laguerre polynomials for m in {0 ... 5} and n in {0 ... 3}

        calcVals   = poly.laguerre(4, 2); // m = 4, n = 2
        targetVals = {  14.405,  13.245,  8.52 ,  5.82 ,  5.,  3.125,  -2.395,  10.005};
        TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);
        
        calcVals   = poly.laguerre(5, 3); // m = 5, n = 3
        targetVals = { 53.23983333,  47.95550000,  27.87200000,  17.5880,
                    14.66666667,   8.39583333,  -0.81183333,  10.1015};
        TS_ASSERT_DELTA(arma::norm(calcVals / targetVals - 1.0), 0.0, 1e-08);    
    }    
    
    void TestBasis(void)
    {
        // Mandatory test #01 - Basis truncation
        //     br = 1.935801664793151, bz = 2.829683956491218, N = 14, Q = 1.3
        Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);

        TS_ASSERT_EQUALS(basis.mMax, 14);

        arma::ivec nMax = {7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1};

        TS_ASSERT((basis.nMax - nMax).is_zero());

        arma::imat n_zMax = {{18, 15, 13, 10, 7, 5, 2}, 
                            {16, 14, 11,  9, 6, 3, 1}, 
                            {15, 13, 10,  7, 5, 2, 0}, 
                            {14, 11,  9,  6, 3, 1, 0}, 
                            {13, 10,  7,  5, 2, 0, 0}, 
                            {11,  9,  6,  3, 1, 0, 0}, 
                            {10,  7,  5,  2, 0, 0, 0}, 
                            { 9,  6,  3,  1, 0, 0, 0}, 
                            { 7,  5,  2,  0, 0, 0, 0}, 
                            { 6,  3,  1,  0, 0, 0, 0}, 
                            { 5,  2,  0,  0, 0, 0, 0}, 
                            { 3,  1,  0,  0, 0, 0, 0}, 
                            { 2,  0,  0,  0, 0, 0, 0}, 
                            { 1,  0,  0,  0, 0, 0, 0}};

        // check if matrices are equal
        TS_ASSERT((basis.n_zMax - n_zMax).is_zero());
    }

};