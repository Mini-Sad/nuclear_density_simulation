#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include "poly.h"

using namespace std;

class Basis{

    private :
    double bz;
    double br;
    int N;
    double Q;
    int max;
    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;
    arma::vec vVals;

    //Basis  r-functions
    arma::vec r;

    //Basis f-functions
    arma::vec f;




    public :
    Basis(double _br, double _bz, int _N, double _Q);


    void calcV(void);

    int getmax(void);

    int getmMax(void);

    arma::ivec getnMax(void);

    arma::imat getnZMax(void);

    arma::vec rPart(arma::vec, int, int);

};

#endif // BASIS_H