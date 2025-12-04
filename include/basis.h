#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include "poly.h"
#include <algorithm>

using namespace std;

class Basis{

    public :
    double br;
    double bz;
    int N;
    double Q;
    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;
    const double pi = M_PI; 

    //Basis  r-functions
    arma::vec r00;

    //Basis f-functions
    arma::vec f;




    public :
    Basis(double _br, double _bz, int _N, double _Q);

    double getnu(int , int , double );

    void setmMax(void);

    void setnMax(void);

    void setnZMax(void);

    int getmMax(void);

    arma::ivec getnMax(void);

    arma::imat getnZMax(void);

    arma::vec rPart(arma::vec, int, int);

};

#endif // BASIS_H