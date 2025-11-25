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




    public :
    Basis(double _br, double _bz, int _N, double _Q);

};

#endif // BASIS_H