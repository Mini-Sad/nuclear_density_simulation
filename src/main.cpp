#include "../include/poly.h"
#include "../include/basis.h"
#include <stdio.h>


int main(){
    Poly poly;
    arma::vec zVals;

    zVals = {-3.1, -2.3, -1.0, -0.3, 0.1, 4.3, 9.2, 13.7};

    poly.calcHermite(6, zVals);

    std::cout << "poly.hermite(4)" << " ========= " ;
    std::cout << std::endl <<  poly.hermite(4) << std::endl;

    std::cout << "-------------------------------------" << std::endl;


    Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);
    std::cout << "mMax" << " ========= " ;
    std::cout << std::endl <<  basis.mMax << std::endl;
    std::cout << "nMax" << " ========= " ;
    std::cout << std::endl <<  basis.nMax << std::endl;
    std::cout << "n_zMax" << " ========= " ;
    std::cout << std::endl <<  basis.n_zMax << std::endl;

    std::cout << "-------------------------------------" << std::endl;

    arma::vec r = {3.1, 2.3, 1.0, 0.0, 0.1, 4.3, 9.2, 13.7};

    poly.calcLaguerre(fabs(8)+1,2+1,r%r/pow(1.935801664793151,2));
    std::cout << "8.2" << " ========= " ;
    std::cout << std::endl <<  poly.laguerre(fabs(8),2)<< std::endl;
    std::cout << "0.0" << " ========= " ;
    std::cout << std::endl <<  poly.laguerre(fabs(0),0)<< std::endl;

    std::cout << "r_function00" << " ========= " ;
    std::cout << std::endl <<  basis.rPart(r,0,0)<< std::endl;
    std::cout << "r_function82" << " ========= " ;
    std::cout << std::endl <<  basis.rPart(r,8,2)<< std::endl;





    

    return 0;
}
