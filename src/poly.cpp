#include "../include/poly.h"

//Constructeur
Poly::Poly(){


    
}



void Poly::calcHermite(int n, arma::vec zVals){
    arma::vec resultat_hermite = arma::ones(zVals.size());
    if(n>0)
    {
        resultat_hermite.col(1)=2*zVals;
    }
    for (int i=2;i<=n;i++)
    {
        resultat_hermite.col(i)=2*zVals*resultat_hermite.col(i-1)-2*(n-1)*resultat_hermite.col(i-2);
    }    
    
}

arma::vec Poly::hermite(int n)
{
    return resultat_hermite.col(n);
}
    