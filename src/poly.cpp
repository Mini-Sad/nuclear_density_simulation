#include "../include/poly.h"

//Constructeur
Poly::Poly(){

}

void calcHermite(int n, arma::vec zVals){
    arma::vec resultat = arma::ones(zVals.size());
    if(n>0)
    {
        resultat.col(1)=2*zVals;
    }
    for (int i=2;i<=n;i++)
    {
        resultat.col(i)=2*zVals*resultat.col(i-1)-2*(n-1)*resultat.col(i-2);
    }    
    
}
/*
arma::vec hermite(int n)
{
    return resultat.col(n);
}
    */