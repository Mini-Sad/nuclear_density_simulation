#include "../include/poly.h"

//Constructeur
Poly::Poly(){
    
}

void Poly::calcHermite(int n, arma::vec zVals){
    resultat_hermite = arma::ones(zVals.size(),n+1);
    if(n>0)
    {
        resultat_hermite.col(1)=2*zVals;
    }
    for (int i=2;i<=n;i++)
    {
        resultat_hermite.col(i)=2*zVals%resultat_hermite.col(i-1)-2*(n-1)*resultat_hermite.col(i-2);
    }    
    
}

arma::vec Poly::hermite(int n)
{
    return resultat_hermite.col(n);
}

void Poly::calcLaguerre(int m, int n, arma::vec zVals)
{
    int k=zVals.n_elem;
    resultat_laguerre=arma::cube(m+1,n+1,k,arma::fill::zeros);
    for (int i=0;i<=m;i++)
    {   
        for(int K=0;K<k;K++)
        {
            resultat_laguerre(i,0,k)=1.0;
        }
        if(n>0)
        {
            for(int K=0;K<k;K++)
            {
                resultat_laguerre(i,1,K)=1.0+i-zVals(K);
            } 
        }
        for (int j=2;j<=n;j++)
        {
            for (int K=0;K<k;k++){
                resultat_laguerre(i,j,K)=(2.0+((i-1-zVals(K))/j))*resultat_laguerre(i,j-1,K)-(1+(i-1)/j)*resultat_laguerre(i,j-2,K);
            }
        }
    }
}

arma::vec Poly::laguerre(int m, int n)
{
    resultat_laguerre.tube(m,n);
}