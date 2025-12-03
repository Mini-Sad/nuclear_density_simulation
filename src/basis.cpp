#include "../include/basis.h"

Basis::Basis(double _br, double _bz, int _N, double _Q) : br(_br), bz(_bz), N(_N), Q(_Q) {
      
    setmMax();
    setnMax();
    setnZMax();
    getmMax();
    getnMax();
    getnZMax();

}
double Basis::getnu(int i, int N, double Q)
{
    return (((N+2.0)*std::pow(Q,2.0/3.0)+(1.0/2.0)-i*Q));
}


void Basis::setmMax(){
    mMax = std::floor(((N+2.0)*std::pow(Q,2.0/3.0)-(1.0/2.0))/Q) ;
    
}

void Basis::setnMax(){
    nMax=arma::zeros<arma::ivec>(mMax);
    for (int m=0;m<mMax;m++)
    {
        nMax(m) = std::floor((mMax-m+1.0)/2.0);
    }
    
}

void Basis::setnZMax()
{
    n_zMax=arma::zeros<arma::imat>(mMax,nMax.size());
    for (int m=0; m<mMax; m++)
    {
        for (int n=0; n<nMax(m);n++)
        {
            int i = m+2*n+1;
            n_zMax(m,n)=std::floor(Basis::getnu(i,N,Q));
        }
    }
}

int Basis::getmMax()
{
    return mMax;
}

arma::ivec Basis::getnMax()
{
    return nMax;
}

arma::imat Basis::getnZMax()
{
    return n_zMax;
}

