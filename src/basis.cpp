#include "../include/basis.h"
#include "../include/poly.h"

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
    return (((N+2.0)*pow(Q,2.0/3.0)+(1.0/2.0)-i*Q));
}


void Basis::setmMax(){
    mMax = std::floor(((N+2.0)*pow(Q,2.0/3.0)-(1.0/2.0))/Q) ;
    
}

void Basis::setnMax(){
    nMax=arma::zeros<arma::ivec>(mMax);
    for (int m=0;m<mMax;m++)
    {
        nMax(m) = floor((mMax-m+1.0)/2.0);
    }
    
}

void Basis::setnZMax()
{
    double nMax_0_float = 0.5 * (mMax - 1) + 1.0;
    int max_n_size = static_cast<int>(floor(nMax_0_float));
    n_zMax=arma::zeros<arma::imat>(mMax,max_n_size);
    for (int m=0; m<mMax; m++)
    {
        for (int n=0; n<nMax(m);n++)
        {
            int i = m+2*n+1;
            n_zMax(m,n)=floor(Basis::getnu(i,N,Q));
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

arma::vec Basis::rPart(arma::vec r, int m, int n)
{
    r00=arma::ones(r.size());
    Poly poly;
    poly.calcLaguerre(fabs(m)+1,n+1,r%r/(br*br));
    double factor =(1.0/(br*sqrt(pi)))*sqrt(tgamma(n + 1.0)/tgamma((n+fabs(m))+1.0));
    arma::vec exp_term=arma::exp(-pow(r,2.0)/(2.0*pow(br,2.0)));
    arma::vec pow_term=arma::pow(r/br,fabs(m));
    arma::vec laguerre_term=poly.laguerre(fabs(m),n);
    r00 = factor * exp_term % pow_term % laguerre_term;
    return r00;

}

