#include "../include/basis.h"

Basis::Basis(double _br, double _bz, int _N, double _Q) : br(_br), bz(_bz), N(_N), Q(_Q) {

    calcV();
    getmax();

}


void Basis::calcV(){
    vVals = arma::ones(N);
    for (int i=0; i<N;i++)
    {
        vVals(i)=(N+2)*pow(Q,2.0/3.0)+(1.0/2.0)-i*Q;
    }

}

int Basis::getmax(){
    for(int i = 0 ; i < N ; i++)
    {
        if(vVals(i)>=1)
        {
            max=i;
        }
    }
    return max;
}

