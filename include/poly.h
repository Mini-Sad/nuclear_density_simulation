#ifndef POLY_H
#define POLY_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace std;


class Poly{

    private:
    arma::vec resultat_hermite;
    arma::vec resultat_laguerre;

    public:
    /**
     * Constructeur
     */
    Poly();
    /**
     * @brief méthide calculant les polynômes d'Hermite pour n dans {0,n-1)
     * @param int : n 
     * @param vec : zVals
     * @return nothing
     */
    void calcHermite(int, arma::vec);

    /**
     * @brief fonction retournant le polynôme d'hermite pour n suite au calcul
     * @param int :n
     * @return resultat du polynôme d'hermite pour n
     */

    arma::vec hermite(int);
    

    /**
     * @brief méthode calculant le polynôme de Laguerre sur des valeurs dépendant de m et n
     * @param int : m
     * @param int : n
     * @param vec : zVals
     * @return  rien
     */
    void calcLaguerre(int, int, arma::vec);

    /**
     * @brief fonction prenant le calcul du polynôme de Lagrange pour un n et m donné
     * @param int : m
     * @param int : n
     * @return le résulat du polynôme de Lagrange pour m et n donné
     */
    arma::vec laguerre(int, int);

};

#endif // POLY_H