#include "../include/poly.h"
#include "../include/basis.h"
#include "../include/solver.h"

#include <armadillo>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

// ---------------------------------------------------------------------
// Convertit un arma::cube en string binaire au format DF3 pour POV-Ray.
// ---------------------------------------------------------------------
std::string cubeToDf3(const arma::cube &m)
{
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);

    int nx = static_cast<int>(m.n_rows);
    int ny = static_cast<int>(m.n_cols);
    int nz = static_cast<int>(m.n_slices);

    // Header DF3 : 3 entiers 16 bits big-endian : nx, ny, nz
    ss.put(static_cast<char>(nx >> 8));
    ss.put(static_cast<char>(nx & 0xff));
    ss.put(static_cast<char>(ny >> 8));
    ss.put(static_cast<char>(ny & 0xff));
    ss.put(static_cast<char>(nz >> 8));
    ss.put(static_cast<char>(nz & 0xff));

    double theMax = m.max();

    // Sécurité pour éviter la division par zéro
    if (theMax <= 1e-16) theMax = 1.0;

    for (arma::uword k = 0; k < m.n_slices; ++k) {       // z
        for (arma::uword j = 0; j < m.n_cols;   ++j) {   // y
            for (arma::uword i = 0; i < m.n_rows; ++i) { // x

                double val = m(i, j, k);

                // Clamp à 0 pour éliminer le "brouillard" (valeurs -1e-15 etc.)
                if (val < 0) val = 0.0;

                // Normalisation [0, 255]
                double norm = val / theMax;
                if (norm > 1.0) norm = 1.0;

                unsigned char v = static_cast<unsigned char>(255.0 * norm);
                ss.put(static_cast<char>(v));
            }
        }
    }

    return ss.str();
}

// ---------------------------------------------------------------------
// Construit un cube 3D rho(x,y,z) à partir d'une densité 2D rho(r_perp,z).
// ---------------------------------------------------------------------
arma::cube buildDensityCubeFrom2D(const arma::mat &density2D,
                                  const arma::vec &rVals,
                                  const arma::vec &zVals)
{
    const int nx = 100;
    const int ny = 100;
    const int nz = 100;

    arma::cube cube(nx, ny, nz, arma::fill::zeros);

    // Domaines physiques : X[-10,10], Y[-10,10], Z[-20,20]
    arma::vec xVals = arma::linspace(-10.0, 10.0, nx);
    arma::vec yVals = arma::linspace(-10.0, 10.0, ny);
    arma::vec zAxisCube = arma::linspace(-20.0, 20.0, nz);

    const double rMax = rVals.max();
    const double zMin = zVals.min();
    const double zMax = zVals.max();

    const int nbR = static_cast<int>(rVals.n_elem);
    const int nbZ = static_cast<int>(zVals.n_elem);

    for (int k = 0; k < nz; ++k) {
        double z = zAxisCube(k);

        double tZ = (z - zMin) / (zMax - zMin);
        int idxZ = static_cast<int>(std::round(tZ * (nbZ - 1)));

        if (idxZ < 0) idxZ = 0;
        if (idxZ >= nbZ) idxZ = nbZ - 1;

        for (int j = 0; j < ny; ++j) {
            double y = yVals(j);
            for (int i = 0; i < nx; ++i) {
                double x = xVals(i);

                double r = std::sqrt(x * x + y * y);
                double value = 0.0;

                if (r <= rMax && nbR > 1) {
                    double tR = r / rMax;
                    int idxR = static_cast<int>(std::round(tR * (nbR - 1)));

                    if (idxR < 0) idxR = 0;
                    if (idxR >= nbR) idxR = nbR - 1;

                    value = density2D(idxZ, idxR);
                }

                cube(i, j, k) = value;
            }
        }
    }

    return cube;
}

int main(int argc, char* argv[])
{
    // -----------------------------------------------------------------
    // 1. SETUP & LOADING
    // -----------------------------------------------------------------

    // Paramètres du sujet (br, bz, N, Q)
    Basis basis(1.935801664793151, 2.829683956491218, 14, 1.3);

    // Fichier de matrice de densité (rho)
    std::string rhoPath = "data/rho.arma";
    if (argc > 1) {
        rhoPath = argv[1];
    }

    std::cout << "Loading density matrix from: " << rhoPath << std::endl;

    arma::mat rho_mat;
    if (!rho_mat.load(rhoPath, arma::arma_ascii)) {
        std::cerr << "Error: Could not load " << rhoPath << std::endl;
        std::cerr << "Usage: ./bin/main [path_to_rho.arma]" << std::endl;
        return 1;
    }

    std::cout << "Loaded rho matrix: " << rho_mat.n_rows
              << " x " << rho_mat.n_cols << std::endl;

    // Grilles (R, Z)
    int nbR = 100;   // r : 0 -> 10
    int nbZ = 200;   // z : -20 -> 20
    arma::vec rVals = arma::linspace(0.0, 10.0, nbR);
    arma::vec zVals = arma::linspace(-20.0, 20.0, nbZ);

    arma::wall_clock timer;

    // -----------------------------------------------------------------
    // 2. OPTIMIZED ALGORITHM
    // -----------------------------------------------------------------

    std::cout << "\nStarting OPTIMIZED Algorithm..." << std::endl;

    Solver solver(basis);
    solver.precomputeBasis(rVals, zVals);

    timer.tic();
    arma::mat densityOptimized = solver.calcDensityOptimized(rho_mat);
    double timeOptimized = timer.toc();

    std::cout << ">>> Optimized Time: " << timeOptimized
              << " seconds" << std::endl;

    // -----------------------------------------------------------------
    // 3. NAIVE ALGORITHM (Pour vérification + export)
    // -----------------------------------------------------------------

    std::cout << "\nStarting NAIVE Algorithm..." << std::endl;

    arma::mat densityNaive = arma::zeros(nbZ, nbR);
    timer.tic();

    int idx_a = 0;
    for (int m_a = 0; m_a < basis.mMax; ++m_a) {
        for (int n_a = 0; n_a < basis.nMax(m_a); ++n_a) {
            for (int nz_a = 0; nz_a < basis.n_zMax(m_a, n_a); ++nz_a) {

                arma::vec R_a = basis.rPart(rVals, m_a, n_a);
                arma::vec Z_a = basis.zPart(zVals, nz_a);
                arma::mat Psi_a = Z_a * R_a.t();

                int idx_b = 0;
                for (int m_b = 0; m_b < basis.mMax; ++m_b) {
                    for (int n_b = 0; n_b < basis.nMax(m_b); ++n_b) {
                        for (int nz_b = 0; nz_b < basis.n_zMax(m_b, n_b); ++nz_b) {

                            double rho_val = rho_mat(idx_a, idx_b);

                            if (std::fabs(rho_val) > 1e-12) {
                                arma::vec R_b = basis.rPart(rVals, m_b, n_b);
                                arma::vec Z_b = basis.zPart(zVals, nz_b);
                                arma::mat Psi_b = Z_b * R_b.t();
                                densityNaive += rho_val * (Psi_a % Psi_b);
                            }
                            ++idx_b;
                        }
                    }
                }
                ++idx_a;
            }
        }
        std::cout << "." << std::flush;
    }

    double timeNaive = timer.toc();
    std::cout << "\n>>> Naive Time: " << timeNaive << " seconds" << std::endl;

    // -----------------------------------------------------------------
    // 4. RESULTATS & SAUVEGARDE (OPT + NAIVE)
    // -----------------------------------------------------------------

    std::cout << "\n--- RESULTS ---" << std::endl;
    std::cout << "Speedup Factor: " << timeNaive / timeOptimized << "x" << std::endl;

    double diffNorm = arma::norm(densityNaive - densityOptimized);
    std::cout << "Difference (Norm): " << diffNorm << std::endl;

    const char* outDir = "output";
    struct stat st;
    if (stat(outDir, &st) != 0) { mkdir(outDir, 0755); }

    std::string prefix(outDir);
    prefix += "/";

    // Grilles
    rVals.save(prefix + "grid_r.arma", arma::arma_ascii);
    zVals.save(prefix + "grid_z.arma", arma::arma_ascii);

    // Densités 2D (arma)
    densityOptimized.save(prefix + "density_optimized.arma", arma::arma_ascii);
    densityNaive.save(prefix + "density_naive.arma", arma::arma_ascii);

    // DF3 Optimized
    std::cout << "Building 3D Cube OPTIMIZED (100x100x100)..." << std::endl;
    arma::cube densityCubeOpt = buildDensityCubeFrom2D(densityOptimized, rVals, zVals);
    std::string df3DataOpt = cubeToDf3(densityCubeOpt);

    std::string df3PathOpt = prefix + "density_optimized.df3";
    std::ofstream df3FileOpt(df3PathOpt, std::ios::binary);
    if (!df3FileOpt) {
        std::cerr << "Error: cannot open " << df3PathOpt << " for writing." << std::endl;
    } else {
        df3FileOpt.write(df3DataOpt.data(), static_cast<std::streamsize>(df3DataOpt.size()));
        df3FileOpt.close();
        std::cout << "DF3 volume saved to " << df3PathOpt << std::endl;
    }

    // DF3 Naive
    std::cout << "Building 3D Cube NAIVE (100x100x100)..." << std::endl;
    arma::cube densityCubeNaive = buildDensityCubeFrom2D(densityNaive, rVals, zVals);
    std::string df3DataNaive = cubeToDf3(densityCubeNaive);

    std::string df3PathNaive = prefix + "density_naive.df3";
    std::ofstream df3FileNaive(df3PathNaive, std::ios::binary);
    if (!df3FileNaive) {
        std::cerr << "Error: cannot open " << df3PathNaive << " for writing." << std::endl;
    } else {
        df3FileNaive.write(df3DataNaive.data(), static_cast<std::streamsize>(df3DataNaive.size()));
        df3FileNaive.close();
        std::cout << "DF3 volume saved to " << df3PathNaive << std::endl;
    }

    return 0;
}
