#include <iostream>
#include <armadillo>
#include "basis.hpp"
#include "solver.hpp"
#include "writer.hpp" 

int main() {
    // 1. Setup the Physics Constants
    // From Project: br=1.9358..., bz=2.8296..., N=14, Q=1.3
    Basis basis(1.935801664793151, 2.829683956491218, 14, 1.3);

    // 2. Load the Density Matrix
    arma::mat rho_mat;
    // Ensure "rho.arma" is in your build/execution directory
    bool status = rho_mat.load("rho.arma", arma::arma_ascii);
    if (!status) {
        std::cerr << "Error: Could not load rho.arma" << std::endl;
        return 1;
    }
    std::cout << "Loaded rho matrix: " << rho_mat.n_rows << "x" << rho_mat.n_cols << std::endl;

    // 3. Define the Grid (r and z axes)
    int nbR = 32;
    int nbZ = 64;
    arma::vec rVals = arma::linspace(0.0, 10.0, nbR);   // 0 to 10 fm
    arma::vec zVals = arma::linspace(-20.0, 20.0, nbZ); // -20 to 20 fm

    // 4. Use the Solver (Optimized)
    Solver solver(basis);
    
    // Pre-calculate everything (this handles the heavy lifting)
    solver.precomputeBasis(rVals, zVals);
    
    // Timer start
    arma::wall_clock timer;
    timer.tic();
    
    // Run the Optimized Calculation
    arma::mat density = solver.calcDensityOptimized(rho_mat);
    
    // Timer stop
    double duration = timer.toc();
    std::cout << "Optimized Time: " << duration << " seconds" << std::endl;
    
    // 5. Save Results
    density.save("density.mat", arma::arma_ascii);

    // 6. Visualization Export (Optional but recommended)
    // Map to 3D Cube for POV-Ray
    std::cout << "Generating 3D Visualization..." << std::endl;
    
    int nCubeX = 32, nCubeY = 32, nCubeZ = 64;
    arma::cube vol(nCubeX, nCubeY, nCubeZ);
    double minX = -10.0, maxX = 10.0;
    double minY = -10.0, maxY = 10.0;
    double minZ = -20.0, maxZ = 20.0;

    for (int k = 0; k < nCubeZ; ++k) {
        double z_phys = minZ + k * (maxZ - minZ) / (nCubeZ - 1);
        int z_idx = (int)std::round( (z_phys - (-20.0)) / (40.0) * (nbZ - 1) );
        if (z_idx < 0) z_idx = 0; 
	if (z_idx >= nbZ) z_idx = nbZ - 1;

        for (int j = 0; j < nCubeY; ++j) {
            double y_phys = minY + j * (maxY - minY) / (nCubeY - 1);
            for (int i = 0; i < nCubeX; ++i) {
                double x_phys = minX + i * (maxX - minX) / (nCubeX - 1);
                double r_phys = std::sqrt(x_phys*x_phys + y_phys*y_phys);
                int r_idx = (int)std::round( (r_phys - 0.0) / (10.0) * (nbR - 1) );

                if (r_idx >= 0 && r_idx < nbR) {
                    vol(i, j, k) = density(z_idx, r_idx);
                } else {
                    vol(i, j, k) = 0.0;
                }
            }
        }
    }
    
    std::string df3Data = Writer::cubeToDf3(vol);
    Writer::saveString("density.df3", df3Data);

    return 0;
}
