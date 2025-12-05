// ... [Previous calculation code] ...
    
    // --- VISUALIZATION STEP ---
    #include "writer.hpp"

    std::cout << "Preparing 3D Cube for Visualization..." << std::endl;

    // 1. Define 3D Grid Parameters (Mandatory)
    int nCubeX = 32;
    int nCubeY = 32;
    int nCubeZ = 64;
    arma::cube vol(nCubeX, nCubeY, nCubeZ);

    // Physical limits for the cube
    double minX = -10.0, maxX = 10.0;
    double minY = -10.0, maxY = 10.0;
    double minZ = -20.0, maxZ = 20.0;

    // 2. Loop over every voxel in the 3D cube
    for (int k = 0; k < nCubeZ; ++k) {
        // Map index k to physical Z
        double z_phys = minZ + k * (maxZ - minZ) / (nCubeZ - 1);
        
        // Find closest index in our calculated 2D density matrix
        // We assume zVals was created with linspace(-20, 20, nbZ)
        // Simple Nearest Neighbor mapping:
        int z_idx = (int)std::round( (z_phys - (-20.0)) / (40.0) * (nbZ - 1) );
        // Bounds check
        if (z_idx < 0) z_idx = 0;
        if (z_idx >= nbZ) z_idx = nbZ - 1;

        for (int j = 0; j < nCubeY; ++j) {
            double y_phys = minY + j * (maxY - minY) / (nCubeY - 1);

            for (int i = 0; i < nCubeX; ++i) {
                double x_phys = minX + i * (maxX - minX) / (nCubeX - 1);

                // 3. Calculate Radius r
                double r_phys = std::sqrt(x_phys*x_phys + y_phys*y_phys);

                // 4. Map r to our 2D density matrix
                // We assume rVals was linspace(0, 10, nbR)
                int r_idx = (int)std::round( (r_phys - 0.0) / (10.0) * (nbR - 1) );

                double val = 0.0;
                // If r is within our calculated range, grab the value
                if (r_idx >= 0 && r_idx < nbR) {
                    val = density(z_idx, r_idx);
                } else {
                    // Outside calculation domain (r > 10)
                    val = 0.0;
                }

                vol(i, j, k) = val;
            }
        }
    }

    // 3. Export
    std::string df3Data = Writer::cubeToDf3(vol);
    Writer::saveString("density.df3", df3Data);
