#include "../include/writer.hpp"
#include <sstream>
#include <cmath>
#include <fstream>
#include <iostream>

std::string Writer::cubeToDf3(const arma::cube &m) {
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    
    // Header: Size in X, Y, Z (2 bytes each, Big Endian)
    int nx = m.n_rows;
    int ny = m.n_cols;
    int nz = m.n_slices;

    ss.put(nx >> 8); ss.put(nx & 0xff);
    ss.put(ny >> 8); ss.put(ny & 0xff);
    ss.put(nz >> 8); ss.put(nz & 0xff);

    double theMin = 0.0;
    double theMax = m.max();

    // Loop Z, Y, X to write voxel data
    // Note: POV-Ray expects integer values 0-255
    for (uint k = 0; k < m.n_slices; k++) {
        for (uint j = 0; j < m.n_cols; j++) {
            for (uint i = 0; i < m.n_rows; i++) {
                double val = m(i, j, k);
                // Clamp and scale
                uint v = 255 * (std::fabs(val) - theMin) / (theMax - theMin);
                ss.put(v);
            }
        }
    }
    return ss.str();
}

std::string Writer::cubeToRaw(const arma::cube &m) {
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    double theMin = 0.0;
    double theMax = m.max();

    for (uint k = 0; k < m.n_slices; k++) {
        for (uint j = 0; j < m.n_cols; j++) {
            for (uint i = 0; i < m.n_rows; i++) {
                 uint v = 255 * (std::fabs(m(i, j, k)) - theMin) / (theMax - theMin);
                ss.put(v);
            }
        }
    }
    return ss.str();
}

void Writer::saveString(const std::string& filename, const std::string& data) {
    std::ofstream out(filename.c_str(), std::ios::binary);
    if (out.is_open()) {
        out << data;
        out.close();
        std::cout << "Saved " << filename << std::endl;
    } else {
        std::cerr << "Error: Could not save " << filename << std::endl;
    }
}
