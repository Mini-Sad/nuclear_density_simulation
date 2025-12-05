#ifndef WRITER_HPP
#define WRITER_HPP

#include <armadillo>
#include <string>

class Writer {
public:
    // Convert a density cube to POV-Ray df3 format string
    static std::string cubeToDf3(const arma::cube &m);
    
    // Convert a density cube to Blender raw format string
    static std::string cubeToRaw(const arma::cube &m);

    // Helper to save string to file
    static void saveString(const std::string& filename, const std::string& data);
};

#endif
