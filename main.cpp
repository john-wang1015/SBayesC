#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "data.hpp"

int main() {
    try {
        BinaryFileReader reader("data/1000G_eur_chr22.ldm.full.bin");

        // Read entire file as raw data
        auto rawData = reader.readAll();
        std::cout << "Raw data size: " << rawData.size() << " bytes\n";

        // Read file as integers
        auto intData = reader.readAsType<int>();
        std::cout << "Read " << intData.size() << " integers from the file:\n" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}