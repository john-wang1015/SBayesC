/*
    This is the data.cpp file used to read the summary statistics data
*/

#include "data.hpp"
#include <chrono>



void readFile::readBinInfo(const std::string& binInfoFilePath){
    /*
        Read .bin.info data for LD matrix
    */

}

void readFile::readBinFullLD(const std::string& binFilePath) {
    /*
        read the .bin file for full LD matrix
        NB, full LD matrix is square matrix, i.e., size of numSNP * numSNP
    */
    auto start = std::chrono::high_resolution_clock::now(); // Start timing
    std::ifstream file(binFilePath, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Error: Unable to open file " + binFilePath);
    }

    // Get file size
    std::streamsize fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    if (fileSize % sizeof(float) != 0) {
        throw std::runtime_error("Error: File size is not a multiple of float size. File: " + binFilePath);
    }

    // Read the binary data into a buffer
    std::vector<char> buffer(fileSize);
    if (!file.read(buffer.data(), fileSize)) {
        throw std::runtime_error("Error: Failed to read file " + binFilePath);
    }

    // Calculate matrix dimensions
    double dimension = std::sqrt(fileSize / sizeof(float));
    int dim = static_cast<int>(dimension);

    if (std::floor(dimension) != dimension) {
        throw std::runtime_error("Error: The LD matrix is not a square matrix. File: " + binFilePath);
    }

    // Map the buffer directly to an Eigen matrix
    const float* floatData = reinterpret_cast<const float*>(buffer.data());
    Eigen::Map<const Eigen::MatrixXf> mappedMatrix(floatData, dim, dim);

    // Copy the mapped matrix to B
    B = mappedMatrix;

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for read .bin file: " << elapsed.count() << " seconds" << std::endl;
}


void readFile::readSummary(const std::string& summaryFilePath) {
    /*
        Read the summary statistics from the .ma file
    */
    std::ifstream file(summaryFilePath);
    if (!file) {
        throw std::runtime_error("Error: Unable to open file " + summaryFilePath);
    }

    // Header for the summary statistics
    std::vector<std::string> header = {
        "Chrom", "ID", "GenPos", "PhysPos", "A1", "A2", "A1Freq", "Index",
        "WindStart", "WindEnd", "WindSize", "WindWidth", "N", "SamplVar", "LDsum"
    };


   


}

