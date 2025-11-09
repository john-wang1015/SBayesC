/*
    This is the data.cpp file used to read the summary statistics data
*/

#include "data.hpp"
#include <chrono>



void Data::readBinInfo(const std::string& binInfoFilePath){
    /*
        Read .bin.info data for LD matrix
    */

}

void Data::readBinFullLD(const std::string& binFilePath) {
    /*
        read the .bin file for full LD matrix
        NB, full LD matrix is square matrix, i.e., size of numSNP * numSNP
    */
    auto start = std::chrono::high_resolution_clock::now(); // Start timing
    std::ifstream file(binFilePath, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Error: Unable to open file " + binFilePath);
    }

    std::streamsize fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    if (fileSize % sizeof(float) != 0) {
        throw std::runtime_error("Error: File size is not a multiple of float size. File: " + binFilePath);
    }

    std::vector<char> buffer(fileSize);
    if (!file.read(buffer.data(), fileSize)) {
        throw std::runtime_error("Error: Failed to read file " + binFilePath);
    }

    numSNP = std::sqrt(fileSize / sizeof(float));
    double dimension = std::sqrt(fileSize / sizeof(float));
    int dim = static_cast<int>(dimension);

    if (std::floor(dimension) != dimension) {
        throw std::runtime_error("Error: The LD matrix is not a square matrix. File: " + binFilePath);
    }

    const float* floatData = reinterpret_cast<const float*>(buffer.data());
    Eigen::Map<const Eigen::MatrixXf> mappedMatrix(floatData, dim, dim);

    B = mappedMatrix;

    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for read .bin file: " << elapsed.count() << " seconds" << std::endl;
}


void Data::readSummary(const std::string &summaryFilePath) {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing
    std::ifstream file(summaryFilePath);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Unable to open file " + summaryFilePath);
    }

    std::string line;
    std::vector<float> b_values;
    std::vector<float> se_values;
    std::vector<float> n_values;

    // Read and ignore the header line
    if (!std::getline(file, line)) {
        throw std::runtime_error("Error: Failed to read header from file " + summaryFilePath);
    }

    // Read each line and parse values
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;

        // Tokenize the line using space or tab as the delimiter
        while (std::getline(iss, token, ' ')) {
            if (!token.empty()) {  // Skip empty tokens
                tokens.push_back(token);
            }
        }

        // Ensure there are enough columns
        if (tokens.size() < 8) {
            throw std::runtime_error("Error: Insufficient columns in line: " + line);
        }

        // Parse b, se, and n
        float b = std::stof(tokens[4]);      // Column 5: b
        float se = std::stof(tokens[5]);     // Column 6: se
        float n = std::stod(tokens[7]);     // Column 8: n

        b_values.push_back(b);
        se_values.push_back(se);
        n_values.push_back(n);
    }

    // Convert vectors to Eigen types
    bhat = Eigen::Map<Eigen::VectorXf>(b_values.data(), b_values.size());
    se = Eigen::Map<Eigen::VectorXf>(se_values.data(), se_values.size());
    n = Eigen::Map<Eigen::VectorXf>(n_values.data(), n_values.size());


    std::cout << "Summary statistics loaded: " << bhat.size() << " entries." << std::endl;
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for read GWAS summary statistics file: " << elapsed.count() << " seconds" << std::endl;
}



