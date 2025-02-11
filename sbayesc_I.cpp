#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include <stdexcept>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

std::random_device rd;
std::mt19937 rng(123);

void readBinFullLD(const std::string& binFilePath, unsigned &numSNP, MatrixXf &B) {
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


void readBinTxtFile(const std::string& binFilePath, unsigned& numSNP, Eigen::MatrixXf& B) {
    std::ifstream file(binFilePath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << binFilePath << std::endl;
        return;
    }

    std::string line;
    numSNP = 1000;  // Expected matrix size
    B = Eigen::MatrixXf(numSNP, numSNP);

    unsigned row = 0;
    while (std::getline(file, line) && row < numSNP) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::vector<double> values;
        double value;

        while (ss >> value) {
            values.push_back(value);
        }

        for (unsigned col = 0; col < std::min(values.size(), static_cast<size_t>(numSNP)); col++) {
            B(row, col) = values[col];
        }

        row++;
    }

    file.close();
}


void readSummary(const std::string& summaryFilePath, Eigen::VectorXf& b, Eigen::VectorXf& se, Eigen::VectorXf& n, unsigned numSNP) {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    std::ifstream file(summaryFilePath);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Unable to open file " + summaryFilePath);
    }

    std::string line;
    std::vector<float> b_values, se_values, n_values;

    // Read header and determine column indices
    if (!std::getline(file, line)) {
        throw std::runtime_error("Error: Failed to read header from file " + summaryFilePath);
    }

    std::istringstream header_stream(line);
    std::vector<std::string> headers;
    std::string header;

    while (header_stream >> header) {  // Tokenizing header (space or tab)
        headers.push_back(header);
    }

    // Find indices of b, se, and n
    int b_idx = -1, se_idx = -1, n_idx = -1;
    for (size_t i = 0; i < headers.size(); i++) {
        if (headers[i] == "b") b_idx = i;
        else if (headers[i] == "se") se_idx = i;
        else if (headers[i] == "n") n_idx = i;
    }

    // Ensure required columns exist
    if (b_idx == -1 || se_idx == -1) {
        throw std::runtime_error("Error: Required column(s) missing in the header!");
    }

    // Read and parse each line
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;

        while (iss >> token) {  // Tokenizing line (space or tab)
            tokens.push_back(token);
        }

        if (tokens.size() <= std::max(b_idx, se_idx)) {
            throw std::runtime_error("Error: Line has insufficient columns: " + line);
        }

        // Parse values
        float b_value = std::stof(tokens[b_idx]);
        float se_value = std::stof(tokens[se_idx]);
        float n_value = (n_idx != -1 && tokens.size() > n_idx) ? std::stof(tokens[n_idx]) : static_cast<float>(numSNP);

        b_values.push_back(b_value);
        se_values.push_back(se_value);
        n_values.push_back(n_value);
    }

    // Convert vectors to Eigen VectorXf
    b = Eigen::Map<Eigen::VectorXf>(b_values.data(), b_values.size());
    se = Eigen::Map<Eigen::VectorXf>(se_values.data(), se_values.size());
    n = Eigen::Map<Eigen::VectorXf>(n_values.data(), n_values.size());

    std::cout << "Summary statistics loaded: " << b.size() << " entries." << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time to read GWAS summary statistics file: " << elapsed.count() << " seconds" << std::endl;
}


void saveMatrixToBinary(const std::string& filename, const Eigen::MatrixXf& matrix) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    int rows = matrix.rows();
    int cols = matrix.cols();
    file.write(reinterpret_cast<char*>(&rows), sizeof(int));
    file.write(reinterpret_cast<char*>(&cols), sizeof(int));
    file.write(reinterpret_cast<const char*>(matrix.data()), rows * cols * sizeof(float));

    file.close();
}

float sample_normal(float mean, float stddev) {
    std::normal_distribution<float> norm_dist(mean, stddev);
    return norm_dist(rng);
}

float sample_chisq(float dof) {
    std::chi_squared_distribution<float> chi_dist(dof);
    return chi_dist(rng);
}

float sample_beta(float alpha, float beta) {
    std::gamma_distribution<float> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<float> gamma_beta(beta, 1.0);
    float x = gamma_alpha(rng);
    float y = gamma_beta(rng);
    return x / (x + y);
}

float sample_bernoulli(float p) {
    std::discrete_distribution<int> dist({ 1 - p, p });
    return dist(rng);
}


int main() {
    unsigned numSNP;
    VectorXf b, se, n;
    MatrixXf LD;

    //std::string binFilePath = "1000G_eur_chr22.ldm.full.bin";
    //std::string phenoFilePath = "sim_1.ma";
    std::string binFilePath = "ldm_data1.ma";
    std::string phenoFilePath = "GWASss.ma";

    //readBinFullLD(binFilePath, numSNP, LD);
    readBinTxtFile(binFilePath, numSNP, LD);
    readSummary(phenoFilePath, b, se, n ,numSNP);

    unsigned n_iter = 10000;
    float pi_init = 0.1;
    float hsq_init = 0.5;

    VectorXf pi = VectorXf::Zero(n_iter+1);               
    VectorXf hsq = VectorXf::Zero(n_iter+1);          
    pi(0) = pi_init; 
    hsq(0) = hsq_init;

    //VectorXf scale = (1.0 / (numSNP * se.array().square())).sqrt(); 
    //VectorXf bhat = b.array() * scale.array();
    VectorXf bhat = b;

    float vary = 1.0;
    float varg = hsq(0);
    float vare = vary;
    VectorXf sigmaSq = VectorXf::Zero(n_iter+1);
    sigmaSq(0) = varg/(numSNP * pi(0));

    float nub = 4.0f, nue = 4.0f;
    float scaleb = (nub - 2) / nub * sigmaSq(0);
    float scalee = (nue - 2) / nue * vare;

    VectorXf beta = VectorXf::Zero(numSNP);
    MatrixXf beta_mcmc = MatrixXf::Zero(n_iter+1, numSNP);
    VectorXf bhatcorr = bhat;

    MatrixXf keptIter = MatrixXf::Zero(n_iter, 4);
    float invSigmaSq,ssq;
    VectorXd nnz = VectorXd::Zero(n_iter);
    float beta_old, rhs, invLhs, uhat;
    float logDelta_active, logDelta_inactive, pi_current, delta;
    float sigma_beta, sigma_epsilon;
    
    std::cout << std::left << std::setw(10) << "Iteration"
         << std::left << std::setw(10) << "pi" 
         << std::setw(10) << "nnz" 
         << std::setw(15) << "sigma_beta" 
         << std::setw(10) << "hsq" << std::endl;
    
    for (int i = 1; i < n_iter; i++){
        invSigmaSq = 1.0 / sigmaSq(i-1);
        ssq = 0;
        VectorXf numSnpDist_current = VectorXf::Zero(2);

        for (int j = 0; j < numSNP; j++) {
            beta_old = beta(j);
            //rhs = (bhatcorr(j) + beta_old) / (vare / numSNP);
            rhs = bhatcorr(j) + beta_old;
            invLhs = 1.0f / (1.0f / (vare / numSNP) + invSigmaSq);
            uhat = invLhs * rhs;

            logDelta_active = 0.5f*(log(invLhs) - log(sigmaSq(i-1)) + uhat*rhs) + log(pi(i-1));
            logDelta_inactive = log(1.0f-pi(i-1)); 
            pi_current  = 1.0 / (1.0 + exp(logDelta_inactive - logDelta_active));

            delta = sample_bernoulli(pi_current);
            if (delta > 0){
                numSnpDist_current(1) += 1;
                beta(j) = sample_normal(uhat, sqrt(invLhs));
                 
                bhatcorr = bhatcorr.array() + LD.col(j).array()*(beta_old - beta(j));
                ssq = ssq + (beta(j)*beta(j));
                nnz(i-1) = nnz(i - 1) + 1;
            }else{
                numSnpDist_current(0) += 1;
                bhatcorr = bhatcorr.array() + LD.col(j).array()*beta_old;
                beta(j) = 0.0;
            }

            vare = (bhatcorr - bhat).squaredNorm() / numSNP;

        }
        
        VectorXf beta_mcmc_sample = beta.array();// / scale.array();
        beta_mcmc.row(i) = beta_mcmc_sample.transpose();

        pi(i) = sample_beta(numSnpDist_current(1), numSnpDist_current(0));

        //sigma_beta = sample_chisq(nnz(i-1) + nub);
        //sigmaSq(i) = (ssq + nub*scaleb)/sigma_beta;

        sigma_beta = sample_chisq(nnz(i-1) + nub);
        sigmaSq(i) = (ssq / numSNP + nub * scaleb) / sigma_beta;

        varg = beta.dot(bhat - bhatcorr);
        hsq(i) = varg /vary;

        if (i % 500 == 0){
        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::left << std::setw(10) << i
            << std::left << std::setw(10) << pi(i) 
            << std::setw(10) << int(nnz(i-1)) 
            << std::setw(15) << sigmaSq(i) 
            << std::setw(10) << hsq(i) << std::endl;
        }

        keptIter.row(i-1) << pi(i), int(nnz(i-1)), sigmaSq(i), hsq(i);

    }

    saveMatrixToBinary("ldm_data1_I_result.bin", beta_mcmc);

    return 0;
};
