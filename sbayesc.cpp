#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <Eigen/Dense>

using namespace Eigen;

std::random_device rd;
std::mt19937 rng(rd());

/*
    Error occur on readBinFullLD method, need to fix
*/

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

void readBinTxtFile(const std::string& binFilePath, unsigned &numSNP, MatrixXf &B){
    
}

void readSummary(const std::string &summaryFilePath, VectorXf &b, VectorXf &se) {
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
    b = Eigen::Map<Eigen::VectorXf>(b_values.data(), b_values.size());
    se = Eigen::Map<Eigen::VectorXf>(se_values.data(), se_values.size());
    //n = Eigen::Map<Eigen::VectorXf>(n_values.data(), n_values.size());

    std::cout << "Summary statistics loaded: " << b.size() << " entries." << std::endl;
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for read GWAS summary statistics file: " << elapsed.count() << " seconds" << std::endl;
}

float sample_normal(float mean, float stddev) {
    std::normal_distribution<float> norm_dist(mean, stddev);
    return norm_dist(rng);
}

float sample_uniform() {
    std::uniform_real_distribution<float> uniform_dist(0.0, 1.0);
    return uniform_dist(rng);
}

float sample_chi_squared(float dof) {
    std::chi_squared_distribution<float> chi_dist(dof);
    return chi_dist(rng);
}

float sample_scaled_inv_chi_squared(float dof, float scale) {
    std::chi_squared_distribution<float> chi_dist(dof);
    float x = chi_dist(rng);
    return (dof * scale) / x;
}

float sample_beta(float alpha, float beta) {
    std::gamma_distribution<float> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<float> gamma_beta(beta, 1.0);
    float x = gamma_alpha(rng);
    float y = gamma_beta(rng);
    return x / (x + y);
}

float sample_bernoulli(float p) {
    std::bernoulli_distribution bernoulli_dist(p);
    return bernoulli_dist(rng);
}

int main() {
    unsigned numSNP;
    VectorXf b, se;
    MatrixXf LD;

    std::string binFilePath = "1000G_eur_chr22.ldm.full.bin";
    std::string phenoFilePath = "sim_1.ma";

    readBinFullLD(binFilePath, numSNP, LD);
    readSummary(phenoFilePath, b, se);

    unsigned n_iter = 10000;
    float pi_init = 0.1;
    float hsq_init = 0.5;
    VectorXf gamma(2);
    gamma << 0, 1;

    VectorXf pi = VectorXf::Zero(n_iter+1);               // store all pi value
    VectorXf hsq = VectorXf::Zero(n_iter+1);              // store all heritability value
    pi(0) = pi_init; 
    hsq(0) = hsq_init;

    VectorXf scale = (1.0 / (numSNP * se.array().square())).sqrt(); 
    // scale SNP effect
    VectorXf bhat = b.array() * scale.array();  

    float vary = 1.0;
    float varg = hsq(0);
    float vare = vary;
    VectorXf sigmaSq = VectorXf::Zero(n_iter+1);
    sigmaSq(0) = varg / (numSNP * pi(0));

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
    
    std::cout << std::left << std::setw(10) << "pi" 
         << std::setw(10) << "nnz" 
         << std::setw(15) << "sigma_beta" 
         << std::setw(10) << "hsq" << std::endl;
    
    for (int i = 1; i < n_iter; i++){
        invSigmaSq = 1.0 / sigmaSq(i-1);
        nnz(i-1) = 0;
        ssq = 0;
        Vector2f numSnpDist_current;
        numSnpDist_current << 0.0, 0.0;

        for (int j = 0; j < numSNP; j++){
            beta_old = beta(j);
            rhs = (bhatcorr(j) + beta_old)/ (vare/numSNP);
            invLhs = 1.0/(1.0/(vare/numSNP) + invSigmaSq);
            uhat = invLhs * rhs;

            logDelta_active = 0.5*(log(invLhs) - log(sigmaSq(i-1)) + uhat*rhs) + log(pi(i-1));
            logDelta_inactive = log(1 - pi(i-1));
            pi_current = 1.0 / (1.0 + exp(logDelta_inactive - logDelta_active));

            delta = sample_bernoulli(pi_current);
            if (delta == 1){
                numSnpDist_current(1) += 1.0;

                beta(j) = sample_normal(uhat, sqrt(invLhs));
                bhatcorr = bhatcorr.array() + LD.col(j).array()*(beta_old - beta(j));
                ssq += (beta(j)*beta(j));
                nnz(i-1) += 1;
            }else{
                numSnpDist_current(0) += 1.0;
                bhatcorr = bhatcorr.array() + LD.col(j).array()*beta_old;
                beta(j) = 0.0;
            }
        }
        VectorXf beta_mcmc_sample = beta.array() / scale.array();
        beta_mcmc.row(i) = beta_mcmc_sample.transpose();

        pi(i) = sample_beta(numSnpDist_current(1) + 1, numSnpDist_current(0) + 1);

        sigma_beta = sample_scaled_inv_chi_squared(1, nnz(i-1)+nub);
        sigmaSq(i) = (ssq + nub*scaleb)/sigma_beta;

        varg = beta.dot(bhat - bhatcorr);
        hsq(i) = varg / vary;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::left << std::setw(10) << pi(i) 
            << std::setw(10) << int(nnz(i-1)) 
            << std::setw(15) << sigmaSq(i) 
            << std::setw(10) << hsq(i) << std::endl;

        keptIter.row(i-1) << pi(i), int(nnz(i-1)), sigmaSq(i), hsq(i);

    }

    return 0;
};