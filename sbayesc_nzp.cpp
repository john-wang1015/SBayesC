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

//std::random_device rd;
//std::mt19937 rng(123);

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
    numSNP = 2000;  // Expected matrix size
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

void recoverMatrix(MatrixXf& B, VectorXf& b, VectorXf& se, VectorXf& n, unsigned numSNP, MatrixXf& XTX, VectorXf& XTy) {
    unsigned n_size = numSNP;
    VectorXf diagonalElements(n_size);

    for (unsigned j = 0; j < n_size; ++j) {
        diagonalElements[j] = 1.0f / (se[j] + (b[j] * b[j] / n[j]));
    }

    MatrixXf D = diagonalElements.asDiagonal();

    VectorXf D_sqrt_vec = D.diagonal().cwiseSqrt();
    DiagonalMatrix<float, Dynamic> D_sqrt(D_sqrt_vec);
    XTX = D_sqrt * B * D_sqrt;

    XTy = D.diagonal().cwiseProduct(b);
}

void writeVectorsToBinary(const std::string& filename,
    const Eigen::VectorXf& sigmaSq,
    const Eigen::VectorXd& nnz) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Write size of sigmaSq
    int sigmaSqSize = sigmaSq.size();
    file.write(reinterpret_cast<const char*>(&sigmaSqSize), sizeof(int));

    // Write data of sigmaSq
    file.write(reinterpret_cast<const char*>(sigmaSq.data()), sigmaSqSize * sizeof(float));

    // Write size of nnz
    int nnzSize = nnz.size();
    file.write(reinterpret_cast<const char*>(&nnzSize), sizeof(int));

    // Write data of nnz
    file.write(reinterpret_cast<const char*>(nnz.data()), nnzSize * sizeof(double));

    file.close();
}

// Sample from normal distribution
float sample_normal(std::mt19937& rng, float mean, float stddev) {
    std::normal_distribution<float> norm_dist(mean, stddev);
    return norm_dist(rng);
}

// Sample from uniform distribution (0,1)
float sample_uniform(std::mt19937& rng) {
    std::uniform_real_distribution<float> uniform_dist(0.0, 1.0);
    return uniform_dist(rng);
}

// Sample from chi-squared distribution
float sample_chisq(std::mt19937& rng, float dof) {
    std::chi_squared_distribution<float> chi_dist(dof);
    return chi_dist(rng);
}

// Sample from scaled inverse chi-squared distribution
float sample_scaled_inv_chi_squared(std::mt19937& rng, float dof, float scale) {
    std::chi_squared_distribution<float> chi_dist(dof);
    float x = chi_dist(rng);
    return (dof * scale) / x;
}

// Sample from beta distribution
float sample_beta(std::mt19937& rng, float alpha, float beta) {
    std::gamma_distribution<float> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<float> gamma_beta(beta, 1.0);
    float x = gamma_alpha(rng);
    float y = gamma_beta(rng);
    return x / (x + y);
}

// Sample from Bernoulli distribution
float sample_bernoulli(std::mt19937& rng, float p) {
    std::discrete_distribution<int> dist({ 1 - p, p });
    return dist(rng);
}

int main(int argc, char* argv[]) {
    cout << "Starting..." << endl;

    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " <LD matrix file> <GWAS summary file> <output beta file> <output nnz file> <random_seed> <n_size> [n_iter] [pi_init] [hsq_init]" << endl;
        return 1;
    }

    // Read input filenames
    string binFilePath = argv[1];
    string phenoFilePath = argv[2];
    string outputBetaFile = argv[3];
    string outputNnzFile = argv[4];

    // Read numerical inputs
    int random_seed = stoi(argv[5]);
    float n_size = stof(argv[6]); // Convert to float (ensuring 1e4 = 10000.0)

    // Optional hyperparameters
    unsigned n_iter = (argc > 7) ? stoi(argv[7]) : 10000;
    float pi_init = (argc > 8) ? stof(argv[8]) : 0.1;
    float hsq_init = (argc > 9) ? stof(argv[9]) : 0.5;

    // Set the random seed
    mt19937 rng(random_seed);
    cout << "Using random seed: " << random_seed << endl;
    cout << "Using n_size: " << n_size << endl;

    // Load data
    unsigned numSNP;
    VectorXf b, se, n;
    MatrixXf LD;
    readBinTxtFile(binFilePath, numSNP, LD);
    readSummary(phenoFilePath, b, se, n, numSNP);

    // Initialize parameters
    VectorXf pi = VectorXf::Zero(n_iter + 1);
    VectorXf hsq = VectorXf::Zero(n_iter + 1);
    pi(0) = pi_init;
    hsq(0) = hsq_init;

    VectorXf scale = (1.0 / (numSNP * se.array().square())).sqrt();
    VectorXf bhat = b.array() * scale.array();

    float vary = 1.0;
    float varg = hsq(0);
    float vare = vary;
    VectorXf sigmaSq = VectorXf::Zero(n_iter + 1);
    sigmaSq(0) = varg / (numSNP * pi(0));

    float nub = 4.0f, nue = 4.0f;
    float scaleb = (nub - 2) / nub * sigmaSq(0);
    float scalee = (nue - 2) / nue * vare;

    VectorXf beta = VectorXf::Zero(numSNP);
    MatrixXf beta_mcmc = MatrixXf::Zero(n_iter + 1, numSNP);
    VectorXf bhatcorr = bhat;

    MatrixXf keptIter = MatrixXf::Zero(n_iter, 4);
    VectorXd nnz = VectorXd::Zero(n_iter);

    cout << left << setw(10) << "Iteration"
         << left << setw(10) << "pi"
         << setw(10) << "nnz"
         << setw(15) << "sigma_beta"
         << setw(10) << "hsq" << endl;

    for (int i = 1; i < n_iter; i++) {
        float invSigmaSq = 1.0 / sigmaSq(i - 1);
        float ssq = 0;
        VectorXf numSnpDist_current = VectorXf::Zero(2);

        for (int j = 0; j < numSNP; j++) {
            float beta_old = beta(j);
            float rhs = (bhatcorr(j) + beta_old) / (vare / numSNP);
            float invLhs = 1.0f / (1.0f / (vare / numSNP) + invSigmaSq);
            float uhat = invLhs * rhs;

            float logDelta_active = 0.5f * (log(invLhs) - log(sigmaSq(i - 1)) + uhat * rhs) + log(pi(i - 1));
            float logDelta_inactive = -0.5f * log((2 * M_PI) / n_size) - beta(j) * beta(j) * n_size / 2 + log(1.0f - pi(i - 1));
            float pi_current = 1.0 / (1.0 + exp(logDelta_inactive - logDelta_active));

            float delta = sample_bernoulli(rng, pi_current);
            if (delta > 0) {
                numSnpDist_current(1) += 1;
                beta(j) = sample_normal(rng, uhat, sqrt(invLhs));

                bhatcorr += LD.col(j) * (beta_old - beta(j));
                ssq += beta(j) * beta(j);
                nnz(i - 1) += 1;
            } else {
                numSnpDist_current(0) += 1;
                beta(j) = sample_normal(rng, 0.0, sqrt(1 / n_size));
                bhatcorr += LD.col(j) * (beta_old - beta(j));
            }
        }

        VectorXf beta_mcmc_sample = beta.array() / scale.array();
        beta_mcmc.row(i) = beta_mcmc_sample.transpose();

        pi(i) = sample_beta(rng, numSnpDist_current(1), numSnpDist_current(0));

        float sigma_beta = sample_chisq(rng, nnz(i - 1) + nub);
        sigmaSq(i) = (ssq + nub * scaleb) / sigma_beta;

        varg = beta.dot(bhat - bhatcorr);
        hsq(i) = varg / vary;

        if (i % 500 == 0) {
            cout << fixed << setprecision(6);
            cout << left << setw(10) << i
                 << left << setw(10) << pi(i)
                 << setw(10) << int(nnz(i - 1))
                 << setw(15) << sigmaSq(i)
                 << setw(10) << hsq(i) << endl;
        }

        keptIter.row(i - 1) << pi(i), int(nnz(i - 1)), sigmaSq(i), hsq(i);
    }

    // Save output files
    saveMatrixToBinary(outputBetaFile, beta_mcmc);
    writeVectorsToBinary(outputNnzFile, sigmaSq, nnz);

    return 0;
}