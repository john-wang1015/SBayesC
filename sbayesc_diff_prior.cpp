#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric> 
#include <vector>
#include <random>
#include <chrono>
#include <stdexcept>
#include <Eigen/Dense>
#include <boost/random.hpp>

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
    cout << "Time for read .bin file: " << elapsed.count() << " seconds" << std::endl;
}


void readBinTxtFile(const std::string& binFilePath, unsigned& numSNP, Eigen::MatrixXf& B) {
    std::ifstream file(binFilePath);
    if (!file.is_open()) { 
        std::cerr << "Error: Cannot open file " << binFilePath << std::endl;
        return;
    }

    std::vector<std::vector<float>> tempData;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::vector<float> values;
        float value;

        while (ss >> value) {
            values.push_back(value);
        }

        if (!values.empty()) {
            tempData.push_back(values);
        }
    }

    file.close();

    // Determine numSNP (assuming square matrix)
    numSNP = tempData.size();

    if (numSNP == 0 || tempData[0].size() != numSNP) {
        std::cerr << "Error: The input file does not contain a valid square matrix." << std::endl;
        return;
    }

    // Initialize Eigen matrix
    B = Eigen::MatrixXf(numSNP, numSNP);

    for (unsigned i = 0; i < numSNP; i++) {
        for (unsigned j = 0; j < numSNP; j++) {
            B(i, j) = tempData[i][j];
        }
    }
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

    cout << "Summary statistics loaded: " << b.size() << " entries." << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Time to read GWAS summary statistics file: " << elapsed.count() << " seconds" << std::endl;
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

void savePosMean(const string& filename, const VectorXf& matrix) {
    ofstream outFile(filename, ios::binary);
    if (!outFile) {
        cerr << "Error: Unable to open file " << filename << " for writing." << endl;
        return;
    }
    int rows = matrix.size();
    outFile.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(matrix.data()), rows * sizeof(float));
    outFile.close();
}

void saveVectorToBinary(const string& filename, const VectorXf& data) {
    ofstream file(filename, ios::binary);
    if (!file) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }
    int size = data.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(int));
    file.write(reinterpret_cast<const char*>(data.data()), size * sizeof(float));
    file.close();
}

// Sample from normal distribution
float sample_normal(boost::random::mt19937& rng, float mean, float stddev) {
    boost::random::normal_distribution<float> norm_dist(mean, stddev);
    return norm_dist(rng);
}

// Sample from uniform distribution (0,1)
float sample_uniform(boost::random::mt19937& rng) {
    boost::random::uniform_real_distribution<float> uniform_dist(0.0, 1.0);
    return uniform_dist(rng);
}

// Sample from chi-squared distribution
float sample_chisq(boost::random::mt19937& rng, float dof) {
    boost::random::chi_squared_distribution<float> chi_dist(dof);
    return chi_dist(rng);
}

// Sample from scaled inverse chi-squared distribution
float sample_scaled_inv_chi_squared(boost::random::mt19937& rng, float dof, float scale) {
    boost::random::chi_squared_distribution<float> chi_dist(dof);
    float x = chi_dist(rng);
    return (dof * scale) / x;
}

float sample_beta(boost::random::mt19937& rng, float alpha, float beta) {
    constexpr float EPSILON = 1e-20f;  // Smallest allowed value

    // Ensure alpha and beta are strictly positive
    alpha = std::max(alpha, EPSILON);
    beta = std::max(beta, EPSILON);

    boost::random::gamma_distribution<float> gamma_alpha(alpha, 1.0);
    boost::random::gamma_distribution<float> gamma_beta(beta, 1.0);

    float x = gamma_alpha(rng);
    float y = gamma_beta(rng);
    
    return x / (x + y);
}


// Sample from Bernoulli distribution (returns int 0 or 1)
int sample_bernoulli(boost::random::mt19937& rng, float p) {
    boost::random::bernoulli_distribution<> bernoulli_dist(p);
    return bernoulli_dist(rng) ? 1 : 0;
}

std::vector<float> sample_dirichlet(boost::random::mt19937& rng, const std::vector<float>& alpha) {
    constexpr float EPSILON = 1e-6f;  // Smallest allowed value
    std::vector<float> gamma_samples(alpha.size());
    float sum = 0.0f;

    // Sample from Gamma distribution for each alpha_i
    for (size_t i = 0; i < alpha.size(); ++i) {
        float alpha_i = std::max(alpha[i], EPSILON); // Ensure positive alpha
        boost::random::gamma_distribution<float> gamma_dist(alpha_i, 1.0);
        gamma_samples[i] = gamma_dist(rng);
        sum += gamma_samples[i];
    }

    // Normalize to get Dirichlet samples
    for (size_t i = 0; i < alpha.size(); ++i) {
        gamma_samples[i] /= sum;
    }

    return gamma_samples;
}

int main(int argc, char* argv[]) {
    cout << "starting...."<<std::endl;
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
    boost::random::mt19937 rng;
    rng.seed(random_seed);
    //mt19937 rng(random_seed);
    cout << "Using random seed: " << random_seed << endl;
    cout << "Using n_size: " << n_size << endl;

    // Load data
    unsigned numSNP;
    VectorXf b, se, n;
    MatrixXf LD;
    readBinTxtFile(binFilePath, numSNP, LD);
    readSummary(phenoFilePath, b, se, n, numSNP);

    // Initialize parameters
    VectorXf pi = VectorXf::Zero(n_iter+1);
    VectorXf hsq = VectorXf::Zero(n_iter+1);
    pi(0) = pi_init;
    hsq(0) = hsq_init;

    VectorXf scale = (1.0 / (n_size * se.array().square())).sqrt();
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
    MatrixXf beta_mcmc = MatrixXf::Zero(n_iter, numSNP);
    VectorXf bhatcorr = bhat;
    MatrixXf keptIter = MatrixXf::Zero(n_iter, 4);
    VectorXd nnz = VectorXd::Zero(n_iter);

    cout << left << setw(10) << "Iteration"
         << left << setw(10) << "pi"
         << setw(10) << "nnz"
         << setw(15) << "sigmaSq"
         << setw(10) << "hsq" << endl;

    for (int i = 0; i < n_iter; i++) {
        float invSigmaSq = 1.0 / sigmaSq(i);
        float logSigmaSq = log(sigmaSq(i));
        float logpi = log(pi(i));
        float logpic = log(1.0f - pi(i));
        float ssq = 0;
        float ssq0 = 0;
        VectorXf numSnpDist_current = VectorXf::Zero(2);
        VectorXf uhat_temp = VectorXf::Zero(10);

        for (int j = 0; j < numSNP; j++) {
            float beta_old = beta(j);
            float rhs = (bhatcorr(j) + beta_old) / (vare / n_size);
            float invLhs1 = 1.0 / (1.0 / (vare / n_size) + invSigmaSq);
            float uhat1 = invLhs1 * rhs;

            float invLhs0 = 1.0 / (1.0 / (vare / n_size) + n_size);
            float uhat0 = invLhs0 * rhs;

            float logDelta_active = 0.5 * (log(invLhs1) - logSigmaSq + uhat1 * rhs) + logpi;
            float logDelta_inactive = 0.5 * (log(invLhs0) - logSigmaSq + uhat0 * rhs) + logpic;
            float pi_current = 1.0 / (1.0 + exp(logDelta_inactive - logDelta_active));

            int delta = sample_bernoulli(rng, pi_current);
            if (delta > 0) {
                numSnpDist_current(1) += 1;
                beta(j) = sample_normal(rng, uhat1, sqrt(invLhs1)); 
                bhatcorr += LD.col(j) * (beta_old - beta(j));
                ssq += beta(j) * beta(j);
                nnz(i) += 1;
            } else {
                numSnpDist_current(0) += 1;
                beta(j) = sample_normal(rng, uhat0, sqrt(invLhs0)); 
                bhatcorr += LD.col(j) * (beta_old - beta(j));  
            }
        }

        VectorXf beta_mcmc_sample = beta.array() / scale.array();
        beta_mcmc.row(i) = beta_mcmc_sample.transpose();

        pi(i+1) = sample_beta(rng, numSnpDist_current(1)+1, numSnpDist_current(0)+1);
        
        float sigma_beta = sample_chisq(rng, nnz(i) + nub);
        sigmaSq(i+1) = (ssq + nub * scaleb) / sigma_beta;
        
        VectorXf scaled_bhat = (bhat.array() - bhatcorr.array()) / scale.array();
        VectorXf beta_scale = beta.array() / scale.array();
        varg = beta_scale.dot(scaled_bhat);
        hsq(i) = varg / vary;

        if (i % 100 == 0) {
            cout << fixed << setprecision(6);
            cout << left << setw(10) << i
                 << left << setw(10) << pi(i)
                 << setw(10) << int(nnz(i))
                 << setw(20) << sigmaSq(i)
                 << setw(10) << hsq(i) << endl;
        }

        keptIter.row(i) << pi(i), int(nnz(i)), sigmaSq(i), hsq(i);
    }

    int start_row = 5000, end_row = n_iter;
    MatrixXf beta_mcmc_subset = beta_mcmc.block(start_row, 0, end_row - start_row, numSNP);
    VectorXf beta_posterior_mean = beta_mcmc_subset.colwise().mean();
    float pi_posterior_mean = pi.segment(start_row, end_row - start_row).mean();
    float hsq_posterior_mean = hsq.segment(start_row, end_row - start_row).mean();
    int nnz_pos_mean = std::round(nnz.segment(start_row, end_row - start_row).mean());

    cout << "posterior mean: "<<endl;
    cout << left << setw(10) << "pi1"
         << left << setw(10) << "pi2"
         << setw(10) << "nnz"
         << setw(20) << "sigmaSq"
         << setw(10) << "h2" << endl;
    cout << fixed << setprecision(6);
            cout << left << setw(10) << 1-pi_posterior_mean
                 << left << setw(10) << pi_posterior_mean
                 << setw(10) << nnz_pos_mean
                 << setw(20) << sigmaSq.segment(start_row, end_row - start_row).mean()
                 << setw(10) << hsq_posterior_mean << endl;

    // Save the posterior mean vectors to binary files
    string outputPosteriorBetaFile = outputBetaFile + "_posterior_mean_beta.bin";
    string outputPosteriorPiFile = outputBetaFile + "_posterior_mean_pi.bin";
    string outputPosteriorHsqFile = outputBetaFile + "_posterior_mean_hsq.bin";

    saveVectorToBinary(outputPosteriorBetaFile, beta_posterior_mean);
    saveVectorToBinary(outputPosteriorPiFile, VectorXf::Constant(1, pi_posterior_mean));
    saveVectorToBinary(outputPosteriorHsqFile, VectorXf::Constant(1, hsq_posterior_mean));

    cout << "Posterior means computed and saved to:\n"
         << outputPosteriorBetaFile << "\n"
         << outputPosteriorPiFile << "\n"
         << outputPosteriorHsqFile << endl;

    // Save other output files
    saveMatrixToBinary(outputBetaFile, beta_mcmc);
    writeVectorsToBinary(outputNnzFile, sigmaSq, nnz);

    return 0;
}
