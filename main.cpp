#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "data.hpp"
#include "inference.hpp"
#include "gibbs.hpp"

// Define log-probability for the 1st conditional
double logProb1(const Eigen::VectorXf& x) {
    return -0.5 * x(0) * x(0);
}

// Define log-probability for the 2nd conditional
double logProb2(const Eigen::VectorXf& x) {
    return -0.5 * (x(1) - 3) * (x(1) - 3);
}

// Define gradient for the 1st conditional
Eigen::VectorXf gradient1(const Eigen::VectorXf& x) {
    Eigen::VectorXf grad(x.size());
    grad.setZero(); // Ensure all dimensions are initialized
    grad(0) = -x(0); // Only update the relevant dimension
    return grad;
}

// Define gradient for the 2nd conditional
Eigen::VectorXf gradient2(const Eigen::VectorXf& x) {
    Eigen::VectorXf grad(x.size());
    grad.setZero(); // Ensure all dimensions are initialized
    grad(1) = -(x(1) - 3); // Only update the relevant dimension
    return grad;
}


int main() {
    // Initialize sampler
    gibbsHMCSampler sampler(1000, 100, 10, 0.1, 10, 2);

    // Set initial state
    Eigen::VectorXf initState(2);
    initState << 0, 0;
    sampler.setInitialState(initState);

    // Define conditionals and gradients
    std::vector<std::function<double(const Eigen::VectorXf&)>> conditionals = {logProb1, logProb2};
    std::vector<std::function<Eigen::VectorXf(const Eigen::VectorXf&)>> gradients = {gradient1, gradient2};

    // Perform sampling
    std::vector<Eigen::VectorXf> samples = sampler.sample(conditionals, gradients);

    // Print the samples
    for (const auto& sample : samples) {
        std::cout << sample.transpose() << std::endl;
    }
    /*
    std::string binFilePath = "1000G_eur_chr22.ldm.full.bin";
    std::string phenoFilePath = "sim_1.ma";
    readFile inputData;
    inputData.readBinFullLD(binFilePath);
    inputData.readSummary(phenoFilePath);

    std::cout << inputData.B(1,1) << std::endl;
    */
    return 0;
}