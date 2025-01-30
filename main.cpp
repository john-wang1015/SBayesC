#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "data.hpp"
#include "inference.hpp"

int main() {
    std::string binFilePath = "1000G_eur_chr22.ldm.full.bin";
    std::string phenoFilePath = "sim_1.ma";
    Data data;
    data.readBinFullLD(binFilePath);
    data.readSummary(phenoFilePath);

    SBayesC::reconstruction recon;
    recon.approximateD(data, data.se, data.bhat, data.n);
    recon.buildXTX(data.B, data);
    recon.buildXTy(data.bhat, data);

    //std::cout << "Diagonal Matrix D:\n" << recon.D.block(0,0,10,10) << std::endl;

    unsigned num_iterations = 10000;  // Example: 10,000 MCMC iterations
    SBayesC model(data, num_iterations);    

    // Define parameters for sampling
    float mean = 0.0;
    float variance = 1.0;
    float pi = 0.5;

    // Call sampleFromPrior
    model.snpEffect.sampleFromPrior(data, model.currentState, model.histMCMCSamples, mean, variance, pi);
    model.snpEffect.initialR(data, model.histMCMCSamples, model.r_current, model.r_hist);

    // Print results
    std::cout << "Sampled first 10 SNP effects:\n" << model.currentState.block(0, 0, 10, 1) << std::endl;
    std::cout << "Sampled first 10 values for r:\n" << model.r_current.block(0, 0, 10, 1) << std::endl;

    return 0;
}