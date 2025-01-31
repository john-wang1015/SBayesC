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
    unsigned num_iterations = 10000;
    

    inferenceSBayesC infer(binFilePath, phenoFilePath, num_iterations);
    infer.initialState();

    return 0;
}


