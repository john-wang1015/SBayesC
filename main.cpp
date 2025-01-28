﻿#include <cmath>
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

    BayesC::reconstruction recon;
    recon.approximateD(data, data.se, data.bhat, data.n);

    std::cout << "Diagonal Matrix D:\n" << recon.D.block(0, 0, 10, 10) << std::endl;

    return 0;
}