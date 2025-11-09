#ifndef DATA_HPP
#define DATA_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

class SnpInfo{
    public:
        const string ID;
        const int chrom;
        const int GenPos;
        const int PhysPos;
        const string A1;
        const string A2;
        const float A1Freq;
        const int Index;

};

class Data {
    public:
        unsigned numSNP;
        
        MatrixXf XTX;       // X'X matrix
        VectorXf XTy;       // X'y matrix
        MatrixXf Ddiag;     // diag matrix need to be approximated
        MatrixXf B;         // LD matrix
        VectorXf bhat;         // estimated marginal effect, i.e., b_hat
        VectorXf se;        // standard error for b_hat
        VectorXf n;         // sample size

        void readSummary(const std::string& summaryFilePath);
        void readBinInfo(const std::string& binInfoFilePath);
        void readBinFullLD(const std::string& binFilePath);
        void readBinSparseLD(const std::string& binFilePath);                                     
        void readBinShrunkLD(const std::string& binFilePath);                                      
};

#endif // DATA_HPP
