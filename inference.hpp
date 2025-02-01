#ifndef INFERENCE_HPP
#define INFERENCE_HPP

#include <iostream>
#include <Eigen/Dense>
#include "model.hpp"
#include "gibbs.hpp"
#include "data.hpp"
#include "stat.hpp"

using namespace std;
using namespace Eigen;

class MCMC {
public:
    unsigned chainLength;
    unsigned chainThin;  // Fixed typo from chinThin
    unsigned numberChain;
    unsigned numberIterations;
    std::string modelUsed;

    virtual void initialState() = 0; // Ensure derived classes implement this
    virtual void runInference() = 0;
};

class inferenceSBayesC : public MCMC {
    public:
        class reconstruction {
            public:
                MatrixXf D;

                void approximateD(const Data &data, const VectorXf &se, const VectorXf &bhat, const VectorXf &n);
                void buildXTX(const MatrixXf &B, Data &data);
                void buildXTy(const MatrixXf &bhat, Data &data);
            };

    public:
        Data data;
        string binFilePath;
        string phenoFilePath;
        MatrixXf histMCMCSamples;  
        MatrixXf r_hist; 
        VectorXf estimatePi_hist;
        VectorXf sigma_beta_hist;
        VectorXf sigma_se_hist;

        reconstruction recon;
        SBayesC model; 

        inferenceSBayesC(const std::string &binFilePath, 
                                   const std::string &phenoFilePath, 
                                   unsigned int num_iterations)
            : binFilePath(binFilePath), 
            phenoFilePath(phenoFilePath),
            model(this->data, num_iterations) 
        {
            this->numberIterations = num_iterations; 
            this->data.readBinFullLD(this->binFilePath);
            this->data.readSummary(this->phenoFilePath);

            this->histMCMCSamples = MatrixXf::Zero(this->numberIterations, this->data.numSNP);
            this->r_hist = MatrixXf::Zero(this->numberIterations, this->data.numSNP);
            this->estimatePi_hist = VectorXf::Zero(this->numberIterations);
            this->sigma_beta_hist = VectorXf::Zero(this->numberIterations);
            this->sigma_se_hist = VectorXf::Zero(this->numberIterations);
        }

        void initialState() override;
        void runInference() override;
};

class inferenceSBayesCadj : public MCMC {
public:
    void initialState(void);
    void runInference(void);
};

class inferenceSBayesCIden : public MCMC {
public:
    void initialState(void);
    void runInference(void);
};

#endif // INFERENCE_HPP
