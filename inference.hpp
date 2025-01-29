#ifndef INFERENCE_HPP
#define INFERENCE_HPP

#include <iostream>
#include <Eigen/Dense>
#include "model.hpp"
#include "gibbs.hpp"
#include "data.hpp"

using namespace std;
using namespace Eigen;

class MCMC{
    /*
         Multiple MCMC algorithm provide:
            1. vanilla Hamiltonian monte carlo (vHMC)-within-Gibbs sampler
            2. NUTS-within-Gibbs sampler (improve Alg 1)
            3. Skinny Gibbs
    */
    public:
        unsigned chainLength;
        unsigned chinThin;
        unsigned numberChain;
        unsigned numberIterations;

        virtual void initialState(void) = 0;
        virtual void sampleSingleStep(void) = 0;
};

class Model{
    /*
        Select the model used for inference, options are:
            1. BayesC
            2. BayesC_adj_prior
            3. SBayesC with R is identity matrix
    */
    public:
        string modelUsed;

        virtual void runInference(void) = 0;
};

class inferenceBayesC: public MCMC, public Model{
    public:
        void initialState() override {

        }

        void sampleSingleStep() override {

        }

        void runInference() override {

        }
};

class inferenceBayesCadj : public MCMC, public Model {
    public:
        void initialState() override {

        }

        void sampleSingleStep() override {

        }

        void runInference() override {

        }
};

class inferenceSBayesCIden : public MCMC, public Model {
    public:
        void initialState() override {

        }

        void sampleSingleStep() override {

        }

        void runInference() override {

        }
};

#endif // INFERENCE_HPP

