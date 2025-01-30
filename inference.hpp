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
        string modelUsed;

        virtual void initialState(void) = 0;
        virtual void sampleSingleStep(void) = 0;
        virtual void runInference(void) = 0;
};


class inferenceSBayesC: public MCMC, public Stat::Bernoulli, public SBayesC {
    public:
        void initialState(void){
            unsigned numSNP = Model::numSnps;
            
            
        }

        void sampleSingleStep(void){

        }

        void runInference(void){

        }
};

class inferenceSBayesCadj : public MCMC{
    public:
        void initialState(void){

        }

        void sampleSingleStep(void){

        }

        void runInference(void){

        }
};

class inferenceSBayesCIden : public MCMC{
    public:
        void initialState(void){

        }

        void sampleSingleStep(void){

        }

        void runInference(void){

        }
};

#endif // INFERENCE_HPP

