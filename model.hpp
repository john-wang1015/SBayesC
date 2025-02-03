#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include "stat.hpp"
#include "data.hpp"

using namespace std;
using namespace Eigen;

class Parameter {
    // base class for a single parameter
    public:
        const string label;
        float value;   // sampled value
        
        Parameter(const string &label): label(label){
            value = 0.0;
        }
};

class ParamSet {
    // base class for a set of parameters of same kind, e.g. fixed effects, snp effects ...
    public:
        const string label;
        const vector<string> &header;
        const unsigned size;
        VectorXf values;
            
        ParamSet(const string &label, const vector<string> &header)
        : label(label), header(header), size(int(header.size())){
            values.setZero(size);
        }
};

class Model {
    public:
        unsigned numSnps;
        
        vector<ParamSet*> paramSetVec;
        vector<Parameter*> paramVec;
        vector<Parameter*> paramToPrint;
        vector<ParamSet*> paramSetToPrint;
    };


/*
    Original SBayesC model
*/
class SBayesC : public Model {
    public:
        class SNPEffect: public ParamSet, public Stat::Normal, public Stat::Bernoulli {
            public:
                SNPEffect() : ParamSet("SNP Effects", vector<string>()), Stat::Normal(), Stat::Bernoulli() {}
            
                void sampleFromPrior(const Data& data, VectorXf& currentState, MatrixXf &histMCMCSamples, const float current_sigma_beta, const float current_pi);
                void initialR(const Data& data, const MatrixXf &histMCMCSamples, VectorXf &r_current, MatrixXf &r_hist);
                //void computeR(const Data& data, const VectorXf &currentState, VectorXf &r_current);
                void scaleBeta();
                void fullconditional(const Data& data, const VectorXf& r_current, VectorXf& currentState, const float sigma_beta2, const float sigma_epsilon2, const unsigned j);
                //void gradient(); // if use HMC-within-Gibbs
        };

        class Pi: public Parameter, public Stat::Beta{
            public:
                float current_pi;

                Pi(): Parameter("Pi"), Stat::Beta() {
                    current_pi = 0.0f;
                }

                void sampleFromPrior();
                void fullconditional(const Data &data, const float numSnpEff, float current_pi);
                void updatePi();
                //void gradient(); // if use HMC-within-Gibbs
        };

        class EffectVar: public Parameter, public Stat::InvChiSq, public Stat::Gamma {
            public:
                float df;
                float scale;
                float current_sigma_beta;
            
                EffectVar() : Parameter("Effect Variance"), Stat::InvChiSq(), Stat::Gamma() {
                    df = 4.0f;
                    scale = 1.0f;
                    current_sigma_beta = 0.0f;
                }
                
                void sampleFromPrior();
                void fullconditional();
                void scaleEffVar();
        };

        class ResidualVar : public Parameter, public Stat::InvChiSq, public Stat::Gamma {
            public:
                float df;
                float scale;
                float current_sigma_se;

                ResidualVar() : Parameter("Residual Variance"), Stat::InvChiSq(), Stat::Gamma() {
                    df = 4.0f;
                    scale = 1.0f;
                    current_sigma_se = 0.0f;
                }
        
                void sampleFromPrior();
                void fullconditional();
                void scaleResVar();
        };

        class Heritability: public Parameter {
            public:
                Heritability() : Parameter("hsq") {}

                void compute();
        };

        class NumNonZeroSNP : public Parameter {
            public:
                unsigned nonZeroEff;
                unsigned zeroEff;

                NumNonZeroSNP() : Parameter("Number of Non-zero snp") {}

                void countZero(const Data &data, const VectorXf &currentState){
                    zeroEff = (currentState.array() == 0).count();
                    nonZeroEff = data.numSNP - zeroEff;
                };
        };

        // need have some method to scale the parameters values

    public:
        const Data &data;
        unsigned num_iterations;
        VectorXf currentState;              // store current samples
        VectorXf betaHat;
        //MatrixXf histMCMCSamples;           // store all samples
        VectorXf r_current;
        MatrixXf r_hist;
        VectorXf estimatePi;                  // pi value
        VectorXf sigma_beta;
        VectorXf sigma_se;

        SNPEffect snpEffect;
        Pi pi;
        EffectVar effectVar;
        ResidualVar residualVar;
        Heritability hsq;
        NumNonZeroSNP numNonZeroSNP;

        SBayesC(const Data& data, unsigned num_iterations) : 
            data(data), 
            snpEffect(), // Explicitly call constructor
            pi(), 
            effectVar(), 
            residualVar(), 
            hsq(), 
            numNonZeroSNP() 
        {
            unsigned beta_size = data.numSNP;
            sigma_beta  = VectorXf::Zero(num_iterations);
            estimatePi  = VectorXf::Zero(num_iterations);

            currentState = VectorXf::Zero(beta_size);
            //histMCMCSamples = MatrixXf::Zero(num_iterations, beta_size);
            r_current = VectorXf::Zero(beta_size);
            //r_hist = MatrixXf::Zero(num_iterations, beta_size);
        }

};


/*
    SBayesC model with zero effect SNP ~ N(0,1/ sample size)
*/
class SBayesC_adj_prior : public SBayesC {
    public:
        class SNPEffect : public SBayesC::SNPEffect {
            public:
                
        };


};

/*
    SBayesC model with R*sigma_beta where R is identity matrix not D^-1*X^T
*/
class SBayesCI: public Model{
    public:
        class SNPEffect: public ParamSet, public Stat::Normal{
            // SNP effect beta_j
            public:

                void fullConditional(const VectorXf &r_adjust,const MatrixXf XTX, VectorXf &beta_current, const float sigma_e, const float sigma_beta);
                void gradient();
        };

        class Pi: public Parameter, public Stat::Beta{
            // mixture component, pi
            public:

                void update();
        };

        class SNPEffectVar: public Parameter, public Stat::InvChiSq{
            // sigma_beta^2
            public:

                void sampleFromPrior();
                void fullConditional();
                void gradient();
        };

        class residualVar: public Parameter, public Stat::InvChiSq{
            // sigma_epsilon^2
            public:

                void sampleFromPrior();
                void fullConditional();
                void gradient();
        };

    public:
        const Data &data;
        VectorXf beta_current; // store current state for beta values
        MatrixXf beta_history; // store all beta values
        VectorXf r_adjust; 



};

/*
    Combine the change of SBayesC_adj_prior and SBayesCI
*/
class SBayesC_new : public SBayesC {
    public:
        class SNPEffect : public SBayesC::SNPEffect {
            public:
                
        };


};


#endif //MODEL_HPP
