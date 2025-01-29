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


class BayesC : public Model {
    public:
        class reconstruction{
            /*
                Class use to reconstruct XTX and XTy based on summary 
                statistics data. Individual level model can also use 
                this as input. 
            */
            public:
                MatrixXf D;

                void approximateD(const Data &data, const VectorXf &se, const VectorXf &bhat, const VectorXf &n);
                void buildXTX(const MatrixXf &B, Data &data);
                void buildXTy(const MatrixXf &bhat, Data &data);
                void recover(Data& data);
        };

        class SNPEffect: public ParamSet, public Stat::Normal, public Stat::Uniform {
            public:
                void sampleFromPrior(const Data& data, VectorXf& currentState, MatrixXf &histMCMCSamples, const float mean, const float variance, const float pi);
                void fullconditional(const VectorXf y, const MatrixXf X, const unsigned index, const unsigned iteration);
                void gradient(); // if use HMC-within-Gibbs
        };

        class Pi: public Parameter, public Stat::Beta, public Stat::Bernoulli{
            public:
                void sampleFromPrior();
                void fullconditional();
                void gradient(); // if use HMC-within-Gibbs
        };

        class EffectVar: public Parameter, public Stat::InvChiSq {

        };

        class ResidualVar : public Parameter, public Stat::InvChiSq {

        };

        class Heritability: public Parameter {

        };

        class NumNonZeroSNP : public Parameter {

        };

        // need have some method to scale the parameters values

    public:
        const Data &data;
        VectorXf currentState;              // store current samples
        MatrixXf histMCMCSamples;           // store all samples
        MatrixXf X;                         // recovered genotype
        MatrixXf y;                         // recovered phenotype

        reconstruction recon;
        SNPEffect snpEffect;
        Pi pi;
        EffectVar effectVar;
        ResidualVar residualVar;
        Heritability hsq;
        NumNonZeroSNP numNonZeroSNP;

};

class BayesC_adj_prior : public BayesC {
    public:
        class SNPEffect : public BayesC::SNPEffect {
            public:
            
        };


};


class SBayesC: public Model{
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

#endif //MODEL_HPP
