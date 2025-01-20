#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include "stat.hpp"
#include "data.hpp"

using namespace std;

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
        
        virtual void sampleUnknowns(void) = 0;
        virtual void sampleStartVal(void) = 0;
    };

//class BayesC: public Model{
//}

class SBayesC: public Model{
    public:
        class fixedEffect: public ParamSet, public Stat::Flat{
            public:

                void fullConditional();
                void gradient();
        };

};

#endif // !MODEL_HPP
