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

};

class model{

};

class data{

};

class estimation: public MCMC, public model, public data{
    /* 
        Perform bayesian inference, i.e., posterior = likelihood * prior
    */


};

#endif // !INFERENCE_HPP

