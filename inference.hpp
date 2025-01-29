#ifndef INFERENCE_HPP
#define INFERENCE_HPP

#include <iostream>
#include <Eigen/Dense>
#include "model.hpp"
#include "gibbs.hpp"
#include "data.hpp"

using namespace std;
using namespace Eigen;

//class MCMC{
    /*
         Multiple MCMC algorithm provide:
            1. vanilla Hamiltonian monte carlo (vHMC)-within-Gibbs sampler
            2. NUTS-within-Gibbs sampler (improve Alg 1)
            3. Skinny Gibbs
    */

//};

//class Model{

//};

//class inference: public MCMC, public Model{
    /* 
        Perform bayesian inference, i.e., posterior = likelihood * prior
    */


//};

#endif // INFERENCE_HPP

