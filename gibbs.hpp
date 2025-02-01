#ifndef GIBBS_HPP
#define GIBBS_HPP

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>

class samples{
    public:
        unsigned numChain;
        

        void transformation();
        void invTransformation();

};

class gibbsSampler: public samples {
    /*
        Simple Gibbs sampler
    */

   public:
        


};

class HMCSampler: public samples {
    /*
        Class for performing Gibbs sampling with Hamiltonian Monte Carlo (HMC) 
        for full conditional probability sampling.
    */
    public:
        // Parameters
        unsigned chainLength; // Total number of samples in the chain
        unsigned burnIn;      // Number of samples to discard as burn-in
        unsigned thin;        // Thinning interval (e.g., save every nth sample)

        Eigen::VectorXf init_state; // Initial state of the Markov chain
        double stepSize;            // Step size for HMC leapfrog integration
        unsigned numLeapfrogSteps;  // Number of leapfrog steps for HMC
        unsigned dim;               // Dimensionality of the problem

        // Constructors
        HMCSampler(
            unsigned chainLength, 
            unsigned burnIn, 
            unsigned thin, 
            double stepSize, 
            unsigned numLeapfrogSteps, 
            unsigned dim
        ) : chainLength(chainLength), burnIn(burnIn), thin(thin), 
            stepSize(stepSize), numLeapfrogSteps(numLeapfrogSteps), dim(dim) {}

        // Set the initial state of the sampler
        void HMCSampler::setInitialState(const Eigen::VectorXf& state) {
            init_state = state;
        }

        // Member functions
        void setInitialState(const Eigen::VectorXf& state); // Set initial state
        std::vector<Eigen::VectorXf> sample(
            const std::vector<std::function<double(const Eigen::VectorXf&)>>& fullConditionals,
            const std::vector<std::function<Eigen::VectorXf(const Eigen::VectorXf&)>>& gradients
        ); // Perform Gibbs sampling with HMC

        // Perform a single HMC step
        Eigen::VectorXf hmcStep(
            const Eigen::VectorXf& current,
            const std::function<double(const Eigen::VectorXf&)>& logProb,
            const std::function<Eigen::VectorXf(const Eigen::VectorXf&)>& gradient
        ); 
        double computeLogProb(const Eigen::VectorXf& state); // Compute log-probability
        Eigen::VectorXf computeGradient(const Eigen::VectorXf& state); // Compute gradient
};




#endif //GIBBS_HPP
