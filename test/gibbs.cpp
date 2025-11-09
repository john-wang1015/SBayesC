#include "gibbs.hpp"

void samples::transformation(){

}

void samples::invTransformation(){
    
}

// Constructor implementation
gibbsHMCSampler::gibbsHMCSampler(
    unsigned chainLength, 
    unsigned burnIn, 
    unsigned thin, 
    double stepSize, 
    unsigned numLeapfrogSteps, 
    unsigned dim
) : chainLength(chainLength), burnIn(burnIn), thin(thin), 
    stepSize(stepSize), numLeapfrogSteps(numLeapfrogSteps), dim(dim) {}

// Set the initial state of the sampler
void gibbsHMCSampler::setInitialState(const Eigen::VectorXf& state) {
    init_state = state;
}

// Perform a single HMC step
Eigen::VectorXf gibbsHMCSampler::hmcStep(
    const Eigen::VectorXf& current,
    const std::function<double(const Eigen::VectorXf&)>& logProb,
    const std::function<Eigen::VectorXf(const Eigen::VectorXf&)>& gradient
) {
    // Check dimensions
    if (current.size() != dim) {
        throw std::runtime_error("Dimension mismatch: 'current' size does not match 'dim'");
    }

    // Initialize position and momentum
    Eigen::VectorXf position = current;
    Eigen::VectorXf momentum = Eigen::VectorXf::Random(dim); // Random momentum
    double initialHamiltonian = -logProb(position) + 0.5 * momentum.squaredNorm();

    // Leapfrog integration
    Eigen::VectorXf velocity = momentum;
    for (unsigned i = 0; i < numLeapfrogSteps; ++i) {
        Eigen::VectorXf grad = gradient(position);
        if (grad.size() != dim) {
            throw std::runtime_error("Dimension mismatch: Gradient size does not match 'dim'");
        }

        // Half-step for velocity
        velocity -= 0.5 * stepSize * grad;
        // Full-step for position
        position += stepSize * velocity;
        // Half-step for velocity
        grad = gradient(position);
        velocity -= 0.5 * stepSize * grad;
    }

    // Compute final Hamiltonian
    double finalHamiltonian = -logProb(position) + 0.5 * velocity.squaredNorm();
    double acceptanceProbability = std::exp(initialHamiltonian - finalHamiltonian);

    // Metropolis accept-reject step
    if (static_cast<double>(rand()) / RAND_MAX < acceptanceProbability) {
        return position; // Accept the new state
    } else {
        return current; // Reject and return the old state
    }
}


// Perform Gibbs sampling with HMC
std::vector<Eigen::VectorXf> gibbsHMCSampler::sample(
    const std::vector<std::function<double(const Eigen::VectorXf&)>>& fullConditionals,
    const std::vector<std::function<Eigen::VectorXf(const Eigen::VectorXf&)>>& gradients
) {
    // Initialize variables
    std::vector<Eigen::VectorXf> samples;
    Eigen::VectorXf currentState = init_state;

    // Perform Gibbs sampling
    for (unsigned iter = 0; iter < chainLength; ++iter) {
        for (unsigned i = 0; i < dim; ++i) {
            // HMC sampling for the i-th variable
            currentState = hmcStep(currentState, fullConditionals[i], gradients[i]);
        }

        // Store samples after burn-in and thinning
        if (iter >= burnIn && (iter - burnIn) % thin == 0) {
            samples.push_back(currentState);
        }
    }

    return samples;
}

// Placeholder for custom log-probability computation
double gibbsHMCSampler::computeLogProb(const Eigen::VectorXf& state) {
    // This function could be defined if the class requires a specific log-probability computation
    throw std::runtime_error("computeLogProb is not implemented.");
}

// Placeholder for custom gradient computation
Eigen::VectorXf gibbsHMCSampler::computeGradient(const Eigen::VectorXf& state) {
    // This function could be defined if the class requires a specific gradient computation
    throw std::runtime_error("computeGradient is not implemented.");
}
