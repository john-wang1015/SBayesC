#include "gibbs.hpp"

void samples::transformation(){

}

void samples::invTransformation(){
    
}


Eigen::VectorXf HMCSampler::hmcStep(
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


std::vector<Eigen::VectorXf> HMCSampler::sample(
    const std::vector<std::function<double(const Eigen::VectorXf&)>>& fullConditionals,
    const std::vector<std::function<Eigen::VectorXf(const Eigen::VectorXf&)>>& gradients
) {
    std::vector<Eigen::VectorXf> samples;
    Eigen::VectorXf currentState = init_state;

    // Perform Gibbs sampling
    for (unsigned iter = 0; iter < chainLength; ++iter) {
        for (unsigned i = 0; i < dim; ++i) {
            currentState = hmcStep(currentState, fullConditionals[i], gradients[i]);
        }

        if (iter >= burnIn && (iter - burnIn) % thin == 0) {
            samples.push_back(currentState);
        }
    }

    return samples;
}

double HMCSampler::computeLogProb(const Eigen::VectorXf& state) {
    throw std::runtime_error("computeLogProb is not implemented.");
}

Eigen::VectorXf HMCSampler::computeGradient(const Eigen::VectorXf& state) {
    throw std::runtime_error("computeGradient is not implemented.");
}
