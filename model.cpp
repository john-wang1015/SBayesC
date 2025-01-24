#include "model.hpp"

void SBayesC::SNPEffect::fullConditional(const VectorXf &r_adjust, const float sigma_e, const float sigma_beta, const MatrixXf XTX, const int j_index, VectorXf beta_current){
    /*
        function of funll conditional probability for beta_j
    */
    float r_j = r_adjust[j_index];
    if (XTX(j_index, j_index) == 0) {
        throw std::runtime_error("Diagonal element of XTX at index " + std::to_string(j_index) + " is zero, causing division by zero.");
    }
    float l_jc = XTX(j_index, j_index) + (sigma_e / sigma_beta);

    float mean_j = r_j / l_jc;
    float stddev_j = std::sqrt(sigma_e / l_jc);

    beta_current[j_index] = Stat::Normal::sample(mean_j, stddev_j);
}

void SBayesC::SNPEffect::gradient(){
    /*
        Since the full conditional probability is Gaussian, not necessary to use autodiff library
    */
   

};