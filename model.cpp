#include "model.hpp"

void SBayesC::SNPEffect::fullConditional(const VectorXf &r_adjust,const MatrixXf XTX, VectorXf &beta_current, const float sigma_e, const float sigma_beta){
    /*
        function of funll conditional probability for beta_j
    */
    for (int i = 0; i < r_adjust.size(); i++){
        float r_j = r_adjust[i];
        float l_jc = XTX(i, i) + (sigma_e / sigma_beta);
        
        float mean_j = r_j / l_jc;
        float stddev_j = std::sqrt(sigma_e / l_jc);

        beta_current[i] = Stat::Normal::sample(mean_j, stddev_j);
    }

}

void SBayesC::SNPEffect::gradient(){
    /*
        Since the full conditional probability is Gaussian, not necessary to use autodiff library
    */


};