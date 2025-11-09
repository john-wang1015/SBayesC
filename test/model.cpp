#include "model.hpp"

using namespace Eigen;
using namespace std;

void SBayesC::SNPEffect::sampleFromPrior(const Data& data, VectorXf &currentState, MatrixXf &histMCMCSamples, const float current_sigma_beta, const float current_pi){
    unsigned beta_size = data.numSNP;
    VectorXf p(2);
    p << 1 - current_pi, current_pi;

    if (currentState.size() != beta_size) {
        currentState = VectorXf::Zero(beta_size);
    }

    for (unsigned i = 0; i < beta_size; i++) {
        unsigned is_included = Stat::Bernoulli::sample(p);

        if (is_included){ // if probability 1 - pi, then follow Normal dist
            currentState[i] = Stat::Normal::sample(0.0, current_sigma_beta);
        }
        else {
            currentState[i] = 0.0;
        }
    }

    if (histMCMCSamples.rows() == 0 || histMCMCSamples.cols() != beta_size) {
        std::cerr << "Error: histMCMCSamples is not properly initialized!" << std::endl;
    } else {
        histMCMCSamples.row(0) = currentState.transpose();
    }
};

void SBayesC::SNPEffect::initialR(const Data& data, const MatrixXf &histMCMCSamples, VectorXf &r_current, MatrixXf &r_hist){
    unsigned beta_size = data.numSNP;

    if (histMCMCSamples.rows() == 0) {
        std::cerr << "Error: histMCMCSamples is empty!" << std::endl;
        return;
    }

    VectorXf beta_init = histMCMCSamples.row(0).transpose();
    r_current.resize(beta_size);
    r_current = data.XTy - data.XTX * beta_init;

    if (r_hist.rows() == 0 || r_hist.cols() != beta_size) {
        std::cerr << "Warning: r_hist is not initialized. Resizing it now.\n";
        r_hist = MatrixXf::Zero(histMCMCSamples.rows(), beta_size);
    }

    r_hist.row(0) = r_current.transpose();
}

/*
void SBayesC::SNPEffect::computeR(const Data& data, const VectorXf &currentState, VectorXf &r_current){
    unsigned beta_size = data.numSNP;
    VectorXf r_old = r_current;
    VectorXf r_new = VectorXf::Zero(beta_size);

    for (unsigned i = 0; i < beta_size; i++) {
        r_new(i) = r_old(i) + data.XTX(i, i) * currentState(i);
    }

    r_current = r_new;
};
*/

void SBayesC::SNPEffect::scaleBeta(){
    
}

void SBayesC::SNPEffect::fullconditional(const Data &data, const VectorXf &r_current, VectorXf &currentState, const float sigma_beta2, const float sigma_epsilon2, const unsigned j) {
    /*
        Computes the full conditional probability:
        P(beta_j | beta_{-j}, pi, sigma_beta^2, sigma_epsilon^2) = Normal(X_j^T * w / l_jc, sigma_epsilon^2 / l_jc)
        with probability 1 - pi
    */
    
    float l_jc = data.XTX(j,j) + (sigma_epsilon2 / sigma_beta2);

    float beta_hat_j = r_current(j) / l_jc;

    currentState(j) = Stat::Normal::sample(beta_hat_j, sigma_epsilon2 / l_jc);
}

//void SBayesC::SNPEffect::gradient() {
    /*
        update it when we need to 
    */

//};

void SBayesC::Pi::sampleFromPrior(){
    current_pi = Stat::Beta::sample(1.0,1.0);
}

void SBayesC::Pi::fullconditional(const Data &data, const float numSnpEff, float current_pi) {
    float alphaTilde = numSnpEff + 1.0;
    float betaTilde  = data.numSNP - numSnpEff + 1.0;
    current_pi = Beta::sample(alphaTilde, betaTilde);
};

//void SBayesC::Pi::gradient() {
//}

void SBayesC::EffectVar::sampleFromPrior(){
    current_sigma_beta = InvChiSq::sample(df, scale);
};


void SBayesC::EffectVar::fullconditional(){

};

void SBayesC::EffectVar::scaleEffVar(){

};


void SBayesC::ResidualVar::sampleFromPrior(){
    current_sigma_se = InvChiSq::sample(df, scale);
};

void SBayesC::ResidualVar::fullconditional(){

};

void SBayesC::ResidualVar::scaleResVar(){

};

void SBayesCI::SNPEffect::fullConditional(const VectorXf &r_adjust,const MatrixXf XTX, VectorXf &beta_current, const float sigma_e, const float sigma_beta){
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

};

//void SBayesCI::SNPEffect::gradient(){
    /*
        Since the full conditional probability is Gaussian, not necessary to use autodiff library
    */
//};