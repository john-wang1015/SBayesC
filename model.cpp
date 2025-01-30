#include "model.hpp"

using namespace Eigen;
using namespace std;

void BayesC::reconstruction::approximateD(const Data &data, const VectorXf &se, const VectorXf &bhat, const VectorXf &n){
    unsigned n_size = data.numSNP;
    cout << "number of sample: " << n_size << endl;
    VectorXf diagonalElements(n_size);

    for (unsigned j = 0; j < n_size; ++j) {
        diagonalElements[j] = 1.0f / (se[j] + (bhat[j] * bhat[j] / n[j]));
    }

    D = diagonalElements.asDiagonal();
};

void BayesC::reconstruction::buildXTX(const MatrixXf &B, Data &data){
    MatrixXf D_sqrt = D.cwiseSqrt();
    data.XTX = D_sqrt * B * D_sqrt;
};

void BayesC::reconstruction::buildXTy(const MatrixXf &bhat, Data &data){
    data.XTy = D * bhat;
};

//void BayesC::reconstruction::recover(Data& data) {
    /*
        Using Cholesky decomposition to recover X and y. If simulate data 
        used, it should recover the exact X and y.
    */

//}

void BayesC::SNPEffect::sampleFromPrior(const Data& data, VectorXf &currentState, MatrixXf &histMCMCSamples,const float mean, const float variance, const float pi){
    unsigned beta_size = data.numSNP;
    VectorXf p(2);
    p << 1 - pi, pi;
    for (unsigned i = 0; i < beta_size; i++) {
        unsigned is_included = Stat::Bernoulli::sample(p);

        if (i <= 10){
            cout << is_included <<endl;
        }

        if (is_included){ // if probability 1 - pi, then follow Normal dist
            currentState[i] = Stat::Normal::sample(mean, variance);
        }
        else {
            currentState[i] = 0.0;
        }
    }
    histMCMCSamples.row(0) = currentState.transpose();
};

void BayesC::SNPEffect::fullconditional(const VectorXf y, const MatrixXf X, const VectorXf &currentState, float current_value, const unsigned index, const float sigma_beta2, const float sigma_epsilon2) {
    /*
        Computes the full conditional probability:
        P(beta_j | beta_{-j}, pi, sigma_beta^2, sigma_epsilon^2) = Normal(X_j^T * w / l_jc, sigma_epsilon^2 / l_jc)
        with probability 1 - pi
    */
    float sample_old = current_value;
    VectorXf X_beta = X * currentState;
    VectorXf Xj_betaj = X.col(index) * currentState(index);
    VectorXf w = y - (X_beta - Xj_betaj);

    float l_jc = X.col(index).dot(X.col(index)) + (sigma_epsilon2 / sigma_beta2);

    // Compute beta_hat_j = X_j^T * w / l_jc
    float beta_hat_j = X.col(index).dot(w) / l_jc;
    // Update the state with the new sample using full conditional probability
    current_value = Stat::Normal::sample(beta_hat_j, sigma_epsilon2 / l_jc);
}

void BayesC::SNPEffect::gradient() {


};

void BayesC::Pi::fullconditional() {

};

void BayesC::Pi::gradient() {


}

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

};

void SBayesC::SNPEffect::gradient(){
    /*
        Since the full conditional probability is Gaussian, not necessary to use autodiff library
    */


};