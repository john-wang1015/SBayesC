#include "model.hpp"

using namespace Eigen;
using namespace std;

void SBayesC::reconstruction::approximateD(const Data &data, const VectorXf &se, const VectorXf &bhat, const VectorXf &n){
    unsigned n_size = data.numSNP;
    cout << "number of sample: " << n_size << endl;
    VectorXf diagonalElements(n_size);

    for (unsigned j = 0; j < n_size; ++j) {
        diagonalElements[j] = 1.0f / (se[j] + (bhat[j] * bhat[j] / n[j]));
    }

    D = diagonalElements.asDiagonal();
};

void SBayesC::reconstruction::buildXTX(const MatrixXf &B, Data &data){
    VectorXf D_sqrt_vec = D.diagonal().cwiseSqrt();
    DiagonalMatrix<float, Dynamic> D_sqrt(D_sqrt_vec);
    data.XTX = D_sqrt * B * D_sqrt;
};

void SBayesC::reconstruction::buildXTy(const MatrixXf &bhat, Data &data){
    data.XTy = D.diagonal().cwiseProduct(bhat);
};

//void BayesC::reconstruction::recover(Data& data) {
    /*
        Using Cholesky decomposition to recover X and y. If simulate data 
        used, it should recover the exact X and y.
    */

//}

void SBayesC::SNPEffect::sampleFromPrior(const Data& data, VectorXf &currentState, MatrixXf &histMCMCSamples,const float mean, const float variance, const float pi){
    unsigned beta_size = data.numSNP;
    VectorXf p(2);
    p << 1 - pi, pi;

    if (currentState.size() != beta_size) {
        currentState = VectorXf::Zero(beta_size);
    }

    for (unsigned i = 0; i < beta_size; i++) {
        unsigned is_included = Stat::Bernoulli::sample(p);

        if (is_included){ // if probability 1 - pi, then follow Normal dist
            currentState[i] = Stat::Normal::sample(mean, variance);
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

void SBayesC::SNPEffect::computeR(const Data& data, const VectorXf currentState, VectorXf &r_current, MatrixXf &r_hist, const unsigned iter){
    unsigned beta_size = data.numSNP;
    VectorXf r_old = r_current;
    VectorXf r_new = VectorXf::Zero(beta_size);

    for (unsigned i = 0; i < beta_size; i++) {
        r_new(i) = r_old(i) + data.XTX(i, i) * currentState(i);
    }

    r_current = r_new;
    r_hist.row(iter) = r_current.transpose();
};

void SBayesC::SNPEffect::fullconditional(const Data &data, const VectorXf &r_current, VectorXf &currentState, const float sigma_beta2, const float sigma_epsilon2, const unsigned j) {
    /*
        Computes the full conditional probability:
        P(beta_j | beta_{-j}, pi, sigma_beta^2, sigma_epsilon^2) = Normal(X_j^T * w / l_jc, sigma_epsilon^2 / l_jc)
        with probability 1 - pi
    */
    
    float l_jc = data.XTX(j,j) + (sigma_epsilon2 / sigma_beta2);

    float beta_hat_j = r_current(j) / l_jc;
    // Update the state with the new sample using full conditional probability
    currentState(j) = Stat::Normal::sample(beta_hat_j, sigma_epsilon2 / l_jc);
}

void SBayesC::SNPEffect::gradient() {


};

void SBayesC::Pi::fullconditional() {
    /*
         Pi ~ Beta(1,1), the full conditional probability is:

    */

};

void SBayesC::Pi::gradient() {


}



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

void SBayesCI::SNPEffect::gradient(){
    /*
        Since the full conditional probability is Gaussian, not necessary to use autodiff library
    */


};