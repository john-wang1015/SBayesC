#include "inference.hpp"
#include "model.hpp"

void inferenceSBayesC::reconstruction::approximateD(const Data &data, const VectorXf &se, const VectorXf &bhat, const VectorXf &n){
    unsigned n_size = data.numSNP;
    cout << "number of sample: " << n_size << endl;
    VectorXf diagonalElements(n_size);

    for (unsigned j = 0; j < n_size; ++j) {
        diagonalElements[j] = 1.0f / (se[j] + (bhat[j] * bhat[j] / n[j]));
    }

    D = diagonalElements.asDiagonal();
};

void inferenceSBayesC::reconstruction::buildXTX(const MatrixXf &B, Data &data){
    VectorXf D_sqrt_vec = D.diagonal().cwiseSqrt();
    DiagonalMatrix<float, Dynamic> D_sqrt(D_sqrt_vec);
    data.XTX = D_sqrt * B * D_sqrt;
};

void inferenceSBayesC::reconstruction::buildXTy(const MatrixXf &bhat, Data &data){
    data.XTy = D.diagonal().cwiseProduct(bhat);
};


void inferenceSBayesC::initialState() {
    this->recon.approximateD(this->data, this->data.se, this->data.bhat, this->data.n);
    this->recon.buildXTX(this->data.B, this->data);
    this->recon.buildXTy(this->data.bhat, this->data);

    //SBayesC model(this->data, this->numberIterations);   
    model.pi.sampleFromPrior();
    model.effectVar.sampleFromPrior();
    model.residualVar.sampleFromPrior();
    model.snpEffect.sampleFromPrior(this->data, model.currentState, this->histMCMCSamples, model.effectVar.current_sigma_beta, model.pi.current_pi);
    model.snpEffect.initialR(this->data, this->histMCMCSamples, model.r_current, this->r_hist); 

    /*
    cout << "Initial Pi is: " << model.pi.current_pi << endl;
    cout << "Sigma_beta value is: " << model.effectVar.current_sigma_beta << endl;
    cout << "Sigma_se value is: " << model.residualVar.current_sigma_se << endl;
    std::cout << "Sampled first 10 SNP effects:\n" << this->histMCMCSamples.block(0, 0, 10, 10) << std::endl;
    std::cout << "Sampled first 10 values for r:\n" << this->r_hist.block(0, 0, 10, 10) << std::endl;
    */

    /*
        Pending: !store all values to history matrix
    */   
}

void inferenceSBayesC::runInference() {

    for (int i = 0; i < numberIterations; i++){
        // store the r_current and beta_current as vectorXf from Matrx of history
        VectorXf beta_old = this->histMCMCSamples.row(i);
        VectorXf r_old = this->r_hist.row(i);
        float pi_old = model.pi.current_pi;

        model.snpEffect.computeR(this->data, beta_old, r_old);
        VectorXf r_new = r_old;
        this->r_hist.row(i+1) = r_new;

        for (int j = 0; j < this-> data.numSNP; j++){
            

        }

        // update effectVar and residualVar


        // update Pi


        // compute hsq

    }

}

