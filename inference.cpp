#include "inference.hpp"
#include "model.hpp"

inferenceSBayesC::inferenceSBayesC(const std::string &binFilePath, 
                                   const std::string &phenoFilePath, 
                                   unsigned int num_iterations)
    : binFilePath(binFilePath), 
      phenoFilePath(phenoFilePath) 
{
    this->numberIterations = num_iterations; 
    this->data.readBinFullLD(this->binFilePath);
    this->data.readSummary(this->phenoFilePath);

    this->histMCMCSamples = MatrixXf::Zero(this->numberIterations, this->data.numSNP);
    this->r_hist = MatrixXf::Zero(this->numberIterations, this->data.numSNP);
}

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

    //std::cout << "Diagonal Matrix D:\n" << this->recon.D.block(0, 0, 10, 10) << std::endl;

    SBayesC model(this->data, this->numberIterations);   

    float mean = 0.0;
    float variance = 1.0;
    float pi = 0.5;

    // Call sampleFromPrior
    model.snpEffect.sampleFromPrior(this->data, model.currentState, this->histMCMCSamples, mean, variance, pi);
    model.snpEffect.initialR(this->data, this->histMCMCSamples, model.r_current, this->r_hist);

    // Print results
    //std::cout << "Sampled first 10 SNP effects:\n" << this->histMCMCSamples.block(0, 0, 10, 10) << std::endl;
    //std::cout << "Sampled first 10 values for r:\n" << this->r_hist.block(0, 0, 10, 10) << std::endl;
}

void inferenceSBayesC::runInference() {

}

