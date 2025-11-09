#include "ridge_logit.h"

logit::ridge_fn::ridge_fn(const arma::mat& X, const arma::vec& Y, const double& lambda)
    : X(X), Y(Y), lambda(lambda) {}

double logit::ridge_fn::EvaluateWithGradient(const arma::mat& x, arma::mat& gradient) {
  arma::vec X_beta = X * x;
  arma::vec exp_term = arma::exp(X_beta);
  arma::vec logit_term = 1 / (1 + exp_term);

  gradient = 2 * lambda * x -
             ((Y + logit_term).t() * X).t() / static_cast<double>(Y.n_elem);

  return lambda * arma::accu(x % x) -
         (arma::dot(Y, X_beta) + arma::accu(arma::log1p(exp_term))) /
             static_cast<double>(Y.n_elem);
}

arma::vec logit::ridge_logit(const arma::mat& X, const arma::vec& Y, arma::vec beta_0,
                      const double& penalty) {
  ens::L_BFGS optim;
  ridge_fn f(X, Y, penalty);

  // note that the result is written into beta_0
  optim.Optimize(f, beta_0);

  return beta_0;
}
