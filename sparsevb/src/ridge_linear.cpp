#include "ridge_linear.h"

linear::ridge_fn::ridge_fn(const arma::mat& X, const arma::vec& Y,
                           const double& lambda)
    : X(X), Y(Y), lambda(lambda) {}

double linear::ridge_fn::EvaluateWithGradient(const arma::mat& x,
                                              arma::mat& gradient) {
  arma::vec YX_beta = Y - X * x;

  gradient =
      2 * lambda * x - 2 * X.t() * YX_beta / static_cast<double>(Y.n_elem);

  return lambda * arma::accu(x % x) -
         std::pow(arma::norm(YX_beta), 2) / static_cast<double>(Y.n_elem);
}

arma::vec linear::ridge_linear(const arma::mat& X, const arma::vec& Y,
                               arma::vec beta_0, const double& penalty,
                               const bool& exact) {
  if (exact) {
    return arma::inv_sympd(X.t() * X +
                           penalty * arma::eye(beta_0.n_elem, beta_0.n_elem)) *
           X.t() * Y;
  } else {
    ens::L_BFGS optim;
    ridge_fn f(X, Y, penalty);

    // note that the result is written into beta_0
    optim.Optimize(f, beta_0);

    return beta_0;
  }
}
