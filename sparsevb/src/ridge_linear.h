#ifndef RIDGE_LINEAR_H
#define RIDGE_LINEAR_H

#include <RcppEnsmallen.h>

// [[Rcpp::depends(RcppEnsmallen)]]

namespace linear {

class ridge_fn {
 private:
  const arma::mat X;
  const arma::vec Y;
  const double lambda;

 public:
  ridge_fn(const arma::mat& X, const arma::vec& Y, const double& lambda);

  double EvaluateWithGradient(const arma::mat& x, arma::mat& gradient);
};

arma::vec ridge_linear(const arma::mat& X, const arma::vec& Y, arma::vec beta_0,
                       const double& penalty, const bool& exact);

}  // namespace linear
#endif
