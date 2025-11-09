#ifndef RIDGE_LOGIT_H
#define RIDGE_LOGIT_H

#include <RcppEnsmallen.h>

// [[Rcpp::depends(RcppEnsmallen)]]

namespace logit {

class ridge_fn {
 private:
  const arma::mat X;
  const arma::vec Y;
  const double lambda;

 public:
  ridge_fn(const arma::mat& X, const arma::vec& Y, const double& lambda);

  double EvaluateWithGradient(const arma::mat& x, arma::mat& gradient);
};

arma::vec ridge_logit(const arma::mat& X, const arma::vec& Y, arma::vec beta_0,
                      const double& penalty);

}  // namespace logit
#endif
