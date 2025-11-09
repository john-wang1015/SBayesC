#ifndef FIT_LINEAR_H
#define FIT_LINEAR_H

#include <RcppEnsmallen.h>

#include "laplace_obj_fn.h"
#include "misc_fn.h"
#include "ridge_linear.h"

// [[Rcpp::depends(RcppEnsmallen)]]

// [[Rcpp::export]]
Rcpp::List fit_linear_laplace(const arma::mat& X, const arma::vec& Y,
                              arma::vec& mu, arma::vec& sigma, arma::vec& gamma,
                              const double& alpha, const double& beta,
                              const double& lambda, arma::uvec& update_order,
                              const bool& prioritized_init,
                              const double& ridge_penalty,
                              const size_t& max_iter, const double& tol,
                              const bool& exact_math);

// [[Rcpp::export]]
Rcpp::List fit_linear_gaussian(
    const arma::mat& X, const arma::vec& Y, arma::vec& mu, arma::vec& sigma,
    arma::vec& gamma, const double& alpha, const double& beta,
    const double& lambda, arma::uvec& update_order,
    const bool& prioritized_init, const double& ridge_penalty,
    const size_t& max_iter, const double& tol, const bool& exact_math);

#endif
