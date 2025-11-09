#ifndef FIT_LOGIT_H
#define FIT_LOGIT_H

#include <RcppEnsmallen.h>

#include "laplace_obj_fn.h"
#include "misc_fn.h"
#include "ridge_logit.h"

// [[Rcpp::depends(RcppEnsmallen)]]

// [[Rcpp::export]]
Rcpp::List fit_logistic_laplace(const arma::mat& X, const arma::vec& Y,
                                arma::vec& mu, arma::vec& sigma,
                                arma::vec& gamma, const double& alpha,
                                const double& beta, const double& lambda,
                                arma::uvec& update_order,
                                const bool& prioritized_init,
                                const double& ridge_penalty,
                                const size_t& max_iter, const double& tol);

// [[Rcpp::export]]
Rcpp::List fit_logistic_gaussian(const arma::mat& X, const arma::vec& Y,
                                 arma::vec& mu, arma::vec& sigma,
                                 arma::vec& gamma, const double& alpha,
                                 const double& beta, const double& lambda,
                                 arma::uvec& update_order,
                                 const bool& prioritized_init,
                                 const double& ridge_penalty,
                                 const size_t& max_iter, const double& tol);

#endif
