#include "fit_linear.h"

Rcpp::List fit_linear_laplace(const arma::mat& X, const arma::vec& Y,
                              arma::vec& mu, arma::vec& sigma, arma::vec& gamma,
                              const double& alpha, const double& beta,
                              const double& lambda, arma::uvec& update_order,
                              const bool& prioritized_init,
                              const double& ridge_penalty,
                              const size_t& max_iter, const double& tol,
                              const bool& exact_math) {
  // initalize entropy loss function
  arma::vec old_entr = entropy(gamma);

  // initialize L-BFGS optimizer from ensmallen
  ens::L_BFGS optim;

  // perform prioritized init
  if (prioritized_init) {
    // compute ridge regression estimator of mu
    mu = linear::ridge_linear(X, Y, mu, ridge_penalty, exact_math);

    // generate prioritized update order
    update_order = arma::sort_index(arma::abs(mu), "descending");
  }

  // pre-process update parameters
  double const_lodds = std::log(alpha) - std::log(beta) + 0.5 +
                       0.5 * std::log(M_PI) + std::log(lambda) - 0.5 * M_LN2;
  arma::rowvec YX_vec = Y.t() * X;
  arma::vec half_diag(mu.n_elem);
  for (arma::uword i = 0; i < half_diag.n_elem; ++i) {
    half_diag(i) = 0.5 * std::pow(arma::norm(X.col(i)), 2);
  }
  arma::vec approx_mean = gamma % mu;
  arma::vec X_appm = X * approx_mean;

  // iteration loop
  for (size_t i = 0; i < max_iter; ++i) {
    // update each coordinate
    for (arma::uword k = 0; k < mu.n_elem; ++k) {
      // check if interrupt signal was sent from R
      Rcpp::checkUserInterrupt();

      // the current update dimension
      arma::uword j = update_order(k);

      // start optimization at previous value
      arma::mat x(2, 1);
      x(0) = mu(j);
      x(1) = sigma(j);

      // delete the j-th column from X * approx_mean
      X_appm -= approx_mean(j) * X.col(j);

      // implements equation (16) of the paper
      laplace_obj_fn f(half_diag(j), arma::dot(X.col(j), X_appm) - YX_vec(j),
                       lambda);

      // optimize and save function value
      double opt = optim.Optimize(f, x);
      mu(j) = x(0);
      sigma(j) = std::abs(x(1));

      // implements equation (17) of the paper
      gamma(j) = sigmoid(const_lodds - opt);

      // add j-th column with updated values
      approx_mean(j) = gamma(j) * mu(j);
      X_appm += approx_mean(j) * X.col(j);
    }

    // check for convergence
    arma::vec new_entr = entropy(gamma);
    if (arma::norm(new_entr - old_entr, "inf") <= tol) {
      break;
    } else {
      old_entr = new_entr;
    }
  }

  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("sigma") = sigma,
                            Rcpp::Named("gamma") = gamma);
}

Rcpp::List fit_linear_gaussian(
    const arma::mat& X, const arma::vec& Y, arma::vec& mu, arma::vec& sigma,
    arma::vec& gamma, const double& alpha, const double& beta,
    const double& lambda, arma::uvec& update_order,
    const bool& prioritized_init, const double& ridge_penalty,
    const size_t& max_iter, const double& tol, const bool& exact_math) {
  // initalize entropy loss function
  arma::vec old_entr = entropy(gamma);

  // initialize L-BFGS optimizer from ensmallen
  ens::L_BFGS optim;

  // perform prioritized init
  if (prioritized_init) {
    // compute ridge regression estimator of mu
    mu = linear::ridge_linear(X, Y, mu, ridge_penalty, exact_math);

    // generate prioritized update order
    update_order = arma::sort_index(arma::abs(mu), "descending");
  }

  // pre-process update parameters
  double const_lodds = std::log(alpha) - std::log(beta) - std::log(lambda);
  arma::rowvec YX_vec = Y.t() * X;
  for (arma::uword i = 0; i < sigma.n_elem; ++i) {
    // implements equation (8) of Carbonetto et al.
    sigma(i) = 1 / std::sqrt(std::pow(arma::norm(X.col(i)), 2) +
                             1 / std::pow(lambda, 2));
  }
  arma::vec approx_mean = gamma % mu;
  arma::vec X_appm = X * approx_mean;

  // iteration loop
  for (size_t i = 0; i < max_iter; ++i) {
    // update each coordinate
    for (arma::uword k = 0; k < mu.n_elem; ++k) {
      // check if interrupt signal was sent from R
      Rcpp::checkUserInterrupt();

      // the current update dimension
      arma::uword j = update_order(k);

      // delete the j-th column from X * approx_mean
      X_appm -= approx_mean(j) * X.col(j);

      // implements equation (9) of Carbonetto et al.
      mu(j) = std::pow(sigma(j), 2) * (YX_vec(j) - arma::dot(X.col(j), X_appm));

      // implements equation (10) of Carbonetto et al.
      gamma(j) = sigmoid(const_lodds + std::log(sigma(j)) +
                         0.5 * std::pow(mu(j) / sigma(j), 2));

      // add j-th column with updated values
      approx_mean(j) = gamma(j) * mu(j);
      X_appm += approx_mean(j) * X.col(j);
    }

    // check for convergence
    arma::vec new_entr = entropy(gamma);
    if (arma::norm(new_entr - old_entr, "inf") <= tol) {
      break;
    } else {
      old_entr = new_entr;
    }
  }

  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("sigma") = sigma,
                            Rcpp::Named("gamma") = gamma);
}
