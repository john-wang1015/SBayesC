#include "fit_logit.h"

Rcpp::List fit_logistic_laplace(const arma::mat& X, const arma::vec& Y,
                                arma::vec& mu, arma::vec& sigma,
                                arma::vec& gamma, const double& alpha,
                                const double& beta, const double& lambda,
                                arma::uvec& update_order,
                                const bool& prioritized_init,
                                const double& ridge_penalty,
                                const size_t& max_iter, const double& tol) {
  // initalize entropy loss function
  arma::vec old_entr = entropy(gamma);

  // initialize L-BFGS optimizer from ensmallen
  ens::L_BFGS optim;

  // perform prioritized init
  if (prioritized_init) {
    // mu will be overwritten with the result of the ridge solver
    mu = logit::ridge_logit(X, Y, mu, ridge_penalty);

    // generate prioritized update order
    update_order = arma::sort_index(arma::abs(mu), "descending");
  }

  // pre-process update parameters
  arma::vec approx_mean = gamma % mu;
  arma::vec X_appm = X * approx_mean;
  arma::vec X_sq = arma::square(X) * arma::square(approx_mean);
  arma::vec X_secm =
      arma::square(X) * (gamma % (arma::square(mu) + arma::square(sigma)));
  double const_lodds = std::log(alpha) - std::log(beta) + std::log(lambda);
  arma::rowvec YX_vec = (Y - 0.5).t() * X;

  // iteration loop
  for (size_t i = 0; i < max_iter; ++i) {
    // implements equation (32) of the paper
    arma::vec eta = arma::sqrt(X_secm + arma::square(X_appm) - X_sq);

    // pre-processing per iteration
    arma::vec eta_hyp = hyperbolic_transform(eta);
    arma::rowvec coef_sq = eta_hyp.t() * arma::square(X);

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

      // delete the j-th column from running sums
      X_appm -= approx_mean(j) * X.col(j);
      X_sq -= std::pow(approx_mean(j), 2) * arma::square(X.col(j));
      X_secm -= gamma(j) * (std::pow(mu(j), 2) + std::pow(sigma(j), 2)) *
                arma::square(X.col(j));

      // implements equation (11) of the paper
      laplace_obj_fn f(coef_sq(j),
                       2 * arma::dot(eta_hyp % X.col(j), X_appm) - YX_vec(j),
                       lambda);

      // optimize and save function value
      double opt = optim.Optimize(f, x);
      mu(j) = x(0);
      sigma(j) = std::abs(x(1));

      // implements equation (12) of the paper
      gamma(j) = sigmoid(const_lodds - opt);

      // add j-th column with updated values
      approx_mean(j) = gamma(j) * mu(j);
      X_appm += approx_mean(j) * X.col(j);
      X_sq += std::pow(approx_mean(j), 2) * arma::square(X.col(j));
      X_secm += gamma(j) * (std::pow(mu(j), 2) + std::pow(sigma(j), 2)) *
                arma::square(X.col(j));
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

Rcpp::List fit_logistic_gaussian(const arma::mat& X, const arma::vec& Y,
                                 arma::vec& mu, arma::vec& sigma,
                                 arma::vec& gamma, const double& alpha,
                                 const double& beta, const double& lambda,
                                 arma::uvec& update_order,
                                 const bool& prioritized_init,
                                 const double& ridge_penalty,
                                 const size_t& max_iter, const double& tol) {
  // initalize entropy loss function
  arma::vec old_entr = entropy(gamma);

  // initialize L-BFGS optimizer from ensmallen
  ens::L_BFGS optim;

  // perform prioritized init
  if (prioritized_init) {
    // compute ridge regression estimator of mu
    mu = logit::ridge_logit(X, Y, mu, ridge_penalty);

    // generate prioritized update order
    update_order = arma::sort_index(arma::abs(mu), "descending");
  }

  // pre-process update parameters
  arma::vec approx_mean = gamma % mu;
  arma::vec X_appm = X * approx_mean;
  arma::vec X_sq = arma::square(X) * arma::square(approx_mean);
  arma::vec X_secm =
      arma::square(X) * (gamma % (arma::square(mu) + arma::square(sigma)));
  double const_lodds = std::log(alpha) - std::log(beta) - std::log(lambda);
  arma::rowvec YX_vec = (Y - 0.5).t() * X;

  // iteration loop
  for (size_t i = 0; i < max_iter; ++i) {
    // implements equation (32) of the paper
    arma::vec eta = arma::sqrt(X_secm + arma::square(X_appm) - X_sq);

    // pre-processing per iteration
    arma::vec eta_hyp = hyperbolic_transform(eta);
    arma::rowvec coef_sq = eta_hyp.t() * arma::square(X);

    // implements equation (25) of Carbonetto et al.
    sigma = 1 / arma::sqrt(2 * coef_sq.t() + 1 / std::pow(lambda, 2));

    // update each coordinate
    for (arma::uword k = 0; k < mu.n_elem; ++k) {
      // check if interrupt signal was sent from R
      Rcpp::checkUserInterrupt();

      // the current update dimension
      arma::uword j = update_order(k);

      // delete the j-th column from running sums
      X_appm -= approx_mean(j) * X.col(j);
      X_sq -= std::pow(approx_mean(j), 2) * arma::square(X.col(j));
      X_secm -= gamma(j) * (std::pow(mu(j), 2) + std::pow(sigma(j), 2)) *
                arma::square(X.col(j));

      // implements equation (26) of Carbonetto et al.
      mu(j) = std::pow(sigma(j), 2) *
              (YX_vec(j) - 2 * arma::dot(eta_hyp % X.col(j), X_appm));

      // implements equation (27) of Carbonetto et al.
      gamma(j) = sigmoid(const_lodds + std::log(sigma(j)) +
                         0.5 * std::pow(mu(j) / sigma(j), 2));

      // add j-th column with updated values
      approx_mean(j) = gamma(j) * mu(j);
      X_appm += approx_mean(j) * X.col(j);
      X_sq += std::pow(approx_mean(j), 2) * arma::square(X.col(j));
      X_secm += gamma(j) * (std::pow(mu(j), 2) + std::pow(sigma(j), 2)) *
                arma::square(X.col(j));
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
