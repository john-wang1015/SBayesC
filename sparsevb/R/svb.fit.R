#' Fit Approximate Posteriors to Sparse Linear and Logistic Models
#'
#' @description Main function of the \code{\link{sparsevb}} package. Computes
#'   mean-field posterior approximations for both linear and logistic regression
#'   models, including variable selection via sparsity-inducing spike and slab
#'   priors.
#'
#' @param X A numeric design matrix, each row of which represents a data point.
#' @param Y An \code{nrow(X)}-dimensional response vector, numeric if
#'   \code{family = "linear"} and binary if \code{family = "logistic"}.
#' @param family A character string selecting the regression model, either
#'   \code{"linear"} or \code{"logistic"}. (\emph{default:} \code{"linear"})
#' @param slab A character string specifying the prior slab density, either
#'   \code{"laplace"} or \code{"gaussian"}. (\emph{default:} \code{"laplace"})
#' @param mu An \code{ncol(X)}-dimensional numeric vector, serving as initial
#'   guess for the variational means. (\emph{default:} \code{rep(0, ncol(X))})
#' @param sigma A positive \code{ncol(X)}-dimensional numeric vector, serving as
#'   initial guess for the variational standard deviations. (\emph{default:}
#'   \code{rep(1, ncol(X))})
#' @param gamma An \code{ncol(X)}-dimensional vector of probabilities, serving
#'   as initial guess for the variational inclusion probabilities.
#'   (\emph{default:} \code{rep(alpha/(alpha+beta), ncol(X))})
#' @param alpha A positive numeric value, used by the Beta-hyperprior.
#'   (\emph{default:} \code{1.0})
#' @param beta A positive numeric value, used by the Beta-hyperprior.
#'   (\emph{default:} \code{ncol(X)})
#' @param lambda A numeric value, controlling the scale parameter of the prior
#'   slab density. Used as the inverted scale parameter when \code{prior =
#'   "laplace"} (default) and the standard deviation if \code{prior =
#'   "gaussian"}. (\emph{default:} \code{1.0})
#' @param update_order A permutation of \code{1:ncol(X)}, giving the update
#'   order of the coordinate-ascent algorithm. Setting this parameter when
#'   \code{prioritized_init = TRUE} has no effect; it will be overwritten.
#'   (\emph{default:} \code{1:ncol(X)})
#' @param prioritized_init A Boolean value, controlling whether the ridge
#'   regression estimator for \code{mu} should be computed during
#'   initialization. When \code{TRUE}, the argument \code{mu} serves as initial
#'   guess for the ridge regression estimator and \code{update_order} is
#'   overwritten by ranking the elements the estimator according to magnitude.
#'   (\emph{default:} \code{TRUE})
#' @param ridge_penalty A positive numerical value, controlling the importance
#'   of the penalty term when computing the ridge regression estimator.
#'   (\emph{default:} \code{0.1})
#' @param exact_math A Boolean variable, controlling if the linear ridge
#'   regression estimator should be computed in closed form or iteratively. Has
#'   no effect when \code{family = "logistic"}. (\emph{default:} \code{FALSE})
#' @param rescale A Boolean variable, controlling if \code{X} and \code{Y}
#'   should be rescaled by the estimated variance of the underlying noise. Has
#'   no effect when \code{family = "logistic"}. (\emph{default:} \code{TRUE})
#' @param max_iter A positive integer, controlling the maximum number of
#'   iterations for the variational update loop. (\emph{default:} \code{1000})
#' @param tol A positive numerical value, controlling the termination criterion
#'   for maximum absolute differences between binary entropies of successive
#'   iterates. (\emph{default:} \code{10e-6})
#'
#' @return The approximate mean-field posterior, given as a named list
#'   containing numeric vectors \code{"mu"}, \code{"sigma"}, and \code{"gamma"}.
#'   In mathematical terms, \deqn{\theta_j\mid \mu_j, \sigma_j, \gamma_j
#'   \sim_{ind.} \gamma_j N(\mu_j, \sigma^2) + (1-\gamma_j) \delta_0.}
#'
#' @examples
#' \dontrun{
#'
#' ### Simulate a linear regression problem of size n times p, with sparsity level s ###
#'
#' n <- 2500
#' p <- 5000
#' s <- 25
#'
#' ### Generate toy data ###
#'
#' X <- matrix(rnorm(n*p), n, p) #standard Gaussian design matrix
#'
#' theta <- numeric(p)
#' theta[sample.int(p, s)] <- runif(s, -3, 3) #sample non-zero coefficients in random locations
#'
#' pos_TR <- as.numeric(theta != 0) #true positives
#'
#' Y <- X %*% theta + rnorm(n) #add standard Gaussian noise
#'
#' ### Run the algorithm in linear mode with Laplace prior and prioritized initialization ###
#'
#' test <- svb.fit(X, Y, family = "linear")
#'
#' posterior_mean <- test$mu * test$gamma #approximate posterior mean
#'
#' pos <- as.numeric(test$gamma > 0.5) #significant coefficients
#'
#' ### Assess the quality of the posterior estimates ###
#'
#' TPR <- sum(pos[which(pos_TR == 1)])/sum(pos_TR) #True positive rate
#'
#' FDR <- sum(pos[which(pos_TR != 1)])/max(sum(pos), 1) #False discovery rate
#'
#' L2 <- sqrt(sum((posterior_mean - theta)^2)) #L_2-error
#'
#' MSPE <- sqrt(sum((X %*% posterior_mean - Y)^2)/n) #Mean squared prediction error
#' }
#'
#' @details Suppose \eqn{\theta} is the \eqn{p}-dimensional true parameter. The
#'   spike-and-slab prior for \eqn{\theta} may be represented by the
#'   hierarchical scheme \deqn{w \sim \mathrm{Beta}(\alpha, \beta)} \deqn{z_j
#'   \mid w \sim_{i.i.d.} \mathrm{Bernoulli}(w)} \deqn{\theta_j\mid z_j
#'   \sim_{ind.} (1-z_j)\delta_0 + z_j g.} As usual, \eqn{\delta_0} represents
#'   the Dirac measure at \eqn{0}. The slab \eqn{g} may be taken either as a
#'   \eqn{\mathrm{Laplace}(0,\lambda^{-1})} or \eqn{N(0,\lambda^2)} density.
#'
#' @export
svb.fit <- function(X,
                    Y,
                    family = c("linear", "logistic"),
                    slab = c("laplace", "gaussian"),
                    mu,
                    sigma,
                    gamma,
                    alpha = 1,
                    beta,
                    lambda = 1,
                    update_order,
                    prioritized_init = TRUE,
                    exact_math = FALSE,
                    rescale = TRUE,
                    ridge_penalty = 0.1,
                    max_iter = 1000,
                    tol = 10e-6) {
    #extract problem dimensions
    n = dim(X)[1]
    p = dim(X)[2]
    
    #initialize missing arguments
    if (missing(mu)) {
        mu = rep(0, p)
    }
    if (missing(sigma)) {
        sigma = rep(1, p)
    }
    if (missing(update_order)) {
        update_order = 1:p - 1
    }
    if (missing(beta)) {
        beta = p
    }
    if (missing(gamma)) {
        gamma = rep(alpha / (alpha + beta), p)
    }
    
    #match internal function call and generate list of arguments
    fn = paste("fit", match.arg(family), match.arg(slab), sep = '_')
    arg = as.list(environment())
    arg = within(arg, rm(family, slab, n, p, fn, rescale))
    if (match.arg(family) == "logistic") {
        arg = within(arg, rm(exact_math))
    }
    
    if(rescale && match.arg(family) == "linear"){
        var = estimateSigma(X, Y)$sigmahat
        X = X/var
        Y = Y/var
    }
    
    #perform chosen computation
    approximate_posterior = do.call(fn, arg)
    
    #convert results to R-style vectors since RcppArmadillo returns in matrix form
    return(lapply(approximate_posterior, as.numeric))
}
