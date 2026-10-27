#' Simulate Fractional Multinomial Logit Data
#'
#' Generate synthetic share data from two data-generating processes
#' satisfying the fractional multinomial logit conditional-mean
#' specification: a Dirichlet distribution and a multinomial
#' distribution normalised by the number of trials.
#'
#' @param N Positive integer. Number of observations to simulate.
#' @param J Integer greater than or equal to 2. Number of outcome
#'   categories (shares), including the baseline category.
#' @param K Positive integer. Number of columns in the design matrix,
#'   including the intercept. Thus, \eqn{K-1} explanatory variables
#'   are generated. Setting \eqn{K=1} produces an intercept-only model.
#' @param concentration Positive numeric scalar. Concentration parameter
#'   of the Dirichlet distribution. Larger values produce shares more
#'   tightly concentrated around their conditional means.
#' @param trials Positive integer. Number of multinomial trials per
#'   observation. The resulting shares are multiples of
#'   \eqn{1/\mathrm{trials}} and may equal exactly zero or one.
#' @param seed Integer random seed passed to \code{\link{set.seed}}
#'   for reproducibility.
#'
#' @details
#' The function generates \eqn{K-1} independent standard normal
#' covariates and an intercept, constructing an
#' \eqn{N \times K}{N x K} design matrix \eqn{X}.
#' The intercept occupies the first column.
#'
#' The true coefficient matrix \eqn{\beta} has dimensions
#' \eqn{J \times K}{J x K}. Its first row is normalised to zero,
#' making the first outcome the baseline category. All remaining
#' coefficients are independently drawn from the uniform
#' distribution on \eqn{[-1,1]}.
#'
#' Following \insertCite{mullahy2015multivariate}{fmlogit}, the
#' conditional expected shares are specified by the multinomial
#' logit transformation:
#' \deqn{
#'   p_{ij} = \mathbb{E}(y_{ij}\mid x_i) =
#'   \frac{\exp(x_i^\top \beta_j)}
#'        {\sum_{h=1}^{J}\exp(x_i^\top \beta_h)},
#'   \qquad \beta_1 \equiv 0.
#' }
#' The probabilities are computed using a numerically stable
#' softmax transformation. Consequently, all conditional means
#' lie in (0, 1) and sum to one for each observation,
#' up to floating-point precision.
#'
#' Two sets of observed shares are generated conditional on the
#' same covariates, coefficient matrix, and expected shares:
#'
#' \describe{
#'   \item{Dirichlet shares (\code{y.dirichlet}):}{
#'     Each observation is drawn from a Dirichlet distribution
#'     with parameter vector
#'     \eqn{(\phi p_{i1},\ldots,\phi p_{iJ})}{(phi p_i1, ..., phi p_iJ)},
#'     where \eqn{\phi}{phi} is \code{concentration}.
#'     The resulting shares are continuously distributed and
#'     strictly between zero and one in exact arithmetic.
#'   }
#'   \item{Multinomial shares (\code{y.multinom}):}{
#'     For each observation, a multinomial count vector is drawn
#'     with \code{trials} trials and category probabilities
#'     \eqn{(p_{i1},\ldots,p_{iJ})}{(p_i1, ..., p_iJ)}.
#'     Dividing the counts by the number of trials gives
#'     fractional shares that may take boundary values
#'     of exactly zero or one.
#'   }
#' }
#'
#' Both processes satisfy
#' \deqn{\mathbb{E}(y_{ij}\mid x_i)=p_{ij},}
#' conditional on the generated coefficients, although their
#' conditional variances and distributions differ.
#'
#' Specifically, their conditional variances are
#' \deqn{
#'   \mathrm{Var}(y_{ij}\mid x_i) =
#'   \begin{cases}
#'     p_{ij}(1-p_{ij})/(\phi+1), & \text{Dirichlet},\\
#'     p_{ij}(1-p_{ij})/m, & \text{Multinomial},
#'   \end{cases}
#' }
#' where \eqn{m} denotes \code{trials}.
#'
#' The function can therefore be used to investigate consistency,
#' finite-sample behaviour, and robustness of the fractional
#' multinomial logit quasi-maximum likelihood estimator under
#' alternative outcome distributions, including boundary values.
#'
#' The function calls \code{\link{set.seed}}, thereby resetting
#' the global random-number generator state.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{x}}{
#'     Numeric \eqn{N \times K}{N x K} design matrix.
#'     The first column is the intercept, followed by
#'     \eqn{K-1} standard normal covariates.
#'   }
#'   \item{\code{beta}}{
#'     Numeric \eqn{J \times K}{J x K} matrix of true
#'     coefficients, with outcomes in rows and regressors
#'     in columns. The first row contains only zeros;
#'     all other entries are drawn from \eqn{U(-1,1)}.
#'   }
#'   \item{\code{p}}{
#'     Numeric \eqn{N \times J}{N x J} matrix of conditional
#'     expected shares obtained from the multinomial
#'     logit transformation.
#'   }
#'   \item{\code{y.dirichlet}}{
#'     Numeric \eqn{N \times J}{N x J} matrix of Dirichlet
#'     shares, with columns named \code{share1}, ...,
#'     \code{shareJ}.
#'   }
#'   \item{\code{y.multinom}}{
#'     Numeric \eqn{N \times J}{N x J} matrix of multinomial
#'     count shares, with the same column names.
#'     Values are multiples of \eqn{1/\mathrm{trials}}.
#'   }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [fmlogit()]
#'
#' @examples
#' # Five outcomes and three design columns
#' sim <- simulatedShares(N = 200L, J = 5L, K = 3L)
#'
#' dim(sim$x)              # 200 x 3
#' dim(sim$beta)           # 5 x 3
#' dim(sim$p)              # 200 x 5
#' dim(sim$y.dirichlet)    # 200 x 5
#' dim(sim$y.multinom)     # 200 x 5
#'
#' # The baseline coefficient vector is zero
#' stopifnot(all(sim$beta[1L, ] == 0))
#'
#' # All shares sum to one
#' stopifnot(
#'   all(abs(rowSums(sim$p) - 1) < 1e-12),
#'   all(abs(rowSums(sim$y.dirichlet) - 1) < 1e-12),
#'   all(abs(rowSums(sim$y.multinom) - 1) < 1e-12)
#' )
#'
#' # Intercept-only model
#' sim0 <- simulatedShares(N = 100L, J = 4L, K = 1L)
#'
#' # Estimate fractional multinomial logit
#' \donttest{
#' fit <- fmlogit(
#'   y = sim$y.multinom,
#'   x = sim$x,
#'   print.level = 0
#' )
#' summary(fit)
#' }
#'
#' @export
simulatedShares <- function(N = 500L, J = 3L, K = 5L,
                            concentration = 10, trials = 10L, seed = 1L) {
  # Validate dimensions and simulation parameters
  positiveInteger <- function(value, name, minimum = 1L) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < minimum || value != floor(value) || value > .Machine$integer.max) {
      stop(sprintf("`%s` must be an integer >= %d.", name, minimum))
    }
    as.integer(value)
  }

  N <- positiveInteger(N, "N")
  J <- positiveInteger(J, "J", minimum = 2L)
  K <- positiveInteger(K, "K")
  trials <- positiveInteger(trials, "trials")

  if (!is.numeric(concentration) || length(concentration) != 1L ||
      !is.finite(concentration) || concentration <= 0) {
    stop("`concentration` must be a positive finite number.")
  }

  set.seed(seed)

  # Design matrix: intercept FIRST
  x <- matrix(1, nrow = N, ncol = K)

  if (K > 1L) x[, -1L] <- matrix(rnorm(N * (K - 1L)), nrow = N)

  colnames(x) <- c("(Intercept)", if (K > 1L) paste0("x", seq_len(K - 1L)))

  # True coefficients: J x K; the first alternative is the baseline
  beta <- matrix(0, nrow = J, ncol = K)
  beta[-1L, ] <- matrix(runif((J-1L)*K, min=-1, max=1), nrow = J-1L, ncol = K)
  rownames(beta) <- paste0("share", seq_len(J))
  colnames(beta) <- colnames(x)

  eta <- tcrossprod(x, beta) # Multinomial-logit conditional means
  eta <- sweep(eta, 1L, apply(eta, 1L, max), FUN = "-")  # Numerically stable softmax
  numerator <- exp(eta)
  p <- numerator / rowSums(numerator)

  colnames(p) <- rownames(beta)

  # Dirichlet DGP
  gamma.draws <- matrix(rgamma(N*J, shape = concentration * as.vector(p)), nrow = N, ncol = J)
  gamma.totals <- rowSums(gamma.draws)

  if (any(!is.finite(gamma.totals)) || any(gamma.totals <= 0))
    stop("Dirichlet simulation failed numerically; try a larger concentration or smaller dimensions.")

  y.dirichlet <- gamma.draws / gamma.totals

  # Multinomial DGP, allowing exact zeros and ones
  y.multinom <- t(vapply(seq_len(N), function(i) {
      as.vector(rmultinom(1L, size = trials, prob = p[i, ])) / trials
    }, numeric(J)))

  colnames(y.dirichlet) <- colnames(y.multinom) <- rownames(beta)

  list(x = x, beta = beta, p = p, y.dirichlet = y.dirichlet, y.multinom = y.multinom)
}
