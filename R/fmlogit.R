#' Find good initial values for fmlogit
#' 
#' Run pairwise fractional logits of \eqn{y_j / (y_1 + y_j)}{y_j / (y_1 + y_j)}.
#' 
#' @details
#' The reasonable initial value comes from `J-1` separate weighted binary-logit
#' regressions, each comparing one alternative with the baseline alternative.
#' Since the normalisation is such that the first `K` elements are zeros (and
#' do not enter the optimisation), we note that
#' \eqn{E(y_{ij}\mid x_i) = p_{ij}}{E(y_ij | x_i) = p_ij} and
#' \eqn{p_{ij}/p_{i1} = \exp (x_i' \beta_j)}{p_ij / p_i1 = exp (x_i' beta_j)}.
#' Now define \eqn{\pi_{ij} = \frac{p_{ij}}{p_{i1} + p_{ij}}}{pi_ij = p_ij / (p_i1 + p_ij)}.
#' Substituting the multinomial-logit probabilities gives
#' \eqn{\pi_{ij} = \frac{\exp(x_i' \beta_j)}{1 + \exp(x_i' \beta_j)} = 
#' \Lambda(x_i' \beta_j)}{pi_ij = exp(x_i' beta_j) / (1 + exp(x_i' beta_j)) =
#' Lambda(x_i' beta_j)}.
#' This is a pairwise application of the fractional quasi-likelihood logic in
#' Papke and Wooldridge (1996) and Mullahy (2015).
#' 
#' 
#' @export
initEst <- function(y, x, weights = NULL) {
  y <- as.matrix(y)
  x <- as.matrix(x)

  j <- ncol(y)
  K <- ncol(x)
  n <- NROW(y)
  
  if (is.null(colnames(y))) colnames(y) <- paste0("share", 1:j)
  if (is.null(colnames(x))) colnames(x) <- paste0("X", 1:K)
  
  if (is.null(weights)) weights <- rep(1, n)
  if (!is.numeric(weights) ||  length(weights) != n0 || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite, non-negative, and match NROW(y).")
  }

  # Estimate each non-baseline alternative against alternative 1
  b <- lapply(2:j, function(h) {
    # Total observed share of baseline and alternative h
    w <- y[, 1] + y[, h]
    # Observations where both shares are zero contribute nothing
    keep <- w > 0
    # Pairwise fractional response
    dat <- data.frame(z = y[keep, h] / w[keep], x[keep, , drop = FALSE])
    
    fit <- glm(z ~ ., data = dat, weights = w[keep] * weights[keep],
               family = quasibinomial(link = "logit"))
    coef(fit)
  })
  names(b) <- colnames(y)[2:j]
  b <- unlist(b)

  which.bad <- which(!is.finite(b))
  if (any(which.bad)) {
    b[1:length(b)] <- 0
    warning("A pairwise GLM did not return finite coefficients; using zeros instead.")
  }
  
  return(b)
}


# Estimation of Fractional Multinomial Logit Models
#
# fmlogit allows you to estimate fractional multinomial logit models with
# Quasi-MLE. It also provides a method for calculating marginal and discrete
# effects, predictions, and aggregate willingness-to-pay calculations.
#
# The main function to call is \code{fmlogit}, which will estimate the model.
# Using \code{effects} to estimate marginal or discrete effects will also be
# very helpful.

#' Estimate Fractional Multinomial Logit Models
#'
#' Used to estimate fractional multinomial logit models using quasi-maximum
#' likelihood estimations following Papke and Wooldridge(1996).
#' @param y the dependent variable (N*J). Can be a matrix or a named data frame.
#'   The first column of the matrix is automatically treated as the baseline.
#' @param x independent variable (N*K). Can be a matrix or a named data frame.
#'   An intercept term is automatically added.
#' @param weights Observation (survey) weights denoting the relative importance
#'   of observations. Defaults to equal weights (a vector of ones).
#' @param beta0 Initial value for the coefficient vector used in optimization.
#'   Default to a vector of smart values built from individual `J-1` estimated
#'   equations (via `glm()`), concatenated one after another.
#' @param method Method of optimization. Passed to `maxLik()`.
#'   Choose from "NR", "BFGS", "CG", "BHHH", "SANN",or "NM". Default to "CG",
#'   the conjugate gradients method. See Details.
#' @param iterlim Maximum number of iterations.
#' @param abstol Absolute convergence tolerance.
#' @param print.level Integer: optimiser verbosity. `0` silences the output,
#'   `1` prints the initial and final details, `2` prints each iteration's
#'   output.
#' @param ... additional parameters passed to \code{maxLik()}
#'
#' @return The function returns an object of class "fmlogit". Use
#'   \code{effects}, \code{predict}, \code{residual}, \code{fitted} to extract
#'   various useful features of the value returned by \code{fmlogit}.
#' @return An object of class "fmlogit" contains the following components:
#' @return \code{estimates}   A list of matrices containing parameter estimates,
#'   standard errors, and hypothesis testing results.
#' @return \code{baseline}    The baseline choice
#' @return \code{likelihood}  The likelihood value
#' @return \code{conv_code}   Convergence diagnostics code.
#' @return \code{convergence} Convergence messages.
#' @return \code{count}       Provides dataset information
#' @return \code{y}           The dependent variable data frame.
#' @return \code{x}           The independent variable data frame. Augmented by
#'   factor dummy transformation, constant term added.
#' @return \code{rowNo}       A vector of row numbers from the original `x` and
#'   `y` that is used for estimation.
#' @return \code{coefficient} Matrix of estimated coefficients. Augmented with
#'   the baseline coefficient (which is a vector of zeros).
#' @return \code{vcov} A list of matrices containing the robust
#'   variance-covariance matrix for each choice variable.
#' @details The fractional multinomial model is the expansion of the multinomial
#' logit to fractional responses. Unlike standard multinomial logit models,
#' which only considers 0-1 respones, fractional multinomial model considers the
#' case where the response variable is fractions that sums up to one. Examples
#' of these type of data are, percentages of budget spent in education, defense,
#' public health; fractions of a population that have middle school, high
#' school, college, or post college education, etc.
#'
#' This function follows Papke and Wooldridge(1996)'s paper, in which they
#' proposed a quasi-maximum likelihood estimator for fractional response data.
#' The likelihood function used here is a standard multinomial likelihood
#' function, see \url{http://maartenbuis.nl/software/likelihoodFmlogit.pdf} for
#' the likelihood used here. Robust standard errors are provided following Papke
#' and Wooldridge (1996), in which they proposed an asymptotically consistent
#' estimator of variance.
#'
#' Maximization is done by calling \code{\link{maxLik}}. `maxLik` is a wrapper
#' function for different maximization methods in R. This include most methods
#' provided by \code{\link{maxLik}},  but also other methods such as BHHH
#' (Berndt--Hall--Hall--Hausman).
#' 
#'
#' MLE convergence can be a problem in R, especially if the data size is large
#' with many explanatory variables. It is recommended to call CG (Conjugate
#' Gradients) or BHHH (Berndt-Hall-Hall-Hausman). Conjugate gradients method is
#' usually faster, but could lead to non-convergence under certain scenarios.
#' BHHH is slower, but has better convergence performance.
#'
#'
#' @examples
#' data <- spending
#' x <- data[,2:5]
#' y <- data[,6:11]
#' results1 <- fmlogit(y,x)
#' print(results1)
#' summary(results1
#' )
#' @references Papke, L. E. and Wooldridge, J. M. (1996), Econometric methods
#'   for fractional response variables with an application to 401(k) plan
#'   participation rates. J. Appl. Econ., 11: 619-632.
#' @export fmlogit


fmlogit <- function(y, x, weights = NULL, beta0 = NULL, method = "CG",
                    iterlim = 5e+05, abstol = 1e-05, print.level = 1, ...) {
  start.time <- proc.time()

  n0 <- NROW(y)
  if (is.null(weights)) weights <- rep(1, n0)
  if (!is.numeric(weights) ||  length(weights) != n0 || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite, non-negative, and match NROW(y).")
  }

  xclass <- sapply(x, class)
  xfac <- which(xclass %in% c("factor", "character"))
  if (length(xfac) > 0) {
    xfacnames <- colnames(x)[xfac]
    strformFac <- paste(xfacnames, collapse = "+")
    xdum <- model.matrix(as.formula(paste("~", strformFac,
                                          sep = "")), data = x)[, -1]
    x <- cbind(x, xdum)
    x <- x[, -xfac]
  }
  xnames <- colnames(x)
  ynames <- colnames(y)
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- dim(x)[1]  # Number of observations
  j <- dim(y)[2]  # Columns of Y, i.e. multiple choices
  k <- dim(x)[2]  # Number of explanatory variables
  xy <- cbind(x, y)
  xy <- na.omit(xy)
  row.remain <- setdiff(1:n, attr(xy, "na.action"))
  x <- xy[, 1:k]
  y <- xy[, (k + 1):(k + j)]
  n <- dim(y)[1]
  remove(xy)

  weights <- weights[row.remain]
  if (!any(weights > 0)) stop("At least one retained observation must have positive weight")

  # adding in the constant term
  if (k == 1) {
    # check if the input x is constant
    if (length(unique(x)) == 1) { # x is constant
      xnames <- "(Intercept)"
      x <- as.matrix(as.numeric(x), nrow = 1)
      colnames(x) <- xnames
      x <- as.matrix(x)
      k <- 0
    } else { # one single variable of input
      xnames <- "x1"
      x <- as.matrix(x)
      k <- dim(x)[2]
      x <- cbind(x, rep(1, n))
      xnames <- c(xnames, "(Intercept)")
      colnames(x) <- xnames
    }
  } else { # normal cases
    x <- x[, apply(x, 2, function(x) length(unique(x)) != 1)]
    xnames <- colnames(x)
    k <- dim(x)[2]
    x <- cbind(x, rep(1, n))
    xnames <- c(xnames, "(Intercept)")
    colnames(x) <- xnames
  }


  collinear <- unique(unlist(testcols(x)))
  while (length(collinear) > 0) {
    if (qr(x)$rank == dim(x)[2])
      print("The model may suffer from multicollinearity problems.")
    break
    if ((k + 1) %in% collinear)
      collinear <- collinear[-length(collinear)]
    x <- x[, -collinear[length(collinear)]]
    xnames <- colnames(x)
    k <- k - 1
    collinear <- unique(unlist(testcols(x)))
  }

  QMLE <- function(betas) {
    betas <- matrix(betas, nrow = j - 1, byrow = TRUE)
    betamat <- rbind(rep(0, k + 1), betas)
    llf <- rep(0, n)
    for (i in 1:j) {
      L <- y[, i] * ((x %*% betamat[i, ]) - log(rowSums(exp(x %*% t(betamat)))))
      llf <- llf + L
    }
    return(weights * llf)
  }
  if (is.null(beta0)) beta0 <- initEst(y = y, x = x, weights = weights)
  if (length(beta0) != (k + 1) * (j - 1)) {
    beta0 <- initEst(y = y, x = x, weights = weights)
    warning("Wrong length of beta0 given. Using smart initial values instead.")
  }
  opt <- maxLik::maxLik(QMLE, start = beta0, method = method,
                control = list(iterlim = iterlim, tol = abstol, print.level = print.level), ...)
  betamat <- matrix(opt$estimate, ncol = k + 1, byrow = T)
  betamat_aug <- rbind(rep(0, k + 1), betamat)
  colnames(betamat_aug) <- xnames
  rownames(betamat_aug) <- ynames
  sigmat <- matrix(nrow = j - 1, ncol = k + 1)
  vcov <- list()

  ## Robust VCOV with weights

  for(i in 1:j) {
      # start calculation
      sum_expxb <- rowSums(exp(x %*% t(betamat_aug))) # sum of the exp(x'b)s
      expxb <- exp(x %*% betamat_aug[i,]) # individual exp(x'b)
      G <- expxb / sum_expxb # exp(x'bj) / sum^J(exp(x'bj))
      g <- (expxb * sum_expxb - expxb^2) / sum_expxb^2 # derivative of the logit function

      # Here the diagonal of A is the 'standard' standard error
      # hat(A) = sum hat(gi)^2 * xi'xi / hat(Gi)(1-hat(Gi))
      # or, xtilde = x * sqrt(g^2/G(1-G)), A = xtilde'xtilde
      x_a <- x * as.vector(sqrt(g^2/(G*(1-G))))
      A <- t(x_a) %*% x_a

      # robust standard error, again following PW(1996)
      mu <- y[,i] - G
      x_b <- x * as.vector(mu * g / G / (1-G))
      B <- t(x_b) %*% x_b
      Var_b <- solve(A) %*% B %*% solve(A)
      std_b <- sqrt(diag(Var_b))
      # std_b= sqrt(diag(solve(A))) is the "unrobust" standard error.
      vcov[[i]] <- Var_b
      if (i > 1) sigmat[i-1,] <- std_b
    }

  ###end of insert

  listmat <- list()
  for (i in 1:(j - 1)) {
    tabout <- matrix(ncol = 4, nrow = k + 1)
    tabout[, 1:2] <- t(rbind(betamat[i, ], sigmat[i, ]))
    tabout[, 3] <- tabout[, 1]/tabout[, 2]
    tabout[, 4] <- 2 * (1 - pnorm(abs(tabout[, 3])))
    colnames(tabout) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    if (length(xnames) > 0) rownames(tabout) <- xnames
    listmat[[i]] <- tabout
  }
  if (length(ynames) > 0) names(listmat) <- ynames[2:j]
  outlist <- list()
  outlist$raw.estimates <- opt$estimate
  outlist$estimates <- listmat
  outlist$baseline <- ynames[1]
  outlist$likelihood <- opt$maximum
  outlist$conv_code <- as.integer(opt$code)
  outlist$convergence <- paste(opt$type,
                               paste(as.character(opt$iterations), "iterations"),
                               opt$message, sep = ",")
  outlist$count <- c(Obs = n, Explanatories = k, Choices = j)
  outlist$y <- y
  outlist$x <- x
  outlist$weights <- weights
  outlist$rowNo <- row.remain
  outlist$coefficient <- betamat_aug
  names(vcov) <- ynames
  outlist$vcov <- vcov

  print(paste("Fractional logit model estimation completed. Time:",
              round(proc.time()[3] - start.time[3], 1), "seconds"))
  return(structure(outlist, class = "fmlogit"))
}

