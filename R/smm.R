Q <- function(x, y, W){
  diff <- x - y
  t(diff)%*%W%*%diff
}

#' Fit a factor copula model
#'
#' @param Y A dataframe or matrix like object
#' @param copFun A copula function estimated by @references factorCopula
#' @param k A vector with length ncol(Y) defining the groups (e.q. equi-dependence or block-euqidependence model)
#' @param S The number of simulations to use
#' @param lower Lower bound for optimazation: Named vector with parameters
#' @param upper Upper bound for optimazation: Named vector with parameters
#' @param seed Fixed random number seed to use during optimazation
#' @param method "DEoptim", "genoud" or "subplex" (see Details)
#' @param control named list of arguments passed to the optimizer
#'
#' @return a list of optimazation results
#' @export
fitFactorCopula <- function(Y, copFun, k = rep(1, ncol(Y)), S, lower, upper, seed = runif(1, 1, .Machine$integer.max), method = c("DEoptim", "genoud", "subplex"), control) {
  method <- match.arg(method)

  Yres <- apply(Y, 2, empDist)

  T <- nrow(Yres)
  N <- ncol(Yres)
  mHat <- moments(Yres, k)
  W <- diag(length(mHat))
  optim <- function(theta){
    names(theta) <- names(lower)
    U <- copFun(theta, S, seed)
    mTilde <- moments(U, k)
    res <- as.vector(Q(mTilde, mHat, W))
    res
  }

  if (method == "DEoptim"){
    checkNamespace("DEoptim")
    res <- DEoptim::DEoptim(optim, lower = lower, upper = upper, control = control)
  }
  if (method == "genoud"){
    checkNamespace("rgenoud")
    control <- c(control, list(fn = optim, nvars = length(lower), Domains = matrix(c(lower, upper), ncol = 2),
                               boundary.enforcement = 2, P9 = 0, BFGS = FALSE, hessian = FALSE,
                               optim.method = "Nelder-Mead"))
    res <- do.call(rgenoud::genoud, control)
  }
  if (method == "subplex"){
    if(is.null(control$runs)){
      runs <- 16
    } else {
      runs <- control$runs
    }
    control$runs <- NULL
    if (is.null(control$cl)){
      res <- lapply(1:runs, function(i) {
        theta0 <- sapply(seq_along(lower), function(j) runif(1, lower[j], upper[j]))
        nloptr::sbplx(theta0, optim, lower = lower, upper = upper, control = control)
      })
    } else {
      cl <- control$cl
      control$cl <- NULL
      res <- parLapply(cl, 1:runs, function(i) {
        theta0 <- sapply(seq_along(lower), function(j) runif(1, lower[j], upper[j]))
        nloptr::sbplx(theta0, optim, lower = lower, upper = upper, control = control)
      })
    }
    min <- which.min(vapply(res, function(x) x$value, numeric(1)))
    res <- res[[min]]
    names(res$par) <- names(lower)
  }
  return(res)
}








