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
fitFactorCopula <- function(Y, copFun, k = rep(1, ncol(Y)), S, lower, upper, seed = runif(1, 1, .Machine$integer.max),
                            method = c("DEoptim", "genoud", "subplex", "two-stage"), control, cl = NULL, trials = NULL) {

  method <- match.arg(method)

  Yres <- apply(Y, 2, empDist)

  T <- nrow(Yres)
  N <- ncol(Yres)
  mHat <- moments(Yres, k)
  W <- diag(length(mHat))
  optim <- function(theta, S){
    names(theta) <- names(lower)
    U <- copFun(theta, S, seed)
    mTilde <- moments(U, k)
    res <- as.vector(Q(mTilde, mHat, W))
    res
  }

  if (method == "two-stage"){
    if (is.null(trials)){
      stop("Please specify the number of trials for the subplex method.")
    }
    # first approximate the global optimum
    checkNamespace("DEoptim")
    controlDE = list(c = 0.4, itermax = 100, reltol = 1e-6, steptol = 20, cl = cl, trace = FALSE)
    res <- DEoptim::DEoptim(optim, lower = lower, upper = upper, control = controlDE, S = 0.2*S)
    theta0 <- res$optim$bestmem
    # the use local optimizer
    res <- fitFactorCopulaSubplex(trials, lower, upper, optim, theta0 = theta0, cl = cl, control = control, S = S)
  }

  if (method == "DEoptim"){
    checkNamespace("DEoptim")
    control$cl <- cl
    res <- DEoptim::DEoptim(optim, lower = lower, upper = upper, control = control, S = S)
  }
  if (method == "genoud"){
    checkNamespace("rgenoud")
    if (is.null(cl)){
      cluster <- FALSE
    } else {
      cluster <- cl
    }
    control <- c(control, list(fn = optim, nvars = length(lower), Domains = matrix(c(lower, upper), ncol = 2),
                               boundary.enforcement = 2, P9 = 0, BFGS = FALSE, hessian = FALSE, cluster = cluster,
                               optim.method = "Nelder-Mead", S = S))
    res <- do.call(rgenoud::genoud, control)
  }

  if (method == "subplex"){
    if (is.null(trials)){
      stop("Please specify the number of trials for the subplex method.")
    }
    res <- fitFactorCopulaSubplex(trials, lower, upper, optim, cl = cl, control = control, S = S)
  }
  return(res)
}

fitFactorCopulaSubplex <- function(trials, lower, upper, optim, S, theta0 = NULL, control, cl){
  res <- parallelLapply(1:trials, function(i){
    if (is.null(theta0)){
      theta0 <- runif(length(lower), lower, upper)
    } else {
      theta0 <- runif(length(theta0), theta0 - abs(theta0*0.1), theta0 + abs(theta0*0.1))
      theta0[theta0 < lower] <- lower[theta0 < lower]
      theta0[theta0 > upper] <- upper[theta0 > upper]
    }
    nloptr::sbplx(theta0, optim, lower = lower, upper = upper, control = control, S = S)
  }, cl = cl)

  theta <- do.call(rbind, lapply(res, function(x) x$par))
  colnames(theta) <- names(lower)
  Qval <- unlist(lapply(res, function(x) x$value))
  converged <- unlist(lapply(res, function(x) x$convergence != 5))
  res <- cbind(theta, Q = Qval, convergence = converged)
  return(res)
}


parallelLapply <- function(x, fun, cl, ...){
  if(is.null(cl)){
    lapply(x, fun, ...)
  } else {
    parLapply(cl, x, fun, ...)
  }
}







