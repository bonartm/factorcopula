getQVal <- function(gVal, W){
  t(gVal)%*%W%*%gVal
}

getGVal <- function(mHat, copFun, theta, S, seed, k){
  U <- copFun(theta, S, seed)
  mTilde <- moments(U, k)
  return(mHat - mTilde)
}

getGHat <- function(theta, copFun, eps, mHat, k, S, seed){
  P <- length(theta)
  M <- length(mHat)

  Gcol <- lapply(1:P, function(j){
    step <- unitVector(P, j)*eps
    ghatplus <- getGVal(mHat, copFun, theta + step, S, seed, k)
    ghatminus <- getGVal(mHat, copFun, theta - step, S, seed, k)
    (ghatplus-ghatminus)/(2*eps)
  })
  G <- do.call(cbind, Gcol)
  stopifnot(nrow(G) == M & ncol(G) == P)
  return(G)
}

getSigmaHat <- function(Y, B, k){
  T <- nrow(Y)
  Ydis <- apply(Y, 2, empDist)
  sigmaB <- replicate(B, {
    b <- sample(1:T, T, replace = TRUE)
    mHatB <- moments(Ydis[b, ], k)
  })
  T*cov(t(sigmaB))
}

getOmegaHat <- function(G, W, sigma){
  inv <- solve(t(G)%*%W%*%G)
  inv%*%t(G)%*%W%*%sigma%*%W%*%G%*%inv
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
fc_fit <- function(Y, copFun, lower, upper, method = c("DEoptim", "genoud", "subplex", "two-stage"), control,
                            cluster = NULL, trials = NULL, S = NULL, k = rep(1, ncol(Y)), seed = runif(1, 1, .Machine$integer.max)) {

  method <- match.arg(method)

  Yres <- apply(Y, 2, empDist)

  T <- nrow(Yres)
  N <- ncol(Yres)
  mHat <- moments(Yres, k)
  W <- diag(length(mHat))


  opti <- function(theta){
    names(theta) <- names(lower)
    gVal <-getGVal(mHat, copFun, theta, S, seed, k)
    as.vector(getQVal(gVal, W))
  }

  if (method == "two-stage"){
    if (is.null(trials)){
      stop("Please specify the number of trials for the subplex method.")
    }
    # first approximate the global optimum
    checkNamespace("DEoptim")
    controlDE = list(c = 0.4, itermax = 100, reltol = 1e-6, steptol = 20, cl = cluster, trace = FALSE)
    res <- DEoptim::DEoptim(opti, lower = lower, upper = upper, control = controlDE)
    theta0 <- res$optim$bestmem
    # the use local optimizer
    res <- fitFactorCopulaSubplex(trials, lower, upper, opti, theta0 = theta0, cl = cluster, control = control)
  }

  if (method == "DEoptim"){
    checkNamespace("DEoptim")
    control$cl <- cluster
    res <- DEoptim::DEoptim(opti, lower = lower, upper = upper, control = control)
  }
  if (method == "genoud"){
    checkNamespace("rgenoud")
    if (is.null(cluster))
      cluster <- FALSE

    control <- c(control, list(fn = opti, nvars = length(lower), Domains = matrix(c(lower, upper), ncol = 2),
                               boundary.enforcement = 2, P9 = 0, BFGS = FALSE, hessian = FALSE, cluster = cluster,
                               optim.method = "Nelder-Mead"))
    res <- do.call(rgenoud::genoud, control)
  }

  if (method == "subplex"){
    if (is.null(trials)){
      stop("Please specify the number of trials for the subplex method.")
    }

    res <- fitFactorCopulaSubplex(trials = trials, lower = lower, upper = upper, fn = opti, cluster = cluster, control = control)

  }
  return(res)
}



fitFactorCopulaSubplex <- function(trials, lower, upper, fn, control, cluster, theta0 = NULL){
  res <- parallelLapply(1:trials, function(i){
    if (is.null(theta0)){
      theta0 <- runif(length(lower), lower, upper)
    } else {
      theta0 <- runif(length(theta0), theta0 - abs(theta0*0.1), theta0 + abs(theta0*0.1))
      theta0[theta0 < lower] <- lower[theta0 < lower]
      theta0[theta0 > upper] <- upper[theta0 > upper]
    }
    nloptr::sbplx(x0 = theta0, fn = fn, lower = lower, upper = upper, control = control)
  }, cl = cluster)

  theta <- do.call(rbind, lapply(res, function(x) x$par))
  colnames(theta) <- names(lower)
  Qval <- unlist(lapply(res, function(x) x$value))
  iter <- unlist(lapply(res, function(x) x$iter))
  res <- list(Q = Qval, iterations = iter, params = theta, best = theta[which.min(Qval),])
  return(res)
}










