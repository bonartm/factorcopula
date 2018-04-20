#' Fit a factor copula model
#'
#' @param Y A dataframe or matrix like object
#' @param factor specification of latent variables, see \link[factorcopula]{config_factor}
#' @param error specification of error term, see \link[factorcopula]{config_error}
#' @param beta specification of parameter matrix, see \link[factorcopula]{config_beta}
#' @param lower Lower bound for optimazation: Named vector with parameters
#' @param upper Upper bound for optimazation: Named vector with parameters
#' @param control named list of arguments passed to the subplex algorithm, see \link[nloptr]{nl.opts}
#' @param S The number of simulations to use
#' @param k A vector with length ncol(Y) defining the groups (e.q. equi-dependence or block-euqidependence model)
#' @param se Wether to estimate standard errors and confidence intervalls
#' @param cl A cluster object, see \link[snow]{makeCluster}
#' @param trials number of model runs with different starting values
#' @param load.balancing if TRUE a load balancing cluster apply is performed
#' @return if recursive a data.frame, else a vector of parameters and model statistics.
#'
#' @section Convergence info:
#'
#' \emph{taken from: \url{http://nlopt.readthedocs.io/en/latest/NLopt_Reference/#Return_values.md}.}
#'
#' Postive integers:
#' \itemize{
#' \item 1 - Generic success return value.
#' \item 2 - Optimization stopped because stopval (above) was reached.
#' \item 3 - Optimization stopped because ftol_rel or ftol_abs (above) was reached.
#' \item 4 - Optimization stopped because xtol_rel or xtol_abs (above) was reached.
#' \item 5 - Optimization stopped because maxeval (above) was reached.
#' \item 6 - Optimization stopped because maxtime (above) was reached.
#' }
#' Negative integers:
#' \itemize{
#' \item -1 - Generic failure code.
#' \item -2 - Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).
#' \item -3 - Ran out of memory.
#' \item -4 - Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)
#' \item-5 - Halted because of a forced termination
#' }
#'
#' @export
fc_fit <- function(Y, factor, error, beta, lower, upper, control, S, k, se = FALSE,
                   cl = NULL, trials = max(1, length(cl)), load.balancing = TRUE) {

  opti <- function(theta, mHat, copFun, seed){
    names(theta) <- names(lower)
    gVal <-optim_g(mHat, copFun, theta, S, seed, k)
    as.vector(optim_q(gVal, W))
  }

  model_estimate <- function(theta, mHat){
    seed <- random_seed()
    copFun <- fc_create(factor, error, beta)
    nloptr::sbplx(x0 = theta, fn = opti, lower = lower, upper = upper,
                  control = control, mHat = mHat, copFun = copFun, seed = seed)
  }

  Yres <- apply(Y, 2, empDist)
  mHat <- moments(Yres, k)
  T <- nrow(Yres)
  N <- ncol(Yres)
  W <- diag(length(mHat))

  if(!is.null(cl)){
    snow::clusterExport(cl, ls(envir = environment()), environment())
  }

  cat("Full model estimation with",trials, "trial(s).\n")
  models <- parallelLapply(x = 1:trials, fun = function(trial){
      model_estimate(stats::runif(length(lower), lower, upper), mHat)
  }, cl = cl, load.balancing = load.balancing)
  best <- model_best(models, names(lower))
  all <- model_theta(models, names(lower))

  res <- list(models = all, best = best)

  if (se){
    cat("Estimating standard errors using the best model.\n")
    sigma <- getSigmaHat(Yres, 1000, k)
    copFun <- fc_create(factor, error, beta)
    seed <- random_seed()
    thetaFull <- best[1:length(lower)]

    G <- getGHat(thetaFull, copFun, 0.1, mHat, k, S, seed)

    omega <- getOmegaHat(G, W, sigma)
    se <- sqrt(omega*(1/T + 1/S))
    lowerCI <- thetaFull - qnorm(1-0.05/2) * se
    upperCI <- thetaFull + qnorm(1-0.05/2) * se
    res$omega = omega
    res$se = se
    res$ci = data.frame(par = thetaFull, lower = lowerCI, upper = upperCI)
  }

  return(res)
}


model_theta <- function(models, namesTheta){
  res <- data.frame(do.call(rbind, lapply(models, function(x) c(x$par, x$value, x$convergence))))
  names(res) <- c(namesTheta, "Q", "convergence")
  res
}

model_best <- function(models, namesTheta){
  Qval <- vapply(models, function(x) x$value, numeric(1))
  best <- models[which.min(Qval)]
  unlist(model_theta(best, namesTheta))
}

random_seed <- function(){
  max <- .Machine$integer.max
  stats::runif(1, -max, max)
}

optim_q <- function(gVal, W){
  t(gVal)%*%W%*%gVal
}

optim_g <- function(mHat, copFun, theta, S, seed, k){
  U <- copFun(theta, S, seed)
  mTilde <- moments(U, k)
  return(mHat - mTilde)
}

getSigmaHat <- function(Ydis, B, k){
  T <- nrow(Ydis)
  sigmaB <- replicate(B, {
    b <- sample(1:T, T, replace = TRUE)
    mHatB <- moments(Ydis[b, ], k)
  })
  T*stats::cov(t(sigmaB))
}

getOmegaHat <- function(G, W, sigma){
  inv <- solve(t(G)%*%W%*%G)
  inv%*%t(G)%*%W%*%sigma%*%W%*%G%*%inv
}


