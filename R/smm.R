#' Fit a factor copula model
#'
#' @param Y A dataframe or matrix like object
#' @param factor specification of latent variables, see \link[factorcopula]{config_factor}
#' @param error specification of error term, see \link[factorcopula]{config_error}
#' @param beta specification of parameter matrix, see \link[factorcopula]{config_beta}
#' @param lower Lower bound for optimazation: Named vector with parameters
#' @param upper Upper bound for optimazation: Named vector with parameters
#' @param control.first.stage named list of options passed to the nloptr method. A full description of all options is shown by the function \code{nloptr::nloptr.print.options()}.
#' @param control.second.stage named list of options passed to the nloptr method. A full description of all options is shown by the function \code{nloptr::nloptr.print.options()}.
#' @param S The number of simulations to use
#' @param k A vector with length ncol(Y) defining the groups (e.q. equi-dependence or block-euqidependence model)
#' @param se Wether to estimate standard errors and confidence intervalls
#' @param B number of bootstrap replication for estimating standard error
#' @return a list of optimization results
#'
#' @export
fc_fit <- function(Y, factor, error, beta, lower, upper, k = rep(1, ncol(Y)), S = 25*nrow(Y), se = FALSE, B = 1000,
                   control.first.stage = list(algorithm = "NLOPT_GN_DIRECT_L", stopval = 0, xtol_rel = 1e-4, maxeval = 200),
                   control.second.stage = list(algorithm = "NLOPT_LN_SBPLX", stopval = 0, xtol_rel = 1e-16, maxeval = 10000)) {

  stopifnot(!is.null(control.first.stage$algorithm))

  opti <- function(theta, mHat, copFun, seed){
    names(theta) <- names(lower)
    gVal <-optim_g(mHat, copFun, theta, S, seed, k)
    as.vector(optim_q(gVal, W))
  }

  model_estimate <- function(theta, mHat, control){
    seed <- random_seed()
    copFun <- fc_create(factor, error, beta)
    nloptr::nloptr(x0 = theta, eval_f = opti, lb = lower, ub = upper,
                  opts = control, mHat = mHat, copFun = copFun, seed = seed)
  }

  Yres <- apply(Y, 2, empDist)
  mHat <- moments(Yres, k)
  T <- nrow(Yres)
  N <- ncol(Yres)
  W <- diag(length(mHat))


  cat("First stage model estimation with", control.first.stage$algorithm, "algorithm.\n")
  model <- model_estimate(stats::runif(length(lower), lower, upper), mHat, control.first.stage)
  names(model$solution) <- names(lower)
  cat("First stage solution:", model$solution, "\n")

  ret <- list(theta.first.stage = model$solution)


  if (!is.null(control.second.stage$algorithm)){
    cat("Second stage model estimation with", control.second.stage$algorithm, "algorithm.\n")
    model <- model_estimate(model$solution, mHat, control.second.stage)
    names(model$solution) <- names(lower)
    cat("Second stage solution:", model$solution, "\n")
    ret$theta.second.stage <- model$solution
  }

  ret$Q <- model$objective
  ret$message <- model$message
  ret$iterations <- model$iterations


  if (se){
    cat("Standard error estimation.\n")
    sigma <- getSigmaHat(Yres, B, k)
    copFun <- fc_create(factor, error, beta)
    seed <- random_seed()
    thetaFull <- model$solution

    G <- getGHat(thetaFull, copFun, 0.1, mHat, k, S, seed)

    omega <- getOmegaHat(G, W, sigma)
    se <- sqrt(diag(omega)*(1/T + 1/S))
    lowerCI <- thetaFull - stats::qnorm(1-0.05/2) * se
    upperCI <- thetaFull + stats::qnorm(1-0.05/2) * se
    Tval <- (thetaFull-0)/se
    pVal <- 2*(stats::pnorm(-abs(Tval)))
    ret$se = se
    ret$ci = data.frame(par = thetaFull, lower = lowerCI, upper = upperCI, p = pVal)
  }

  return(ret)
}

model_best <- function(models){
  Qval <- vapply(models, function(x) x$objective, numeric(1))
  models[[which.min(Qval)]]
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


