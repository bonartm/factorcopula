#' Fit a factor copula model
#'
#' @param Y A dataframe or matrix like object
#' @param factor specification of latent variables, see \link[factorcopula]{config_factor}
#' @param error specification of error term, see \link[factorcopula]{config_error}
#' @param beta specification of parameter matrix, see \link[factorcopula]{config_beta}
#' @param lower Lower bound for optimazation: Named vector with parameters
#' @param upper Upper bound for optimazation: Named vector with parameters
#' @param recursive Wether to estimate recursive or full model
#' @param control named list of arguments passed to the subplex algorithm, see \link[nloptr]{nl.opts}
#' @param S The number of simulations to use
#' @param k A vector with length ncol(Y) defining the groups (e.q. equi-dependence or block-euqidependence model)
#' @param cl A cluster object, see \link[snow]{makeCluster}
#' @param trials number of model runs with different starting values
#' @param load.balancing if TRUE a load balancing cluster apply is performed
#' @return if recursive a data.frame, else a vector of parameters and model statistics
#' @export
fc_fit <- function(Y, factor, error, beta, lower, upper, recursive, control, S, k,
                   cl = NULL, trials = max(1, length(cl)), load.balancing = TRUE) {

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

  opti <- function(theta, mHat, copFun, seed){
    names(theta) <- names(lower)
    gVal <-optim_g(mHat, copFun, theta, S, seed, k)
    as.vector(optim_q(gVal, W))
  }

  if(!is.null(cl)){
    snow::clusterExport(cl, ls(envir = environment()), environment())
  }

  if (recursive){
    if (T <= 300){
      stop("Number of observations is to small.")
    }
    tSeq <- 300:T
    ## estimate some starting values by slicing the input into chunks of 500 observations
    t_start <- slice(1:nrow(Y), 500)
    cat("Recursive model estimation with", length(t_start), "starting value(s) and", length(tSeq), "time periods\n")

    theta_start <- lapply(t_start, function(tValues){
      if(!is.null(cl)){
        snow::clusterExport(cl, "t", environment())

      }
      start <- parallelLapply(x = 1:length(cl), fun = function(trial){
        mHat <- moments(Yres[tValues, ], k)
        theta <- stats::runif(length(lower), lower, upper)
        model_estimate(theta = theta, mHat = mHat)
      }, cl = cl, load.balancing = load.balancing)
      start <- model_best(start)
      cat("Estimated starting value(s) from t =", min(tValues), "to t =", max(tValues), ":", round(start$par,4), "- Q:",round(start$value,4), "\n")
      start$par
    })

    if(!is.null(cl)){
      snow::clusterExport(cl, c("theta_start"), environment())
    }

    result <- parallelLapply(x = tSeq, fun = function(t){
      models <- lapply(theta_start, model_estimate,
                       mHat = moments(Yres[1:t, ], k))
      model_best(models)
    }, cl = cl, load.balancing = load.balancing)

    theta <- model_theta(result)
    theta$t <- tSeq

  } else {

    cat("Full model estimation with",trials, "trial(s)\n")
    full <- parallelLapply(x = 1:trials, fun = function(trial){
      model_estimate(stats::runif(length(lower), lower, upper), mHat)
    }, cl = cl, load.balancing = load.balancing)
    best <- model_best(full)
    theta <- c(best$par, best$value, best$convergence, T)
  }

  names(theta) <- c(names(lower), "Q", "convergence", "t")
  return(round(theta, 4))
}


model_theta <- function(models){
  data.frame(do.call(rbind, lapply(models, function(x) c(x$par, x$value, x$convergence))))
}

model_best <- function(models){
  Qval <- vapply(models, function(x) x$value, numeric(1))
  models[[which.min(Qval)]]
}

random_seed <- function(){
  max <- .Machine$integer.max
  stats::runif(1, -max, max)
}

slice <- function(x, n) split(x, as.integer((seq_along(x) - 1) / n))

optim_q <- function(gVal, W){
  t(gVal)%*%W%*%gVal
}

optim_g <- function(mHat, copFun, theta, S, seed, k){
  U <- copFun(theta, S, seed)
  mTilde <- moments(U, k)
  return(mHat - mTilde)
}


