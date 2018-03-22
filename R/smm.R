#' Fit a factor copula model
#'
#' @param Y A dataframe or matrix like object
#' @param copFun A copula function estimated by @references factorCopula
#' @param k A vector with length ncol(Y) defining the groups (e.q. equi-dependence or block-euqidependence model)
#' @param S The number of simulations to use
#' @param lower Lower bound for optimazation: Named vector with parameters
#' @param upper Upper bound for optimazation: Named vector with parameters
#' @param recursive Wether to estimate recursive or just full model
#' @param control named list of arguments passed to the Subplex Algorithm (see nl.opts for help)
#'
#' @return a matrix or vector of parameters
#' @export
fc_fit <- function(Y, config_factor, config_error, config_beta, lower, upper, recursive, control, S, k, cl = NULL, trials = length(cl), load.balancing = TRUE) {

  Yres <- apply(Y, 2, empDist)
  mHat <- moments(Yres, k)

  T <- nrow(Yres)
  if (T <= 300){
    stop("Number of observations is to small.")
  }

  tSeq <- 300:T
  N <- ncol(Yres)
  W <- diag(length(mHat))

  opti <- function(theta, mHat, copFun, seed){
    names(theta) <- names(lower)
    gVal <-optim_g(mHat, copFun, theta, S, seed, k)
    as.vector(optim_q(gVal, W))
  }

  snow::clusterExport(cl, ls(envir = environment()), environment())

  if (recursive){
    ## estimate some starting values by slicing the input into chunks of 500 observations
    t_start <- slice(1:nrow(Y), 500)
    cat("Recursive model estimation with", length(t_start), "starting value(s) and", length(tSeq), "time periods\n")

    theta_start <- lapply(t_start, function(tValues){
      snow::clusterExport(cl, "t", environment())
      start <- parallelLapply(x = 1:length(cl), fun = function(trial){
        mHat <- moments(Yres[tValues, ], k)
        theta <- runif(length(lower), lower, upper)
        model_estimate(theta = theta, mHat = mHat, fc_create(config_factor, config_error, config_beta))
      }, cl = cl, load.balancing = load.balancing)
      start <- model_best(start)
      cat("Estimated starting value(s) from t =", min(tValues), "to t =", max(tValues), ":", round(start$par,4), "- Q:",round(start$value,4), "\n")
      start$par
    })

    snow::clusterExport(cl, c("theta_start"), environment())

    result <- parallelLapply(x = tSeq, fun = function(t){
      models <- lapply(theta_start, model_estimate,
                       mHat = moments(Yres[1:t, ], k),
                       copFun = fc_create(config_factor, config_error, config_beta))
      model_best(models)
    }, cl = cl, load.balancing = load.balancing)

    theta <- model_theta(result)
    theta$t <- tSeq

  } else {

    cat("Full model estimation with",trials, "trials\n")
    full <- parallelLapply(x = 1:trials, fun = function(trial){
      model_estimate(runif(length(lower), lower, upper), mHat, fc_create(config_factor, config_error, config_beta))
    }, cl = cl, load.balancing = load.balancing)
    best <- model_best(full)
    theta <- c(best$par, best$value, best$convergence, T)
  }

  names(theta) <- c(names(lower), "Q", "convergence", "t")
  return(theta)
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
  runif(1, -max, max)
}

slice <- function(x, n) split(x, as.integer((seq_along(x) - 1) / n))

model_estimate <- function(theta, mHat, copFun){
  seed <- random_seed()
  nloptr::sbplx(x0 = theta, fn = opti, lower = lower, upper = upper,
                control = control, mHat = mHat, copFun = copFun, seed = seed)
}

optim_q <- function(gVal, W){
  t(gVal)%*%W%*%gVal
}

optim_g <- function(mHat, copFun, theta, S, seed, k){
  U <- copFun(theta, S, seed)
  mTilde <- moments(U, k)
  return(mHat - mTilde)
}


