getQVal <- function(gVal, W){
  t(gVal)%*%W%*%gVal
}

getGVal <- function(mHat, copFun, theta, S, seed, k){
  U <- copFun(theta, S, seed)
  mTilde <- moments(U, k)
  return(mHat - mTilde)
}

test <- function(a){
  f <- function(){
    a <<- a+1
    return(a)
  }
}

getGHat <- function(theta, copFun, eps, mHat, k, S){
  P <- length(theta)
  M <- length(mHat)

  Gcol <- lapply(1:P, function(j){
    step <- unitVector(P, j)*eps
    ghatplus <- getGVal(mHat, copFun, theta + step, S, TRUE, k)
    ghatminus <- getGVal(mHat, copFun, theta - step, S, TRUE, k)
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
#' @param recursive Wether to estimate recursive or just full model
#' @param control named list of arguments passed to the Subplex Algorithm (see nl.opts for help)
#'
#' @return a matrix or vector of parameters
#' @export
fc_fit <- function(Y, config_factor, config_error, config_beta, lower, upper, recursive, control, S, k, cl) {

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
    gVal <-getGVal(mHat, copFun, theta, S, seed, k)
    as.vector(getQVal(gVal, W))
  }

  snow::clusterExport(cl, ls(envir = environment()), environment())

  if (recursive){
    ## estimate some starting values by slicing the input into chunks of 500 observations

    t_start <- slice(1:T, 500)

    theta_start <- lapply(t_start, function(t){
      start <- snow::clusterApplyLB(cl, 1:length(cl), function(trial){
        seed <- random_seed()
        copFun <- fc_create(config_factor, config_error, config_beta)
        theta0 <- runif(length(lower), lower, upper)
        names(theta0) <- names(lower)
        nloptr::sbplx(x0 = theta0, fn = opti, lower = lower, upper = upper,
                      control = control, mHat = moments(Yres[t, ], k), copFun = copFun, seed = seed)
      })
      start <- model_best(start)
      cat("Estimated starting value(s) from t =", min(t), "to t =", max(t), ":", round(start$par,4), "- Q:",round(start$value,4), "\n")
      start$par
    })

    snow::clusterExport(cl, c("theta_start"), environment())

    result <- snow::clusterApplyLB(cl, tSeq, function(t){
      models <- lapply(theta_start, function(theta){
        seed <- random_seed()
        copFun <- fc_create(config_factor, config_error, config_beta)
        names(theta) <- names(lower)
        nloptr::sbplx(x0 = theta, fn = opti, lower = lower, upper = upper,
                      control = control, mHat = moments(Yres[1:t, ], k), copFun = copFun, seed = seed)
      })
      model_best(models)
    })

    theta <- model_theta(result)
    theta$t <- tSeq

  } else {
    cat("Full model estimation with",length(cl), "trials\n")

    full <- snow::clusterApplyLB(cl, 1:length(cl), function(trial){
      seed <- random_seed()
      copFun <- fc_create(config_factor, config_error, config_beta)
      theta0 <- runif(length(lower), lower, upper)
      names(theta0) <- names(lower)
      nloptr::sbplx(x0 = theta0, fn = opti, lower = lower, upper = upper,
                    control = control, mHat = mHat, copFun = copFun, seed = seed)
    })
    best <- model_best(full)
    cat("Best model value(s):", round(best$par, 4), "- Q:", round(best$value, 4), "\n")
    theta <- round(c(best$par, best$value, best$convergence, T), 4)
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

# model_estimate <- function(cl, trials, theta0, from, to){
#   snow::clusterApplyLB(cl, 1:trials, function(trial){
#     seed <- random_seed()
#     copFun <- fc_create(config_factor, config_error, config_beta)
#     nloptr::sbplx(x0 = theta0[[trial]], fn = opti, lower = lower, upper = upper,
#                   control = control, mHat = moments(Yres[from:to, ], k), copFun = copFun, seed = seed)
#   })
# }


