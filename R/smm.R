getQVal <- function(gVal, W){
  t(gVal)%*%W%*%gVal
}

getGVal <- function(mHat, copFun, theta, S, fixed, k){
  U <- copFun(theta, S, fixed)
  mTilde <- moments(U, k)
  return(mHat - mTilde)
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
fc_fit <- function(Y, copFun, lower, upper, recursive, control, S, k, cl) {

  Yres <- apply(Y, 2, empDist)
  mHat <- moments(Yres, k)

  T <- nrow(Yres)
  if (T <= 300){
    stop("Number of observations is to small.")
  }
  tSeq <- 300:T
  N <- ncol(Yres)
  W <- diag(length(mHat))

  opti <- function(theta, mHat){
    names(theta) <- names(lower)
    gVal <-getGVal(mHat, copFun, theta, S, TRUE, k)
    as.vector(getQVal(gVal, W))
  }

  if (recursive){
    cat("Recursive model estimation\n")
    snow::clusterExport(cl, ls(envir = environment()), environment())
    start <- snow::clusterApplyLB(cl, 1:length(cl), function(trial){
      theta0 <- runif(length(lower), lower, upper)
      names(theta0) <- names(lower)
      nloptr::sbplx(x0 = theta0, fn = opti, lower = lower, upper = upper,
                    control = control, mHat = moments(Yres[1:min(tSeq), ], k))
    })
    start <- model_best(start)
    names(start$theta) <- names(lower)

    cat("Estimated starting value(s):", start$theta, "- Q:",start$Q, "\n")
    snow::clusterExport(cl, list("start"), environment())

    result <- snow::clusterApplyLB(cl, tSeq, function(t){
      nloptr::sbplx(x0 = start$theta, fn = opti, lower = lower, upper = upper,
                    control = control, mHat = moments(Yres[1:t, ], k))
    })
    theta <- data.frame(t = tSeq, model_theta(result))
  } else {
    cat("Full model estimation with",length(cl), "trials\n")
    snow::clusterExport(cl, ls(envir = environment()), environment())

    full <- snow::clusterApplyLB(cl, 1:length(cl), function(trial){
      theta0 <- runif(length(lower), lower, upper)
      names(theta0) <- names(lower)
      nloptr::sbplx(x0 = theta0, fn = opti, lower = lower, upper = upper,
                    control = control, mHat = mHat)
    })
    theta <- model_theta(full)
    for (i in 1:nrow(theta)){
      cat(i, "Estimated value(s):", theta[i, ], "\n")
    }
    best <- model_best(full)
    cat("Best model value(s):", best$theta, "- Q:",best$Q, "\n")
    theta <- data.frame(t = T, best$theta)
  }
  names(theta) <- c("t", names(lower))
  return(theta)
}

model_theta <- function(models){
  do.call(rbind, lapply(models, function(x) round(x$par,4)))
}

model_best <- function(models){
  theta <- model_theta(models)
  Qval <- vapply(models, function(x) round(x$value,4), numeric(1))
  list(theta = theta[which.min(Qval), , drop = FALSE], Q = Qval[which.min(Qval)])
}
