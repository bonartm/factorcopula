#' Return the p statistics of a recursive factor copula model
#'
#' @param theta A numeric matrix with length(tSeq) rows of recursive parameters
#' @param tSeq A vector of positive integers
#'
#' @return a vector of P statistics for each t
#' @export
fc_pstat <- function(theta, tSeq){
  stopifnot(!is.null(dim(theta)))
  theta <- as.matrix(theta)
  T <- max(tSeq)
  thetaFull <- theta[nrow(theta), ]

  leftPart <- vapply(1:nrow(theta), function(i){
    diff <- theta[i, ] - thetaFull
    t(diff)%*%diff
  }, numeric(1))

  P <- (tSeq/T)^2*T*leftPart
  return(P)
}

#' Calculate recursive test statistics for the moments based break test
#' @param Y matrix of observed values
#' @param tSeq A vector of positive integers
#' @param k a vector defining the groups of the variables
#' @param cl a cluster object, see \link[snow]{makeCluster}
#'
#' @export
fc_mstat <- function(Y, tSeq, k, cl = NULL){
  Ydis <- apply(Y, 2, empDist)
  mFull <- moments(Ydis, k)
  T <- nrow(Y)
  mStats <- parallelLapply(tSeq, function(t, mFull, T, Ydis){
    diff <- moments(Ydis[1:t, ], k) - mFull
    (t/T)^2*T*(t(diff) %*% diff)
  }, cl, mFull = mFull, T = T, Ydis = Ydis)
  return(unlist(mStats))
}


#' Simulate critival values for either the copula or the moments based break test
#' @param type either moments or copula
#' @param Y a numeric matrix of observed values
#' @param B the number of bootstrap replications
#' @param tSeq a sequence with positive integers
#' @param k a vector defining the groups of the variables
#' @param factor specification of latent variables, see \link[factorcopula]{config_factor}
#' @param error specification of error term, see \link[factorcopula]{config_error}
#' @param beta specification of parameter matrix, see \link[factorcopula]{config_beta}
#' @param theta named vector of full model parameter estimates
#' @param cl an optional cluster object, see \link[snow]{makeCluster}
#'
#' @export
fc_critval <- function(type = c("moments", "copula"), Y, B, tSeq, k,
                       factor = NULL, error = NULL, beta = NULL, theta = NULL, cl = NULL){
  #resY: matrix of standardized residuals from empirical data Y
  #B: number of bootstrap samples
  #moments: function which generates a vector of dependence measures
  type <- match.arg(type)
  Ydis <- apply(Y, 2, empDist)
  mHat <- moments(Ydis, k)
  T <- nrow(Y)

  if (type %in% c("both", "copula")){
    stopifnot(!is.null(config_factor) & !is.null(config_error) & !is.null(config_beta))
    copFun <- fc_create(config_factor, config_error, config_beta)
    theta <- unlist(c(theta))
    seed <- random_seed()
    G <- getGHat(theta, copFun, 0.1, mHat, k, S = 25000, seed)
    W <- diag(length(mHat))
    leftPart <- solve(t(G)%*%W%*%G)%*%t(G)%*%W
  } else {
    leftPart <- NULL
  }

  Kb <- parallelLapply(1:B, function(x, T, tSeq, Ydis, k, leftPart, type) {
    b <- sample(1:T, T, replace = TRUE)

    A <- lapply(tSeq, function(t) {
      mHatBt <- moments(Ydis[b, ][1:t, ], k)
      A <- t/T*sqrt(T)*(mHatBt - mHat)
      if (type == "copula"){
        A <- leftPart %*% A
      }
      return(A)
    })

    Kbt <- vapply(seq_along(tSeq), function(i){
      t <- tSeq[i]
      diff <- A[[i]] - t/T*A[[length(tSeq)]]
      t(diff)%*%diff
    }, numeric(1))

    return(max(Kbt))
  }, cl = cl, T = T, tSeq = tSeq, Ydis = Ydis, k = k, leftPart = leftPart, type = type, load.balancing = FALSE)
  return(unlist(Kb))
}

getGHat <- function(theta, copFun, eps, mHat, k, S, seed){
  P <- length(theta)
  M <- length(mHat)

  Gcol <- lapply(1:P, function(j){
    step <- unitVector(P, j)*eps
    ghatplus <- optim_g(mHat, copFun, theta + step, S, seed, k)
    ghatminus <- optim_g(mHat, copFun, theta - step, S, seed, k)
    (ghatplus-ghatminus)/(2*eps)
  })
  G <- do.call(cbind, Gcol)
  stopifnot(nrow(G) == M & ncol(G) == P)
  return(G)
}

unitVector <- function(size, k){
  vec <- rep(0, size)
  vec[k] <- 1
  return(vec)
}




