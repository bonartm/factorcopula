#' Return the p statistics of a recursive factor copula model
#'
#' @param theta A numeric matrix with max(tSeq) rows of recursive parameters
#' @param tSeq A vector of positive integers
#'
#' @return for each entray in tSeq a P statistic
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
#'
#' @export
fc_mstat <- function(Y, tSeq, k = rep(1, ncol(Y)), cl = NULL){
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
#'
#' @export
fc_critval <- function(type = c("moments", "copula"), Y, B, tSeq, k, copFun = NULL, theta = NULL, cl = NULL){
  #resY: matrix of standardized residuals from empirical data Y
  #B: number of bootstrap samples
  #moments: function which generates a vector of dependence measures
  type <- match.arg(type)
  Ydis <- apply(Y, 2, empDist)
  mHat <- moments(Ydis, k)
  T <- nrow(Y)

  if (type %in% c("both", "copula")){
    seed <- runif(1, min = 1, max = .Machine$integer.max)
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
  }, cl = cl, T = T, tSeq = tSeq, Ydis = Ydis, k = k, leftPart = leftPart, type = type)
  return(unlist(Kb))
}

getGHat <- function(theta, copFun, eps, mHat, k, S){
  P <- length(theta)
  M <- length(mHat)

  Gcol <- lapply(1:P, function(j){
    step <- unitVector(P, j)*eps
    ghatplus <- otim_g(mHat, copFun, theta + step, S, TRUE, k)
    ghatminus <- optim_g(mHat, copFun, theta - step, S, TRUE, k)
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




