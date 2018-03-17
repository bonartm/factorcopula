#' Return the p statistics of a recursive factor copula model
#'
#' @param theta A matrix with max(tSeq) rows of recursive parameters
#' @param tSeq A vector of positive integers
#'
#' @return for each entray in tSeq a P statistic
#' @export
fc_pstat <- function(theta, tSeq){
  stopifnot(!is.null(dim(theta)))
  T <- max(tSeq)
  diff <- apply(theta, 1, function(x) {
    diff <- (x - theta[nrow(theta), ])
    t(diff)%*%diff
  })
  P <- (tSeq/T)^2*T*diff
  return(P)
}

#' Calculate recursive test statistics for the moments based break test
#'
#' @export
fc_mstat <- function(Y, tSeq, k = rep(1, ncol(Y))){
  Ydis <- apply(Y, 2, empDist)
  mFull <- moments(Ydis, k)
  T <- nrow(Y)
  mStats <- vapply(tSeq, function(t){
    diff <- moments(Ydis[1:t, ], k) - mFull
    (t/T)^2*T*(t(diff) %*% diff)
  }, numeric(1))
  return(mStats)
}


#' Simulate critival values for either the copula or the moments based break test
#'
#' @export
fc_critval <- function(type = c("moments", "copula"), Y, B, tSeq, copFun = NULL, theta = NULL, cl = NULL, k = rep(1, ncol(Y))){
  #resY: matrix of standardized residuals from empirical data Y
  #B: number of bootstrap samples
  #moments: function which generates a vector of dependence measures
  type <- match.arg(type)
  Ydis <- apply(Y, 2, empDist)
  mHat <- moments(Ydis, k)
  T <- nrow(Y)

  if (type == "copula"){
    seed <- runif(1, min = 1, max = .Machine$integer.max)
    G <- getGHat(thetaFull, copFun, 0.1, mHat, k, S = 25000, seed)
    W <- diag(length(mHat))
    leftPart <- solve(t(G)%*%W%*%G)%*%t(G)%*%W
  } else {
    leftPart <- 1
  }

  Kb <- parallelLapply(1:B, function(x) {
    b <- sample(1:T, T, replace = TRUE)

    Astar <- lapply(tSeq, function(t) {
      mHatBt <- moments(Ydis[b, ][1:t, ], k)
      A <- t/T*sqrt(T)*(mHatBt - mHat)
      matrix(as.vector(leftPart %*% A), nrow = length(mHat))
    })

    Kbt <- vapply(seq_along(tSeq), function(i){
      t <- tSeq[i]
      diff <- Astar[[i]] - t/T*Astar[[length(tSeq)]]
      t(diff)%*%diff
    }, numeric(1))

    return(max(Kbt))
  }, cl = cl)
  return(unlist(Kb))
}




