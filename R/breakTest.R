#' Title
#'
#' @param theta
#' @param tSeq
#'
#' @return
#' @export
getPStat <- function(theta, tSeq){
  stopifnot(!is.null(dim(theta)))
  T <- max(tSeq)
  diff <- apply(theta, 1, function(x) {
    diff <- (x - theta[nrow(theta), ])
    t(diff)%*%diff
  })
  P <- (tSeq/T)^2*T*diff
  return(P)
}



# Structural Break Test ---------------------------------------------------

critVal <- function(Y, B, copFun, theta, tSeq, cl = NULL, k = rep(1, ncol(Y))){
  #resY: matrix of standardized residuals from empirical data Y
  #B: number of bootstrap samples
  #moments: function which generates a vector of dependence measures
  Ydis <- apply(Y, 2, empDist)
  mHat <- moments(Ydis, k)
  T <- nrow(Y)
  seed <- runif(1, min = 1, max = .Machine$integer.max)

  G <- getGHat(theta, copFun, 0.1, mHat, k, S = 25000, seed)
  W <- diag(length(mHat))


  leftPart <- solve(t(G)%*%W%*%G)%*%t(G)%*%W

  Kb <- replicate(B, {
    b <- sample(1:T, T, replace = TRUE)

    Astar <- lapply(tSeq, function(t) {
      mHatBt <- moments(Ydis[b, ][1:t, ], k)
      A <- t/T*sqrt(T)*(mHatBt - mHat)
      leftPart %*% A
    })

    Kbt <- vapply(seq_along(tSeq), function(i){
      t <- tSeq[i]
      diff <- Astar[[i]] - t/T*Astar[[length(tSeq)]]
      t(diff)%*%diff
    }, numeric(1))

    return(max(Kbt))
  })
  return(Kb)
}




