#' Title
#'
#' @param theta
#' @param tSeq
#'
#' @return
#' @export
#'
#' @examples
getPStat <- function(theta, tSeq){
  stopifnot(!is.null(dim(theta)))

  T <- max(tSeq)
  diff <- apply(theta, 1, function(x) sum(x - theta[nrow(theta), ])^2)
  P <- tSeq^2/T*diff
  estBrkPnt <- tSeq[which.max(P)]
  return(list(P = P, breakPoint = estBrkPnt))
}

unitVector <- function(size, k){
  vec <- rep(0, size)
  vec[k] <- 1
  return(vec)
}

# Structural Break Test ---------------------------------------------------

critVal <- function(Y, k = rep(1, ncol(Y)), B, copFun, theta, eps, cl){
  #resY: matrix of standardized residuals from empirical data Y
  #B: number of bootstrap samples
  #moments: function which generates a vector of dependence measures
  #eps: lower bound for the fraction of observations s to use
  mHead <- moments(Y, k)
  T <- nrow(Y)
  tSeq <- seq(floor(eps*T), T, 1)
  W <- diag(length(mHead))

  derivG <- G(theta, copFun, 0.1, S = 25000, mHead, k)
  AstarLeft <- solve(t(derivG) %*%W %*% derivG) %*% t(derivG) %*% W

  K <- parLapply(1:B, function(b){
    sampleInd <- sample(1:T, T, replace = TRUE)
    Astar <- lapply(tSeq, function(t){
      mHeadB <- moments(Y[sampleInd[1:t], ], k)
      A <- (mHeadB - mHead)*sqrt(T)*t/T
      AstarLeft %*% A
    })
    Kb <- lapply(seq_along(Astar), function(i) {
      diff <- Astar[[i]] - tSeq[i]/T*Astar[[length(Astar)]]
      t(diff) %*% diff
    })
    max(unlist(Kb))
  }, cl = cl)
  unlist(K)
}

G <- function(theta, copFun, eps, mHead, k, S){
  Gcol <- lapply(seq_along(theta), function(i){
    upper <- theta + unitVector(length(theta), i)*eps
    lower <- theta - unitVector(length(theta), i)*eps
    gUpper <- mHead - moments(copFun(upper, S), k)
    gLower <- mHead - moments(copFun(lower, S), k)
    (gUpper - gLower)/2*eps
  })
  return(do.call(cbind, Gcol))
}
