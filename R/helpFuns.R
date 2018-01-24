checkNamespace <- function(name){
  if (!requireNamespace(name, quietly = TRUE)) {
    stop(paste(name, "needed for this function to work. Please install it.", call. = FALSE))
  }
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
getStandResiduals <- function(x){
  checkNamespace("rugarch")
  spec = rugarch::ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(1,0), include.mean=TRUE),
                             distribution.model="norm")
  garch <- rugarch::ugarchfit(spec = spec, data = x)
  res <- as.vector(rugarch::residuals(garch, standardize = TRUE))
  return(res)
}

empDist <- function(x){
  data.table::frank(x)/length(x)
}

rankCor <- function(u, v){
  # t <- length(u)
  # #12/t*sum(u*v) - 3
  # 12*(t/(t^2-1))*sum(u*v)-3*(t+1)/(t-1)
  cor(u, v)
}

rst <- function(n, nuInv = 1e-10, lambda = 0){
  # Generates random numbers from the skewed t distribution by Hansen (1994)
  # nuInv: Inverse nu parameter (degrees of freedom, q) (0, )
  # lambda: skewness parameter (-1, 1)
  nu <- 1/nuInv
  stopifnot(nu >= 2 & nuInv > 0 & lambda > -1 & lambda < 1)

  u <- runif(n)
  if (is.infinite(suppressWarnings(gamma(nu/2)))){
    c <- 0
    a <- 0
  } else {
    c <- gamma((nu+1)/2)/(sqrt(pi*(nu-2))*gamma(nu/2))
    a <- 4*lambda*c*((nu-2)/(nu-1));
  }

  b <- sqrt(1 + 3*lambda^2 - a^2);

  f1 <- u < (1-lambda)/2

  inv1 <- (1-lambda)/b*sqrt((nu-2)/nu)*qt(u[f1]/(1-lambda), nu)-a/b
  inv2 <- (1+lambda)/b*sqrt((nu-2)/nu)*qt(0.5+1/(1+lambda)*(u[!f1]-(1-lambda)/2), nu)-a/b

  inv <- rep(0, n)
  inv[f1]  <- inv1
  inv[!f1]  <- inv2
  return(inv)
}

rtInv <- function(n, dfInv){
  rt(n, 1/dfInv)
}

quantDep <- function(u, v, qSeq){
  t <- length(u)
  vapply(qSeq, function(q){
    if (q <= 0.5){
      sum(u <= q & v <= q)/(t*q)
    } else {
      sum(u > q & v > q)/(t*(1-q))
    }
  }, numeric(1))
}


#' Title
#'
#' @param k
#'
#' @return
#' @export
genBetaParMat <- function(k){
  kTab <- table(k)
  M <- max(k)
  first <- rep(paste0("beta", 1:M), times = kTab)
  last <- lapply(1:M, function(m){
    if (m == 1) timesLow <- 0 else timesLow <- sum(kTab[1:m-1])
    if (m == M) timesUp <- 0 else timesUp <- sum(kTab[(m+1):M])
    c(rep(0, timesLow), rep(paste0("beta", m+M), kTab[m]), rep(0, timesUp))
  })
  cbind(first, do.call(cbind, last))
}

