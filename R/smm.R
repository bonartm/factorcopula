Q <- function(x, y, W){
  diff <- x - y
  t(diff)%*%W%*%diff
}

fitFactorCopula <- function(Y, copFun, k = rep(1, ncol(Y)), S, lower, upper, seed = runif(1, 1, .Machine$integer.max), method = c("DEoptim", "genoud", "subplex"), control) {
  method <- match.arg(method)

  T <- nrow(Y)
  N <- ncol(Y)
  mHat <- moments(Y, k)
  W <- diag(length(mHat))
  optim <- function(theta){
    names(theta) <- names(lower)
    U <- copFun(theta, S, seed)
    mTilde <- moments(U, k)
    res <- as.vector(Q(mTilde, mHat, W))
    res
  }

  if (method == "DEoptim"){
    res <- DEoptim(optim, lower = lower, upper = upper, control = control)
  }
  if (method == "genoud"){
    control <- c(control, list(fn = optim, nvars = length(lower), Domains = matrix(c(lower, upper), ncol = 2),
                               boundary.enforcement = 2, P9 = 0, BFGS = FALSE, hessian = FALSE,
                               optim.method = "Nelder-Mead"))
    res <- do.call(genoud, control)
  }
  if (method == "subplex"){
    if(is.null(control$runs)){
      runs <- 16
    } else {
      runs <- control$runs
    }
    control$runs <- NULL
    if (is.null(control$cl)){
      res <- lapply(1:runs, function(i) {
        theta0 <- sapply(seq_along(lower), function(j) runif(1, lower[j], upper[j]))
        sbplx(theta0, optim, lower = lower, upper = upper, control = control)
      })
    } else {
      cl <- control$cl
      control$cl <- NULL
      res <- parLapply(cl, 1:runs, function(i) {
        theta0 <- sapply(seq_along(lower), function(j) runif(1, lower[j], upper[j]))
        sbplx(theta0, optim, lower = lower, upper = upper, control = control)
      })
    }
  }
  return(res)
}








