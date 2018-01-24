loadPackages <- function(cl, packages){
  clusterExport(cl, "packages", envir = environment())
  res <- clusterEvalQ(cl, invisible(lapply(packages, library, character.only = TRUE, logical.return = TRUE)))
  all(unlist(res))
}



getStandResiduals <- function(x){
  spec = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                    mean.model=list(armaOrder=c(1,0), include.mean=TRUE),
                    distribution.model="norm")
  garch <- ugarchfit(spec = spec, data = x)
  res <- as.vector(residuals(garch, standardize = TRUE))
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


# microbenchmark(rankCor(empDist(X[,1]), empDist(X[,2])), times = 20)
# microbenchmark(cor(X[,1], X[,2], method = "spearman"), times = 20)


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

dependence <- function(u, v, q = c(0.05, 0.10, 0.90, 0.95)){
  rank <- rankCor(u, v)
  quantile <- quantDep(u, v, q)
  c(rank, quantile)
}

# library(MASS)
# X <- mvrnorm(1000, rep(0, 100), diag(100))
# U <- apply(X, 2, empDist)
# microbenchmark(.momentsIntra(U), times = 10, unit = "s")
# microbenchmark(cMomentsIntra(U), times = 10, unit = "s")

momentsIntra <- function(U){
  N <- ncol(U)
  m <- rep(0, 5)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      #cat(i, j, "\n")
      m <- m + dependence(U[,i], U[,j])
    }
  }
 2/(N*(N-1))*m
  #cat("\n")
}

momentsInter <- function(Xdist, k){
  N <- ncol(Xdist)
  m <- rep(0, 5)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (k[i] != k[j]){
        m <- m + dependence(Xdist[,i], Xdist[,j])
        #cat(i, j, "\n")
      }
    }
  }
  1/prod(table(k))*m
  #cat("\n")
}


moments <- function(U, k = rep(1, ncol(U))){
  #k: same length as ncol(U), defines groups for which measures are averaged
  M <- max(k)
  N <- ncol(U)

  if (M == N){# unrestricted model
    dep <- matrix(NA, ncol = 5, nrow = 0.5*N*(N-1))
    ind <- 1
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        dep[ind, ] <- dependence(U[,i], U[,j])
        ind <- ind + 1
      }
    }
    dep <- as.vector(dep)
  } else {# equidependence or block equidependence model
    dep <- matrix(0, ncol = 5, nrow = M)
    for (r in 1:M){
      for (s in r:M){
        #cat(r, s, "\n")
        if (r == s){
          m <- momentsIntra(U[, k == r, drop = FALSE])
          dep[r, ] <- dep[r, ] + m/M
        } else {
          cols <- k %in% c(r, s)
          m <- momentsInter(U[, cols, drop = FALSE], k[cols])
          dep[r, ] <- dep[r, ] + m/M
          dep[s, ] <- dep[s, ] + m/M
        }
      }
    }
    dep <- as.vector(dep/M)
  }
  return(dep)
}

matchArgs <- function(args, theta){
  lapply(args, function(arg){
    ind <- match(arg, names(theta))
    if(is.na(ind))
      arg
    else
     theta[ind]
  })
}

simZMat <- function(Z, theta, S, seed = NULL){
  # Z: named list with random number generators and corresponding additional parameters (empty list if no parameters)
  set.seed(seed)
  vapply(names(Z), function(funName) {
    args <- Z[[funName]]
    args$n <- S
    args <- matchArgs(args, theta)
    do.call(funName, args)
  }, numeric(S))
}

simEpsMat <- function(eps, theta, S, N, seed = NULL){
  # eps: named list of length 1 with random numer generator for the disturbance term and corresponding parameters
  set.seed(seed)
  args <- eps[[1]]
  args$n <- S*N
  args <- matchArgs(args, theta)
  matrix(do.call(names(eps[1]), args), ncol = N)
}

genBetaMat <- function(beta, theta){
  ind <- match(beta, names(theta))
  filter <- !is.na(ind)
  beta[seq_along(beta)[filter]] <- theta[ind[filter]]
  return(matrix(as.numeric(beta), ncol = ncol(beta)))
}

factorCopula <- function(beta, N, Z, eps, zFixed = FALSE, epsFixed = zFixed, S = NULL){
  # Returns a function which can be used to simulate values from a factor copula model
  # beta: a character matrix of size [NxK] indicating the names and position of the beta parameters (see example)
  # Z: a named list of size K and names equal to a random number generator function which takes an argument n and possibly other arguments passed as numeric values or characters (if parameter)
  # eps: a named list of size 1 (see Z)
  # zFixed: if TRUE zMatrix is only simulated once and stays fixed
  # epsFixed: see zFixed
  # value: a function with arguments:
    # theta: a named vector of parameters which is matched to the parameters in Z, eps and beta
    # S: the number of simulations
    # seed: possibly a seed to keep the rng fixed during different simulations
  if (zFixed|epsFixed)
    stopifnot(!is.null(S))
  if (zFixed) {
    zMat <- simZMat(Z, NULL, S, NULL)
  }
  if (epsFixed) {
    epsMat <- simEpsMat(eps, NULL, S, N, NULL)
  }

  function(theta, S = NULL, seed = NULL){
    if (!zFixed)
      zMat <- simZMat(Z, theta, S, seed)
    if (!epsFixed)
    epsMat <- simEpsMat(eps, theta, S, N, seed)
    beta <- genBetaMat(beta, theta)
    X <- zMat%*%t(beta) + epsMat
    apply(X, 2, empDist)
  }
}

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




