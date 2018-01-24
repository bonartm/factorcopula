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
