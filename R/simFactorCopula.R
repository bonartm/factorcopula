#' Configurate the latent variables of a factorcopula model
#'
#' @param ... One or more expresssions separated by commas. The name of the expression arguments have to be
#' a valid random number generators, the expressions have to be lists of named unquoted arguments.#'
#' @param par a character vector of parameter names used in the specification of the factor matrix
#'
#' @return A list of \link[rlang]{quosure}s
#' @export
config_factor <- function(..., par = c()){
  factorspec <- rlang::exprs(...)
  fc_check(names(factorspec))
  list(spec = factorspec, par = par, fixed = (length(par) == 0))
}

#' Configurate the error part of a factorcopula model
#'
#' @param ... One named expresssion. The name has to be a
#' valid random number generator, the expression has to be a list of named unquoted arguments.
#'
#' @return A list of \link[rlang]{quosure}s
#' @export
config_error <- function(..., par = c()){
  if (length(rlang::exprs(...)) > 1)
    stop("Only one error function allowed.")
  config_factor(..., par = par)
}

#' Configurate the loadings of a factorcopula model
#'
#' @param type Either unrestrictibe, equidependence or bloc-equidependence
#' @param k Numeric vector of length N with eventually positive increasing integers from 1 to N.
#' @param Z Number of latent variables
#'
#' @return A character matrix of parameters which can be used in \link[factorcopula]{fc_create}
#' @export
config_beta <- function(k, Z = NULL){
  M <- max(k)
  N <- length(k)
  tab <- table(k)
  stopifnot(sum(tab) == N)
  stopifnot(length(tab) == M)
  stopifnot(all(as.numeric(names(tab)) == 1:M))

  if (M == N){# unrestrictive model
    stopifnot(!is.null(Z))
    return(matrix(paste0("beta", 1:(N*Z)), ncol = Z))
  }

  if (all(M == 1)){# equidependence
    stopifnot(!is.null(Z))
    return(matrix(rep(paste0("beta", 1:Z), each = N), ncol = Z, nrow = N))
  }

  # bloc- equidependence
  return(genBetaParMat(k))
}

#' Simulate values from a factor copula model
#'
#' @param beta a character matrix of size [NxK] indicating the names and position of the beta parameters as character strings
#' @param factor a configuration specified by \link[factorcopula]{config_factor}
#' @param eps a configuation specified by \link[factorcopula]{config_error}
#'
#' @return a function which can be used to simulate values from a factor copula model. It has the parameters theta, S and seed
#' @export
fc_create <- function(factor, error, beta){
  force(factor)
  force(error)
  force(beta)

  N <- nrow(beta)
  Z <- length(factor$spec)
  stopifnot(ncol(beta) == Z, is.matrix(beta))

  state <- list(theta = -99, S = -99, seed = -99, zMat = matrix(-99), epsMat = matrix(-99))

  function(theta, S, seed = NULL){
      if(state_changed(state, theta, S, seed, factor$par))
        zMat <- rand_restore(fc_sim, factor$spec, S, theta, seed)
      else
        zMat <- state$zMat

      if(state_changed(state, theta, S, seed, error$par)){
        epsMat <- rand_restore(fc_sim, error$spec, S*N, theta, seed)
        epsMat <- matrix(epsMat, ncol = N)
      } else
        epsMat <- state$epsMat

      betaMat <- eval_beta(beta, theta)

      state <<- list(theta = theta, S = S, seed = seed, zMat = zMat, epsMat = epsMat)

      X <- zMat%*%t(betaMat) + epsMat
      apply(X, 2, empDist)
    }
}


state_changed <- function(state, theta, S, seed, parnames){
  if (is.null(seed) || is.null(state$seed))
    res <- TRUE
  else {
    res <- any(state$seed != seed, state$S != S, length(parnames) != 0, state$theta[parnames] != theta[parnames])
  }
  return(res)
}

rand_restore <- function(fun, ...){
  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed
  else
    oldseed <- NULL

  res <- fun(...)

  if (!is.null(oldseed))
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)
  return(res)
}

fc_sim <- function(config, S, theta, seed){
  vapply(names(config), function(funName) {
    set.seed(seed)
    args <- config[[funName]]
    args <- rlang::eval_tidy(args, as.list(theta))
    args$n <- S
    do.call(funName, args)
  }, numeric(S))
}

eval_beta <- function(beta, theta){
  pos <- theta[beta]
  beta[!is.na(pos)] <- pos[!is.na(pos)]
  matrix(as.numeric(beta), ncol = ncol(beta))
}

fc_check <- function(names){
  if(any(names == ""))
    stop("At least one unnamed factor or error matrix config provided")
  for(name in names){
    if (!exists(name))
      stop("function '", name, "' does not exist")
  }
}

config_fixed <- function(spec){
  for (sp in spec){
    catch <- try(rlang::eval_bare(sp), silent = TRUE)
    if (class(catch) == "try-error"){
      return(FALSE)
    }
  }
  return(TRUE)
}





