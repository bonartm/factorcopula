checkNamespace <- function(name){
  if (!requireNamespace(name, quietly = TRUE)) {
    stop(paste(name, "needed for this function to work. Please install it.", call. = FALSE))
  }
}

parallelLapply <- function(x, fun, cl, load.balancing = TRUE, ...){
  if(is.null(cl)){
    lapply(x, fun, ...)
  } else {
    if (load.balancing){
      snow::clusterApplyLB(cl, x, fun, ...)
    } else {
      snow::clusterApply(cl, x, fun, ...)
    }
  }
}

# getStandResiduals <- function(x){
#   checkNamespace("rugarch")
#   spec = rugarch::ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
#                              mean.model=list(armaOrder=c(1,0), include.mean=TRUE),
#                              distribution.model="norm")
#   garch <- rugarch::ugarchfit(spec = spec, data = x)
#   res <- as.vector(rugarch::residuals(garch, standardize = TRUE))
#   return(res)
# }

#' Load packages on a cluster
#'
#' @param cl cluster object created by \link[snow]{makeCluster}
#' @param packages list or vector of package names
#'
#' @return TRUE if all packages could be loaded on all cluster nodes
#' @export
cluster_library <- function(cl, packages){
  snow::clusterExport(cl, "packages", envir = environment())
  res <- snow::clusterEvalQ(cl, invisible(lapply(packages, library, character.only = TRUE, logical.return = TRUE)))
  all(unlist(res))
}

#' Generate random numbers from the skewed t distribution
#'
#' @param n number of observations
#' @param nu degree of freedoms
#' @param lambda skewness factor between (-1:1)
#'
#' @return a vactor of length n
#' @export
rst <- function(n, nu = 1e9, lambda = 0){
  # Generates random numbers from the skewed t distribution by Hansen (1994)
  # nuInv: Inverse nu parameter (degrees of freedom, q) (0, )
  # lambda: skewness parameter (-1, 1)
  stopifnot(nu >= 2 & lambda > -1 & lambda < 1)

  u <- stats::runif(n)
  if (is.infinite(suppressWarnings(gamma(nu/2)))){
    c <- 0
    a <- 0
  } else {
    c <- gamma((nu+1)/2)/(sqrt(pi*(nu-2))*gamma(nu/2))
    a <- 4*lambda*c*((nu-2)/(nu-1));
  }

  b <- sqrt(1 + 3*lambda^2 - a^2);

  f1 <- u < (1-lambda)/2

  inv1 <- (1-lambda)/b*sqrt((nu-2)/nu)*stats::qt(u[f1]/(1-lambda), nu)-a/b
  inv2 <- (1+lambda)/b*sqrt((nu-2)/nu)*stats::qt(0.5+1/(1+lambda)*(u[!f1]-(1-lambda)/2), nu)-a/b

  inv <- rep(0, n)
  inv[f1]  <- inv1
  inv[!f1]  <- inv2
  return(inv)
}







