#' Load packages on a cluster
#'
#' @param cl cluster object created by @references snow::makeCluster
#' @param packages list or vector of package names
#'
#' @return TRUE if all packages could be loaded on all cluster nodes
#' @export
loadPackagesOnCluster <- function(cl, packages){
  snow::clusterExport(cl, "packages", envir = environment())
  res <- snow::clusterEvalQ(cl, invisible(lapply(packages, library, character.only = TRUE, logical.return = TRUE)))
  all(unlist(res))
}
