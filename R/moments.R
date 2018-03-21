dependence <- function(u, v, q = c(0.05, 0.10, 0.90, 0.95)){
  rank <- rankCor(u, v)
  quantile <- quantDep(u, v, q)
  c(rank, quantile)
}

.moments <- function(U, index){
  dep <- apply(index, 1, function(x) {
    dependence(U[,x[1]], U[, x[2]])
  })
  rowMeans(dep)
}


moments <- function(U, k = rep(1, ncol(U))){
  #k: same length as ncol(U), defines groups for which measures are averaged
  M <- max(k)
  N <- ncol(U)

  elem <- t(combn(1:N, 2))
  elem <- cbind(elem, k[elem[,1]], k[elem[,2]])

  colnames(elem) <- c("i", "j", "ki", "kj")


  if (M == N){# unrestrictive model
    moments <- matrix(NA, ncol = 5, nrow = 0.5*N*(N-1))
    i <- 1
    for (r in 1:(M-1)){
      for (s in (r+1):M){
        moments[i, ] <- .moments(U, matrix(c(r, s), ncol = 2))
        i <- i+1
      }
    }
    moments <- as.vector(moments)

  } else {# restrictive model: equidependence or bloc-equidependence
    moments <- matrix(0, ncol = 5, nrow = M)
    for (r in 1:M){
      for (s in r:M){
        m <- .moments(U, elem[elem[,3] == r & elem[,4] == s, 1:2, drop = FALSE])
        m <- m/M
        if (r == s){#intra dependence, add only to one group
          moments[r, ] <- moments[r, ] + m
        } else {# inter dependence, add to both groups
          moments[r, ] <- moments[r, ] + m
          moments[s, ] <- moments[s, ] + m
        }
      }
    }
  }
  return(c(t(moments)))
}
