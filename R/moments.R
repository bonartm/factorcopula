dependence <- function(u, v, q = c(0.05, 0.10, 0.90, 0.95)){
  rank <- rankCor(u, v)
  quantile <- quantDep(u, v, q)
  c(rank, quantile)
}

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
