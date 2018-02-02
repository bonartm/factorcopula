context("Structural Break Test")

library(factorcopula)
library(parallel)
library(cheopsr)

options(cheopsr.account = "AG-Wied")
options(cheopsr.username = "bonartm")

N <- 2
Z <- list(rst = list(lambda = -0.7, nuInv = 0.3))
eps <- list(rnorm = list())
beta <- matrix("beta1", nrow = 2)
cop <- factorCopula(beta, N, Z, eps)

theta0 <- c(beta1 = 0.5)
theta1 <- c(beta1 = 3)
lower <- c(beta1 = 0)
upper <- c(beta1 = 10)

Y <- qnorm(rbind(cop(theta0, 700), cop(theta1, 700)))

cop <- factorCopula(beta, N, Z, eps, S = 10000, zFixed = TRUE)

tSeq <- seq(300, nrow(Y), 1)

opt <- cheops_slurmcontrol(nodes = 2, tasks = 8, mem = "2gb", time = "00:30:00", partition = "devel")
cheops_lapply(tSeq, function(t, Y, cop, lower, upper){
  fitFactorCopula(Y[1:t, ], cop, lower = lower, upper = upper, method = "subplex",
                  control = list(stopval = 0, xtol_rel = 1e-7, maxeval = 1000), trials = 4)
}, options = opt, args = list(Y = Y, cop = cop, lower = lower, upper = upper), packages = c("factorcopula"))

res <- do.call(parallel::parLapplyLB, c(list(cl = cl, X = l$x, fun = l$fun), l$args))


cheops_jobs()
cheops_getlog("tmp")

res <- unlist(lapply(res, function(x) x[which.min(x[,"Q"]),"beta1"]))
theta <- matrix(res, ncol = 1)
P <- getPStat(theta, tSeq)
plot(tSeq, P, type = "l")
abline( v = 700)

K <- critVal(Y, B = 10, copFun = cop, theta = theta[length(tSeq), 1] , tSeq = tSeq)
abline(h = quantile(K, 1-0.05))

test_that("derivative", {
  X <- MASS::mvrnorm(1000, c(0, 0), matrix(c(1, 0.8, 0.8, 1), ncol = 2))
  expect_equal(cor(X[,1], X[,2], method = "spearman"), rankCor(empDist(X[,1]), empDist(X[,2])))

  X[,1] <- sample(1:10, 1000, replace = TRUE)
  X[,2] <- sample(1:10, 1000, replace = TRUE)
  expect_equal(cor(X[,1], X[,2], method = "spearman"), rankCor(empDist(X[,1]), empDist(X[,2])))

  x <- sample(1:10, 10, replace = TRUE)
  y <- sample(1:10, 10, replace = TRUE)
  expect_equal(cor(x, y, method = "spearman"), rankCor(empDist(y), empDist(x)))


})

