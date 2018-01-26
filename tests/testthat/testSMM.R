context("SMM")

library(factorcopula)

N <- 2
Z <- list(rst = list(nuInv = "nuInv", lambda = "lambda"))

eps <- list(rtInv = (list(dfInv = "nuInv")))
beta <- matrix("beta1", nrow = N)

cop <- factorCopula(beta, N, Z, eps, zFixed = FALSE)

theta <- c(beta1 = 2.5, nuInv = 0.25, lambda = -0.8)

lower <- c(beta1 = 0, nuInv = 0.01, lambda = -0.99)
upper <- c(beta1 = 10, nuInv = 0.24, lambda = 0.99)

T <- 1000

Y <- qnorm(cop(theta, T))

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "two-stage",
                control = list(stopval = 0, xtol_rel = 1e-12, maxeval = 1000), trials = 2)

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "DEoptim",
                control = list(c = 0.4, itermax = 100, reltol = 1e-6, steptol = 20))

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "genoud",
                control = list(max.generations = 100, wait.generations = 50, solution.tolerance = 1e-8))

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "subplex",
                control = list(stopval = 0, xtol_rel = 1e-12, maxeval = 3000), trials = 2)

cl <- makeCluster(4)
loadPackagesOnCluster(cl, c("nloptr", "factorcopula", "DEoptim", "rgenoud"))
clusterExport(cl, ls())

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "two-stage",
                control = list(stopval = 0, xtol_rel = 1e-12, maxeval = 1000), trials = 4, cl = cl)

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "DEoptim",
                control = list(c = 0.4, itermax = 100, reltol = 1e-6, steptol = 50), cl = cl)

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "genoud",
                control = list(max.generations = 100, wait.generations = 50, solution.tolerance = 1e-8), cl = cl)

fitFactorCopula(Y, cop, S = 2000, lower = lower, upper = upper, method = "subplex",
                control = list(stopval = 0, xtol_rel = 1e-12, maxeval = 3000), trials = 4, cl = cl)

stopCluster(cl)

test_that("fit factor copula", {
  expect_equal(1, 1)
})

