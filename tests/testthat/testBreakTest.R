context("Structural Break Test")

library(factorcopula)
library(parallel)
library(cheopsr)
library(ggplot2)

#cheops_install_github("bonartm/factorcopula", ref = "dev")

options(cheopsr.account = "AG-Wied")
options(cheopsr.username = "bonartm")

N <- 3
k <- c(1, 2, 3)
beta <- config_beta(k, 1)

Z <- config_factor(rst = list(nu = 1/0.25, lambda = -0.8))
eps <- config_error(rt = list(df = 1/0.25))

cop <- fc_create(Z, eps, beta)

theta0 <- c(beta1 = 0.5, beta2 = 1.5, beta3 = 1.5)
theta1 <- c(beta1 = 1.5, beta2 = 1.5, beta3 = 1.5)
lower <- c(beta1 = 0, beta2 = 0, beta3 = 0)
upper <- c(beta1 = 6, beta2 = 6, beta3 = 6)

Y <- qnorm(rbind(cop(theta0, 1000),cop(theta1, 1000)))
brk <- 1000
tSeq <- 300:nrow(Y)

U <- apply(Y, 2, factorcopula:::empDist)
mHat <- factorcopula:::moments(U[1:1000, ], k)
round(mHat - factorcopula:::moments(cop(theta0, 25000, NULL), k), 2)

cl <- makeCluster(4)
cluster_library(cl, "factorcopula")
res <- fc_fit(Y = Y[1:1000, ], Z, eps, beta, lower = lower, upper = upper, k = k, recursive = FALSE,
            control = list(stopval = 0, xtol_rel = 1e-13, maxeval = 1000), S = 25000, cl = cl)
stopCluster(cl)

opt <- cheops_slurmcontrol(nodes = 60, tasks = 1, mem = "2gb", time = "02:00:00")
job <- cheops_run(fc_fit, opt, "re-test2",
           args = list(Y = Y, config_factor = Z, config_error = eps, config_beta = beta,
                       lower = lower, upper = upper, k = k, recursive = TRUE,
                       control = list(stopval = 0, xtol_rel = 1e-10, maxeval = 2000), S = 25000),
           packages = "factorcopula")
cheops_jobs()
cat(cheops_getlog("re-test2"), sep = "\n")
#cheops_cancel("rec-test")
res <- cheops_readRDS('./re-test2/res.rds')
res$p <- fc_pstat(res[,(1:length(lower)), drop = FALSE], res$t)
res$m <- fc_mstat(Y, tSeq, k)

ggplot(res, aes(x = t, y = p)) +
  geom_line() +
  geom_smooth()

