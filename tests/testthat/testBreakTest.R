context("Structural Break Test")

library(factorcopula)
library(parallel)
library(cheopsr)
library(ggplot2)

#cheops_install_github("bonartm/factorcopula", ref = "dev")

options(cheopsr.account = "AG-Wied")
options(cheopsr.username = "bonartm")

N <- 6
k <- c(1, 1, 1, 2, 2, 2)
beta <- config_beta(k)

Z <- config_factor(rst = list(nu = 1/0.25, lambda = -0.8), rt = list(df = 1/0.25), rt = list(df = 1/0.25))
eps <- config_error(rt = list(df = 1/0.25))

cop <- fc_create(Z, eps, beta)

theta0 <- c(beta1 = 1.5, beta2 = 0, beta3 = 1.5, beta4 = 0)
theta1 <- c(beta1 = 1.5, beta2 = 1.5, beta3 = 1.5, beta4 = 1.5)
lower <- c(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0)
upper <- c(beta1 = 6, beta2 = 6, beta3 = 6, beta4 = 6)

cop(theta0, 10)
cop(theta0, 10)

Y <- qnorm(rbind(fc_create(Z, eps, beta)(theta0, 1000),fc_create(Z, eps, beta)(theta1, 1000)))
brk <- 1000
tSeq <- 300:nrow(Y)

# U <- apply(Y, 2, factorcopula:::empDist)
# mHat <- factorcopula:::moments(U[1:1000, ], k)
# round(mHat - factorcopula:::moments(cop(c(beta1 = 2, beta2 = -1, beta3 = 1.5, beta4 = 0, lambda = -0.9), 25000, 2), k), 4)

cl <- makeCluster(4)
cluster_library(cl, "factorcopula")
fc_fit(Y = Y[1:1000, ], Z, eps, beta, lower = lower, upper = upper, k = k, recursive = FALSE,
            control = list(stopval = 0, xtol_rel = 1e-9, maxeval = 1000), S = 2500, cl = cl)
stopCluster(cl)

opt <- cheops_slurmcontrol(nodes = 60, tasks = 1, mem = "2gb", time = "01:00:00")
job <- cheops_run(fc_fit, opt, "re-test2",
           args = list(Y = Y, copFun = cop, lower = lower, upper = upper, k = k, recursive = TRUE,
                       control = list(stopval = 0, xtol_rel = 1e-14, maxeval = 3000), S = 25000),
           packages = "factorcopula")
cheops_jobs()
cat(cheops_getlog("re-test2"), sep = "\n")
#cheops_cancel("rec-test")
res <- cheops_readRDS('./re-test2/res.rds')
res$p <- fc_pstat(res[,c(2, 4)], res$t)


ggplot(res, aes(x = t, y = p)) +
  geom_line() +
  geom_smooth()

m <- fc_mstat(Y, tSeq, k)
plot(tSeq, m, type = "l")
abline(v = brk)
