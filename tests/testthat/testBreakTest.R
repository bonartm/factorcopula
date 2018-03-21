context("Structural Break Test")

library(factorcopula)
library(parallel)
library(cheopsr)

#cheops_install_github("bonartm/factorcopula", ref = "dev")

options(cheopsr.account = "AG-Wied")
options(cheopsr.username = "bonartm")

N <- 6
k <- c(1, 1, 1, 2, 2, 2)
beta <- config_beta(k)

Z <- config_factor(rst = list(nu = 1/0.25, lambda = -0.8), rt = list(df = 1/0.25), rt = list(df = 1/0.25))
eps <- config_error(rt = list(df = 1/0.25))

cop <- fc_create(Z, eps, beta)

theta0 <- c(beta1 = 0.5, beta2 = 1, beta3 = 1.5, beta4 = 2)
theta1 <- c(beta1 = 1.5, beta2 = 1.5, beta3 = 1.5, beta4 = 1.5)
lower <- c(beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0)
upper <- c(beta1 = 6, beta2 = 6, beta3 = 6, beta4 = 6)

Y <- qnorm(rbind(cop(theta0, 2000), cop(theta1, 1000)))
brk <- 2000
tSeq <- 300:nrow(Y)

opt <- cheops_slurmcontrol(nodes = 40, tasks = 4, mem = "2gb", time = "01:00:00")
job <- cheops_run(fc_fit, opt, "rec-test",
           args = list(Y = Y, copFun = cop, lower = lower, upper = upper, k = k, recursive = TRUE,
                       control = list(stopval = 0, xtol_rel = 1e-13, maxeval = 1000), S = 25000),
           packages = "factorcopula")
cheops_jobs()
cat(cheops_getlog("rec-test"), sep = "\n")
#cheops_cancel(9297698)
res <- cheops_readRDS('./rec-test/res.rds')
p <- fc_pstat(res[,-1], res$t)
plot(res$t, p, type = "l")

m <- fc_mstat(Y, tSeq, k)
plot(tSeq, m, type = "l")
abline(v = brk)
