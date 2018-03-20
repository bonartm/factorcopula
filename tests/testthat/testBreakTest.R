context("Structural Break Test")

library(factorcopula)
library(parallel)
library(cheopsr)

#cheops_install_github("bonartm/factorcopula", ref = "dev")

options(cheopsr.account = "AG-Wied")
options(cheopsr.username = "bonartm")

N <- 2
k <- c(1, 1)
Z <- config_factor(rst = list(nu = 1/0.25, lambda = -0.8))
eps <- config_error(rt = list(df = 1/0.25))
beta <- config_beta(k, 1)

cop <- fc_create(Z, eps, beta)

theta0 <- c(beta1 = 0.5)
theta1 <- c(beta1 = 1.5)
lower <- c(beta1 = 0)
upper <- c(beta1 = 5)

Y <- qnorm(rbind(cop(theta0, 1000), cop(theta1, 500)))
brk <- 1000

cl <- makeCluster(4)
cluster_library(cl, "factorcopula")
res <- fc_fit(Y, cop, lower, upper, recursive = FALSE, S = 25000, k = k, cl = cl,
       control = list(stopval = 0, xtol_rel = 1e-15, maxeval = 1000))
stopCluster(cl)

opt <- cheops_slurmcontrol(nodes = 40, tasks = 4, mem = "2gb", time = "00:20:00")
job <- cheops_run(fc_fit, opt, "test",
           args = list(Y = Y, copFun = cop, lower = lower, upper = upper, k = k, recursive = TRUE,
                       control = list(stopval = 0, xtol_rel = 1e-13, maxeval = 1000), S = 25000),
           packages = "factorcopula")

cheops_jobs()
cat(cheops_getlog("test"), sep = "\n")
#cheops_cancel(9297698)
res <- cheops_readRDS('./test/res.rds')

