context("Structural Break Test")

library(factorcopula)
library(parallel)
library(cheopsr)

#cheops_install_github("bonartm/factorcopula", ref = "dev")

options(cheopsr.account = "AG-Wied")
options(cheopsr.username = "bonartm")

N <- 2
Z <- config_factor(rst = list(nu = 1/0.25, lambda = -0.8))
eps <- config_error(rt = list(df = 1/0.25))
beta <- config_beta("equidependence", N = 2, Z = 1)

cop <- fc_create(Z, eps, beta)

theta0 <- c(beta1 = 0.5)
theta1 <- c(beta1 = 1.5)
lower <- c(beta1 = 0)
upper <- c(beta1 = 5)

Y <- qnorm(rbind(cop(theta0, 1000), cop(theta1, 500)))
brk <- 1000
tSeq <- seq(300, nrow(Y), 1)


opt <- cheops_slurmcontrol(nodes = 60, tasks = 6, mem = "2gb", time = "00:10:00")
id <- cheops_lapply(tSeq, function(t, Y, cop, lower, upper){
  fc_fit(Y[1:t, ], cop, lower = lower, upper = upper, method = "two-stage",
                  control = list(stopval = 0, xtol_rel = 1e-13, maxeval = 1000), trials = 4, S = 25000)
}, options = opt, args = list(Y = Y, cop = cop, lower = lower, upper = upper), packages = c("factorcopula"))

cheops_jobs()
cat(cheops_getlog("tmp"), sep = "\n")
#cheops_cancel(9297698)

res <- cheops_readRDS('./tmp/res.rds')
res <- unlist(lapply(res, function(x) x[which.min(x[,"Q"]),"beta1"]))
theta <- matrix(res, ncol = 1)
P <- fc_pstat(theta, tSeq)
M <- fc_mstat(Y, tSeq)


cl <- makeCluster(4)
clusterExport(cl, ls())
loadPackagesOnCluster(cl, "factorcopula")
thetaFull <- fc_fit(Y, cop, lower = lower, upper = upper, method = "subplex",
       control = list(stopval = 0, xtol_rel = 1e-11, maxeval = 1000), trials = 8, S = 25000, cluster = cl)
thetaFull <- thetaFull[which.min(thetaFull[,"Q"]),"beta1"]
Pk <- fc_critval("copula", Y, 1000, tSeq, cop, theta = thetaFull, cl = NULL)
Mk <- fc_critval("moment", Y, 1000, tSeq, cl = NULL)

stopCluster(cl)


library(ggplot2)
library(ggrepel)
result <- data.frame(t = tSeq, M = M, P = P)

result$labelM[result$t == brk] <- "Real breakpoint"
result$labelM[which.max(result$M)] <- "Estimated breakpoint"

result$labelP[result$t == brk] <- "Real breakpoint"
result$labelP[which.max(result$P)] <- "Estimated breakpoint"

subtitle <- paste0("Illustration of a change of beta from 0.5 to 1.5. ",
                  "The dashed line indicates the estimated critical value. ",
                  "The real breakpoint lies at t = 1000. The estimated at t = ", tSeq[which.max(result$M)], ".")

ggplot(result, aes(x = t, y = M, label = labelM)) +
  geom_line() +
  geom_hline(yintercept = quantile(Mk, 1-0.05), linetype = 2) +
  geom_label_repel(na.rm = TRUE, box.padding = 0.3, point.padding = 1, alpha = 0.8) +
  labs(title = "Moments based structural break test",
       subtitle = subtitle) +
  theme_minimal()

subtitle <- paste0("Illustration of a change of beta from 0.5 to 1.5. ",
                   "The dashed line indicates the estimated critical value. ",
                   "The real breakpoint lies at t = 1000. The estimated at t = ", tSeq[which.max(result$P)], ".")

ggplot(result, aes(x = t, y = P, label = labelP)) +
  geom_line() +
  geom_hline(yintercept = quantile(Pk, 1-0.05), linetype = 2) +
  geom_label_repel(na.rm = TRUE, box.padding = 0.3, point.padding = 1, alpha = 0.8) +
  labs(title = "Copula based structural break test",
       subtitle = subtitle) +
  theme_minimal()


