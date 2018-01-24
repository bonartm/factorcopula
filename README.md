# factorcopula - r package for high dimensional specification, simulation and estimation of factor copula models
_____
## features
- simulation from user specified factor copula models
- estimation using the [simulated methods of moments](https://pdfs.semanticscholar.org/cc9f/124d25111430f4f2e977869daef6f403e24a.pdf)
- implementation of a [strucutral break test](http://www.wisostat.uni-koeln.de/sites/statistik/abstracts/Manner_Stark_Wied_2017.pdf)


## installation
```R
install.packages("devtools")
devtools::install_github("bonartm/factorcopula")
````

## usage
```R
library(factorcopula)
N <- 2
Z <- list(rst = list(nuInv = "nuInv", lambda = "lambda"))
eps <- list(rtInv = (list(dfInv = "nuInv")))
beta <- matrix("beta1", nrow = N)
copFun <- factorCopula(beta, N, Z, eps)
theta <- c(beta1 = 2, nuInv = 0.25, lambda = -0.8)
lower <- c(0, 0.01, -0.99)
upper <- c(10, 0.49, 0.99)
names(lower) <- names(theta)
names(upper) <- names(theta)
U <- copFun(theta, 1000)
Y <- qnorm(U)

m <- fitFactorCopula(Y, copFun, S = 10000, lower = lower, upper = upper, method = "subplex", 
                      control = list(stopval = 0, xtol_rel = 0, maxeval = 2000, ftol_abs = 1e-5, runs = 4))
plot(Y, pch = 20)
points(qnorm(copFun(m$par, 2000)), col = "red", pch = 20)
````
