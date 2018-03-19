# factorcopula - r package for high dimensional specification, simulation and estimation of factor copula models

## features
- [x] simulation from user specified [factor copula models](http://www.tandfonline.com/doi/full/10.1080/07350015.2015.1062384) (e.g. factors and error terms from the skew t, normal, t distribution)
- [x] (parallel) estimation of high dimensional dependence structures using the [simulated methods of moments](https://pdfs.semanticscholar.org/cc9f/124d25111430f4f2e977869daef6f403e24a.pdf)
- [x] implementation of a [structural break test](http://www.wisostat.uni-koeln.de/sites/statistik/abstracts/Manner_Stark_Wied_2017.pdf)
- [ ] estimation of asymptotic variance and implementation of the J test for overidentifying restrictions



## installation from github master branch
```R
install.packages("devtools")
devtools::install_github("bonartm/factorcopula")
````

## usage
```R
library(factorcopula)

# define a one factor skew-t copula
N <- 2
Z <- list(rst = list(nuInv = "nuInv", lambda = "lambda"))
eps <- list(rtInv = (list(dfInv = "nuInv")))
beta <- matrix("beta1", nrow = N)

# define the set of true parameters and lower and upper bounds
theta <- c(beta1 = 2, nuInv = 0.25, lambda = -0.8)
lower <- c(0, 0.01, -0.99)
upper <- c(10, 0.49, 0.99)
names(lower) <- names(theta)
names(upper) <- names(theta)

# generate the copula function and simulate values from the copua model
copFun <- factorCopula(beta, N, Z, eps)
U <- copFun(theta, 1000)

# use some marginal distributions (here normal distribution) to simulate some Y values
Y <- qnorm(U)

# fit the copula and plot real values aggainst fitted values
m <- fitFactorCopula(Y, copFun, S = 10000, lower = lower, upper = upper, method = "subplex", 
                      control = list(stopval = 0, xtol_rel = 0, maxeval = 2000, ftol_abs = 1e-5, runs = 4))
plot(Y, pch = 20)
points(qnorm(copFun(m$par, 2000)), col = "red", pch = 20)
````
