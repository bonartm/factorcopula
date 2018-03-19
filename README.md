# factorcopula - r package for high dimensional specification, simulation and estimation of factor copula models

## features
- [x] simulation from user specified [factor copula models](http://www.tandfonline.com/doi/full/10.1080/07350015.2015.1062384) (e.g. factors and error terms from the skew t, normal, t distribution)
- [x] (parallel) estimation of high dimensional dependence structures using the [simulated methods of moments](https://pdfs.semanticscholar.org/cc9f/124d25111430f4f2e977869daef6f403e24a.pdf)
- [x] unrestricted, equi-dependence and block-equi-dependence model specification

- [x] implementation of a [structural break test](http://www.wisostat.uni-koeln.de/sites/statistik/abstracts/Manner_Stark_Wied_2017.pdf)
- [ ] estimation of asymptotic variance and implementation of the J test for overidentifying restrictions
- [ ]  implement faster C++ code with Rcpp


## installation from github master branch
```R
install.packages("devtools")
devtools::install_github("bonartm/factorcopula")
````

## usage
```R
library(factorcopula)

# define a one factor skew-t copula

t <- 1500
k <- c(1, 1) # all variables in the same groups for an equidependence model
beta <- config_beta(k = k, Z = 1)
Z <- config_factor(rst = list(nu = 1/nuInv, lambda = lambda), par = c("nuInv", "lambda"))
eps <- config_error(rt = list(df = 1/nuInv), par = c("nuInv"))

# define the vector of true parameters
theta0 <- c(beta1 = 2, nuInv = 0.25, lambda = -0.8)

# generate the copula function and simulate values from the copua model
cop <- fc_create(Z, eps, beta)
U <- cop(theta0, t)

# use some marginal distributions (here normal distribution) to simulate some Y values
Y <- qnorm(U)

# define boundaries for optimzation
lower <- c(0, 0.01, -0.99)
upper <- c(5, 0.49, 0.99)
names(lower) <- names(theta0)
names(upper) <- names(theta0)

# fit the copula and plot real together with simulated values
m <- fc_fit(Y, cop, lower = lower, upper = upper, method = "subplex",
         control = list(stopval = 0, xtol_rel = 1e-9, maxeval = 3000), trials = 4, S = 10000, k = k)
plot(Y, pch = 20)
points(qnorm(cop(m$best, 2000)), col = "red", pch = 20)
````
