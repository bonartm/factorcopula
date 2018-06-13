# factorcopula - r package for high dimensional specification, simulation and estimation of factor copula models

[![Travis-CI Build Status](https://travis-ci.org/bonartm/factorcopula.svg?branch=master)](https://travis-ci.org/bonartm/factorcopula)

## features
- [x] simulation from user specified [factor copula models](http://www.tandfonline.com/doi/full/10.1080/07350015.2015.1062384) (e.g. factors and error terms from the skew t, normal, t distribution)
- [x] estimation of high dimensional dependence structures using the [simulated methods of moments](https://pdfs.semanticscholar.org/cc9f/124d25111430f4f2e977869daef6f403e24a.pdf)
- [x] unrestricted, equi-dependence and block-equi-dependence model specification
- [x] implementation of a [structural break test](http://www.wisostat.uni-koeln.de/sites/statistik/abstracts/Manner_Stark_Wied_2017.pdf)
- [x] estimation of asymptotic variance and confidence intervalls
- [ ] implementation of the J test for overidentifying restrictions
- [ ] implementation of faster C++ code with Rcpp



## installation from github master branch
```R
install.packages("devtools")
devtools::install_github("bonartm/factorcopula")
````

## usage
```R
library(factorcopula)
help(package = "factorcopula")

# define a one factor skew-t copula
t <- 1500
k <- c(1, 1) # all variables in the same groups for an equidependence model
beta <- config_beta(k = k, Z = 1)
Z <- config_factor(rst = list(nu = 1/0.25, lambda = lambda), par = c("lambda"))
eps <- config_error(rt = list(df = 1/0.25))

# define the vector of true parameters
theta0 <- c(beta1 = 1.5, lambda = -0.8)

# generate the copula function and simulate values from the copua model
cop <- fc_create(Z, eps, beta)
U <- cop(theta0, t)


# use some marginal distributions (here normal distribution) to simulate some Y values
Y <- qnorm(U)

# define boundaries for optimzation
lower <- c(beta1 = 0, lambda = -0.9)
upper <- c(beta1 = 5, lambda =  0.9)


# fit the copula 


m <- fc_fit(Y, Z, eps, beta, lower, upper, S = 20000, se = TRUE)
m$theta.second.stage
m$Q
m$message

# plot observed and simulated values
plot(Y, pch = 20)
points(qnorm(cop(m$theta.second.stage, 2000)), col = "red", pch = 20)

# confidence intervalls and p-values
round(m$ci, 4)

````
