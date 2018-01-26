context("Structural Break Test")

library(factorcopula)

Z <- list(rst())


test_that("derivative", {
  X <- MASS::mvrnorm(1000, c(0, 0), matrix(c(1, 0.8, 0.8, 1), ncol = 2))
  expect_equal(cor(X[,1], X[,2], method = "spearman"), rankCor(empDist(X[,1]), empDist(X[,2])))

  X[,1] <- sample(1:10, 1000, replace = TRUE)
  X[,2] <- sample(1:10, 1000, replace = TRUE)
  expect_equal(cor(X[,1], X[,2], method = "spearman"), rankCor(empDist(X[,1]), empDist(X[,2])))

  x <- sample(1:10, 10, replace = TRUE)
  y <- sample(1:10, 10, replace = TRUE)
  expect_equal(cor(x, y, method = "spearman"), rankCor(empDist(y), empDist(x)))


})

