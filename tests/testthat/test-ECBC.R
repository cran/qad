context("Empirical checkerboard copula")

test_that("ECBC returns a matrix of correct dimension", {
  x <- runif(10)
  y <- runif(10)
  expect_type(qad::ECBC(x,y), "double")
  expect_equal(dim(qad::ECBC(x,y)), c(3,3))
  expect_equal(dim(qad::ECBC(x,y, resolution = 2)), c(2,2))
  expect_equal(dim(qad::ECBC(x,y, resolution = 1)), c(1,1))
  expect_equal(dim(qad::ECBC(x,y, resolution = 10)), c(10, 10))
  expect_warning(qad::ECBC(x,y, resolution = 11), "Resolution cannot exceed sample size")
})


test_that("ECBC returns uniform margins", {
  x <- runif(100)
  y <- runif(100)
  expect_equal(colSums(qad::ECBC(x,y)), rep(1/10, 10))
  expect_equal(rowSums(qad::ECBC(x,y)), rep(1/10, 10))
  expect_equal(colSums(qad::ECBC(x,y, resolution = 4)), rep(1/4, 4))
  expect_equal(rowSums(qad::ECBC(x,y, resolution = 4)), rep(1/4, 4))
  expect_equal(colSums(qad::ECBC(x,y, resolution = 100)), rep(1/100, 100))
  expect_equal(rowSums(qad::ECBC(x,y, resolution = 100)), rep(1/100, 100))

  x <- c(1,1,1,2,2,3,3,3)
  y <- c(1,2,3,4,5,6,6,7)
  expect_equal(rowSums(qad::ECBC(x,y, resolution = 3)), rep(1/3, 3))
  expect_equal(colSums(qad::ECBC(x,y, resolution = 3)), rep(1/3, 3))
})


test_that("ECBC.eval returns correct results", {
  #Diagonal matrix (Minimum copula)
  M <- diag(1/5, 5)
  expect_equal(qad::ECBC.eval(M, eval.points = cbind(c(0,0.2,0.3,0.5,0.7,1), c(0.3,0.5,0.6,0.8,0,1))),
               c(0,0.2,0.3,0.5,0,1))
})



