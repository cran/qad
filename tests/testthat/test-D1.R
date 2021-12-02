context("D1 distance")

test_that("D1 distance works correctly", {
  n <- 100
  x1 <- 1:n
  y1 <- x1
  x2 <- 1:n
  y2 <- 1-x2
  expect_equal(qad::D1(x1,y1,x2,y2), 0.5)
  expect_equal(qad::D1(x1,y1,x2,y2, resolution = 1), 0)
  expect_equal(qad::D1(x1,y1,x2,y2, resolution = 100), 0.5)
})


test_that("D1 distance for matrices", {
  x1 <- runif(100)
  y1 <- runif(100)
  x2 <- runif(60)
  y2 <- runif(60)
  C1 <- qad::ECBC(x1,y1)
  C2 <- qad::ECBC(x2,y2)
  expect_warning(qad::D1.ECBC(C1,C2), "The dimensions of the checkerboard copulas have to coincide")
  C1 <- qad::ECBC(x1,y1, resolution = 1)
  C2 <- qad::ECBC(x2,y2, resolution = 1)
  expect_equal(qad::D1.ECBC(C1,C2), 0)
})




