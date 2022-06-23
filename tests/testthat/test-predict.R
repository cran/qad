#context("Prediction function of qad")

test_that("predictions sum to one", {
  n <- 8
  x <- 1:n
  y <- x^2
  sample <- data.frame(x, y)
  qad.fit <- qad::qad(sample, print = FALSE, resolution = 4)
  pred <-qad:: predict.qad(qad.fit, values = c(1, 3, 5, 7), conditioned = "x1", pred_plot = FALSE)
  expect_equal(pred$prediction$`x1=1`, c(1,0,0,0))
  expect_equal(pred$prediction$`x1=3`, c(0,1,0,0))
  expect_equal(pred$prediction$`x1=5`, c(0,0,1,0))
  expect_equal(pred$prediction$`x1=7`, c(0,0,0,1))

  n <- 250
  x <- runif(n, -1, 1)
  y <- x^2 + runif(n, -0.1,0.1)
  qad.fit <- qad::qad(x,y, print = FALSE)
  pred <- qad::predict.qad(qad.fit, values = c(0, 0.5, -0.5), conditioned = "x1", pred_plot = FALSE)
  expect_equal(sum(pred$prediction$`x1=0`), 1)
  expect_equal(sum(pred$prediction$`x1=0.5`), 1)
  expect_equal(sum(pred$prediction$`x1=-0.5`), 1)
})


test_that("length of prediction intervals is correct", {
  n <- 100
  x <- rnorm(n)
  y <- rnorm(n)
  fit <- qad::qad(x,y, print = FALSE, resolution = 11)
  pred <- qad::predict.qad(fit, values = c(0), conditioned = "x1", pred_plot = FALSE)
  expect_equal(length(pred$prediction$Interval), 11)
})


test_that("plot is returned by prediction", {
  n <- 100
  x <- rnorm(n)
  y <- rnorm(n)
  fit <- qad::qad(x,y, print = FALSE, resolution = 11)
  pred <- qad::predict.qad(fit, values = c(0), conditioned = "x1", pred_plot = FALSE)
  expect_equal(class(pred[[2]]), c("gg", "ggplot"))
})
