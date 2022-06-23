#context("Utility functions")

n <- 100
x <- runif(n)
y <- runif(n)
fit <- qad(x,y, print = FALSE, resolution = 10)

test_that("coef.qad", {
  expect_s3_class(fit, "qad")
  expect_equal(unname(coef(fit)),c(fit$results$coef, fit$results$p.values))
})


test_that("summary.qad", {
  test <- summary(fit)
  expect_equal(test$SampleSize, n)
  expect_equal(test$resolution, 10)
  expect_equal(test$dependence_values, fit$results)
  expect_output(summary(fit))
})


test_that("plot.qad", {
  expect_equal(class(plot(fit)), c("gg", "ggplot"))
  expect_equal(class(plot(fit, addSample = TRUE)), c("gg", "ggplot"))
  expect_equal(class(plot(fit, copula = TRUE)), c("gg", "ggplot"))
  expect_equal(class(plot(fit, density = TRUE)), c("gg", "ggplot"))
  expect_equal(class(plot(fit, margins = TRUE)), c("gg", "ggplot"))
  expect_equal(class(plot(fit, panel.grid = FALSE)), c("gg", "ggplot"))
})

test_that("plot_density", {
  expect_equal(class(plot_density(fit$mass_matrix)), c("gg", "ggplot"))
})



test_that("pqad", {
  expect_type(pqad(0.2, 100, R = 10, resolution = 10), "double")
  expect_type(qqad(0.2, 100, R = 10, resolution = 10), "double")
})
