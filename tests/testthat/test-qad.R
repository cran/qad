#context("qad and zeta1")

test_that("zeta1 returns exact values", {
  #Example 1
  x <- c(1,2,3,4,5)
  y <- c(1,2,3,4,5)
  expect_equal(qad::zeta1(x, y, resolution = 5), 0.9)
  expect_equal(qad::zeta1(x, y, resolution = 2), 0.6)
  expect_equal(qad::zeta1(x, y, resolution = 1), 0)
  expect_equal(qad::zeta1(x, y, resolution = 10), 0.9)
  expect_equal(qad::zeta1(x,y), qad::zeta1(y,x))

  #Example 2
  x <- c(1,2,3,4,5,6,7,8,9,10)
  y <- c(1,1,1,1,3,3,1,1,1,1)
  expect_equal(qad::zeta1(x, y, resolution = 10), 0.48)
  expect_equal(qad::zeta1(y, x, resolution = 10), 0.24)

  #Example (checkmin copula)
  NN <- c(2,10,22,44,122,200)
  for(N in NN){
    x <- 1:N
    y <- 1:N
    expect_equal(qad::zeta1(x,y, resolution = N), 1-1/(2*N))
    expect_equal(qad::zeta1(y,x, resolution = N), 1-1/(2*N))
  }
})


test_that("class of qad output and correct data output", {
  n <- 20
  x <- runif(n)
  y <- runif(n)
  expect_output(qad::qad(x,y))
  fit <- qad::qad(x,y, print = FALSE)
  expect_identical(class(fit), "qad")
  expect_equal(fit$n, n)
  expect_equal(fit$resolution, floor(sqrt(n)))
  expect_equal(colSums(fit$mass_matrix), rowSums(fit$mass_matrix))
  expect_equal(NROW(fit$data), n)
  expect_equal(NCOL(fit$data), 2)
})




n <- 100
x <- runif(n)
y <- runif(n)
z <- runif(n)
df <- data.frame(x,y,z)
pw_fit <- pairwise.qad(df)

test_that("pairwise qad", {
  expect_equal(class(pw_fit), "list")
  expect_equal(names(pw_fit), c("q",
                                "max.dependence",
                                "asymmetry",
                                "q_p.values",
                                "q_p.values.adjusted",
                                "max.dependence_p.values",
                                "max.dependence_p.values.adjusted",
                                "asymmetry_p.values",
                                "asymmetry_p.values.adjusted",
                                "resolution",
                                "n_removed_00",
                                "p.adjust.method")
)
  expect_equal(class(pw_fit$q), "data.frame")
  expect_equal(class(pw_fit$max.dependence), "data.frame")
  expect_equal(class(pw_fit$asymmetry), "data.frame")
  expect_equal(class(pw_fit$q_p.values), "data.frame")
  expect_equal(class(pw_fit$max.dependence_p.values), "data.frame")
  expect_equal(class(pw_fit$asymmetry_p.values), "data.frame")
  expect_equal(class(pw_fit$resolution), "data.frame")
})

test_that("heatmap.qad", {
  expect_equal(class(qad::heatmap.qad(pw_fit, select = "dependence")), c("ggplot2::ggplot","ggplot","ggplot2::gg","S7_object","gg" ))
  expect_equal(class(qad::heatmap.qad(pw_fit, select = "dependence",
                                      fontsize = 5, significance = T, scale = "rel", color = "viridis", title = "test")), 
               c("ggplot2::ggplot","ggplot","ggplot2::gg","S7_object","gg" ))
})



