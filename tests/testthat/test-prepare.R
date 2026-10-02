test_that("prepare_bugs_data keeps Normal hosts and first-appearance strain order", {
  root <- testthat::test_path("..", "..")
  sys.source(file.path(root, "R", "config.R"), envir = environment())
  sys.source(file.path(root, "R", "prepare_data.R"), envir = environment())

  raw <- data.frame(
    t = c("S.a", "S.b", "S.a", "S.c"),
    S = c("Normal", "Susceptible", "Normal", "Normal"),
    log10dose = c(1, 2, 3, 4),
    N = c(10, 10, 10, 10),
    Y = c(1, 2, 3, 4)
  )
  out <- prepare_bugs_data(raw)
  expect_equal(out$bugs$N_total, 3)
  expect_equal(as.character(out$serovars), c("S.a", "S.c"))
  expect_equal(out$bugs$K, 2)
  expect_equal(out$bugs$d, 10^c(1, 3, 4))
  expect_equal(out$bugs$strain, c(1L, 1L, 2L))
})

test_that("infection_prob is 0 at dose 0 when beta is positive", {
  root <- testthat::test_path("..", "..")
  sys.source(file.path(root, "R", "config.R"), envir = environment())
  sys.source(file.path(root, "R", "prepare_data.R"), envir = environment())
  expect_equal(infection_prob(1, 1, 0), 0)
})
