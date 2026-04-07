library(testthat)
library(wavFeatExt)

test_that("simulateCNA returns correct structure", {
  sim <- simulateCNA(
    n.obs = 20,
    p = 32,
    n.sim = 1,
    n.block = 8,
    verbose = FALSE
  )
  
  expect_true(is.list(sim))
  expect_length(sim, 1)
  expect_true(is.matrix(sim[[1]]))
})

test_that("wavFeatExt works correctly", {
  sim <- simulateCNA(
    n.obs = 20,
    p = 32,
    n.sim = 1,
    n.block = 8,
    verbose = FALSE
  )
  
  res <- wavFeatExt(sim, type = "detail")
  
  expect_true(is.list(res))
  expect_true(is.list(res[[1]]))
})

test_that("classification runs without error", {
  sim <- simulateCNA(
    n.obs = 20,
    p = 32,
    n.sim = 1,
    n.block = 8,
    verbose = FALSE
  )
  
  det <- wavFeatExt(sim, type = "detail")
  sca <- wavFeatExt(sim, type = "scaling")
  
  y <- factor(rep(c("Group1", "Group2"), each = 10))
  
  model <- classifyWavFeatExt(
    sim,
    y,
    det,
    sca,
    method = "KNN",
    k = 5,
    ite = 1
  )
  
  expect_s3_class(model, "wavFeatExtClassifier")
})