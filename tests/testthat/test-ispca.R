test_that("ispca works",{
  set.seed(1)
  Sigma <- bdsvd.cov.sim(p = 100, b = 4, design = "c")
  Sigma <- Sigma[1:80, 1:80]
  p <- ncol(Sigma)
  X <- mvtnorm::rmvnorm(500, mean = rep(0, p), sigma = Sigma)
  colnames(X) <- seq_len(p)


  ispca.obj <- ispca(X)
  expect_equal(length(ispca.obj), 5)
  expect_equal(length(ispca.obj$l), 4)
  expect_equal(nrow(ispca.obj$X.b[[1]]), 500)


  bdsvd.obj <- bdsvd(X, anp = "1", standardize = FALSE)
  spca.obj <- ispca(X, block.structure = bdsvd.obj)
  expect_equal(length(ispca.obj), 5)
  expect_equal(length(ispca.obj$l), 4)
  expect_equal(nrow(ispca.obj$X.b[[1]]), 500)


  #Wrong class:
  expect_error(ispca(X, block.structure = bdsvd.obj[[1]]))

  #From wrong data set. Warning and error is expected
  block.structure <- detect.blocks(Sigma[5:80, 5:80])
  expect_error(suppressWarnings(ispca(X, block.structure = block.structure)))
  expect_warning(try(ispca(X, block.structure = block.structure)))

  #Variables do not match. Warning and error is expected
  block.structure <- detect.blocks(Sigma)
  block.structure[[1]]$features <- as.character(101:125)
  expect_error(suppressWarnings(ispca(X, block.structure = block.structure)))
  expect_warning(try(ispca(X, block.structure = block.structure)))
})



test_that("prmats works",{
  set.seed(1)
  Sigma <- bdsvd.cov.sim(p = 100, b = 4, design = "c")
  Sigma <- Sigma[1:80, 1:80]
  p <- ncol(Sigma)
  X <- mvtnorm::rmvnorm(500, mean = rep(0, p), sigma = Sigma)
  colnames(X) <- seq_len(p)
  bdsvd.obj <- bdsvd(X, anp = "1", standardize = FALSE)


  prmats.obj <- prmats(X, bdsvd.obj, value = 1)
  expect_equal(length(prmats.obj$prmats), 4)
  expect_equal(sum(vapply(prmats.obj$prmats, `[[`, numeric(1), 1)), 1)
  expect_true(inherits(prmats.obj, "prmats"))

  prmats.obj <- prmats(X, bdsvd.obj, value = 0.6)
  expect_equal(length(prmats.obj$prmats), 2)

  prmats.obj <- prmats(X, bdsvd.obj, rule = "enrich", value = 1)
  expect_equal(length(prmats.obj$prmats), 2)

  expect_warning(prmats.obj <- prmats(X, bdsvd.obj, rule = "enrich", value = 2))
  expect_equal(prmats.obj$prmats[[1]], NULL)
  expect_equal(prmats.obj$prmats$expl.var, NULL)
})


