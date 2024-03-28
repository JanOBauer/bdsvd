test_that("block detection works",{
  V <- matrix(c(1,0,
                2,0,
                0,3,
                0,3), 4, 2, byrow = TRUE)

  rownames(V) <- c("A", "B", "C", "D")
  detected.blocks <- detect.blocks(V)

  expect_equal(detected.blocks[[1]]@features, c("A", "B"))
  expect_equal(detected.blocks[[1]]@block.columns, 1)
  expect_equal(detected.blocks[[2]]@features, c("C", "D"))
  expect_equal(detected.blocks[[2]]@block.columns, 2)

  expect_error(detect.blocks(V, 1.5) )
  expect_error(detect.blocks(V, 2.5) )
})



require(mvtnorm)
set.seed(1)
Sigma <- bdsvd.cov.sim(p = 200, b = 2, design = "c")
X <- rmvnorm(500, mean=rep(0, 200), sigma=Sigma)
colnames(X) <- 1:200

ht <- bdsvd.ht(X)
test_that("hyperparameter tuning for BD-SVD works",{
  expect_equal(ht$dof, 100)
})



set.seed(1)
Sigma <- bdsvd.cov.sim(p = 200, b = 2, design = "c")
X <- rmvnorm(500, mean=rep(0, 200), sigma=Sigma)
colnames(X) <- 1:200

bdsvd.obj <- bdsvd(X, standardize = FALSE)
test_that("BD-SVD detects the correct number of blocks",{
  expect_equal(length(bdsvd.obj), 2)
  expect_equal(length(bdsvd.obj[[1]]), 100)
  expect_equal(length(bdsvd.obj[[2]]), 100)
})



set.seed(1)
Sigma <- bdsvd.cov.sim(p = 200, b = 4, design = "c")
X <- rmvnorm(500, mean=rep(0, 200), sigma=Sigma)
colnames(X) <- 1:200

bdsvd.obj <- bdsvd(X, standardize = FALSE)
test_that("BD-SVD detects the correct number of blocks",{
  expect_equal(length(bdsvd.obj), 4)
  expect_equal(length(bdsvd.obj[[1]]), 50)
  expect_equal(length(bdsvd.obj[[2]]), 50)
  expect_equal(length(bdsvd.obj[[3]]), 50)
  expect_equal(length(bdsvd.obj[[4]]), 50)
})



set.seed(1)
Sigma <- bdsvd.cov.sim(p = 200, b = 2, design = "d")
X <- rmvnorm(500, mean=rep(0, 200), sigma=Sigma)
colnames(X) <- 1:200

bdsvd.obj <- bdsvd(X, standardize = FALSE)
test_that("BD-SVD detects the correct number of blocks",{
  expect_equal(length(bdsvd.obj), 2)
  expect_equal(length(bdsvd.obj[[1]]), 100)
  expect_equal(length(bdsvd.obj[[2]]), 100)
})



set.seed(1)
Sigma <- bdsvd.cov.sim(p = 200, b = 4, design = "d")
X <- rmvnorm(500, mean=rep(0, 200), sigma=Sigma)
colnames(X) <- 1:200

bdsvd.obj <- bdsvd(X, standardize = FALSE)
test_that("BD-SVD detects the correct number of blocks",{
  expect_equal(length(bdsvd.obj), 4)
  expect_equal(length(bdsvd.obj[[1]]), 50)
  expect_equal(length(bdsvd.obj[[2]]), 50)
  expect_equal(length(bdsvd.obj[[3]]), 50)
  expect_equal(length(bdsvd.obj[[4]]), 50)
})
