#' @importFrom methods new
create.block <- function(feature.names, selected.features, block.columns) {
  if (length(feature.names) > 0) {
    selected.features <- feature.names[selected.features]
  }
  return(new("block", features = selected.features, block.columns = block.columns))
}

get.blocks <- function(threshold.matrix, feature.names) {
  p <- nrow(threshold.matrix)
  k <- ncol(threshold.matrix)
  columns <- 1:k
  blocks <- list()

  if (any(colSums(threshold.matrix == 0) == p)) {
    stop("Zero column. Reduce threshold value.\n")
  }

  if (any(rowSums(threshold.matrix == 0) == k)) {
    if (k < p) {
      stop("Zero row. Add more loadings, or reduce threshold value.\n")
    } else {
      stop("Zero row. Reduce threshold value.\n")
    }

  }

  while (length(columns) != 0) {
    block.columns <- columns[1]
    while (!identical(block.columns, find.blocks(threshold.matrix, block.columns))) {
      block.columns <- find.blocks(threshold.matrix, block.columns)
    }

    if (length(block.columns) == 1) {
      block.idx <- which(threshold.matrix[, block.columns] != 0)
    } else {
      block.idx <- which(rowSums(threshold.matrix[, block.columns] != 0) > 0)
    }

    blocks[[length(blocks) + 1]] <- create.block(
      feature.names = feature.names,
      selected.features = block.idx,
      block.columns = block.columns
    )

    columns <- columns[!columns %in% block.columns]
  }

  return(blocks)
}

get.threshold.matrix <- function(loadings, threshold) {
  loadings[which(abs(loadings) <= threshold)] <- 0
  loadings[which(abs(loadings) != 0)] <- 1

  return(loadings)
}

find.blocks <- function(matrix, column.idx) {
  if (length(column.idx) == 1) {
    rows <- which(matrix[, column.idx] != 0)
  } else {
    rows <- which(rowSums(matrix[, column.idx] != 0) > 0)
  }

  if (length(rows) == 1) {
    cols <- which(matrix[rows, ] != 0)
  } else {
    cols <- which(colSums(matrix[rows, ] != 0) > 0)
  }

  return(cols)
}



#' @title Block
#'
#' @description Class used within the package to store the structure and
#' information about the detected blocks.
#' @slot features numeric vector that contains the the variables
#' corresponding to this block.
#' @slot block.columns numeric vector that contains the indices of the
#' singular vectors corresponding to this block.
#' @export
setClass("block", slots = c(features = "vector", block.columns = "vector"))



#' @title Block Detection Using Singular Vectors (BD-SVD).
#'
#' @description Performs BD-SVD iteratively to reveal the block structure. Splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (which is a sparse approximation of the
#' first right singular vector, i.e., a vector with many zero values) that mirrors the shape of the covariance matrix. This
#' procedure is continued iteratively until the block diagonal structure has been revealed.
#'
#' The data matrix ordered according to this revealed block diagonal structure can be obtained by \link{bdsvd.structure}.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param dof.lim Interval limits for the number of non-zero components in the sparse loading (degrees of freedom).
#' If \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom. Default is \code{dof.lim <- c(0, p-1)} which is highly recommended to check for
#' all levels of sparsity.
#'
#' @param anp Which regularization function should be used for the HBIC. \code{anp = "1"} implements \eqn{a_{np} = 1} which corresponds
#' to the BIC, \code{anp = "2"} implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (202Xa), and \code{anp = "3"}
#' implements \eqn{a_{np} = log(log(np))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#'
#'
#'@param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @param trace Print out progress as iterations are performed. Default is \code{TRUE}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package.
#'
#' @return
#' A list containing the feature names of the submatrices of \code{X}. The length of the list equals
#' the number of submatrices.
#'
#' @seealso \link{bdsvd.structure}, \link{bdsvd.ht}, \link{single.bdsvd}
#'
#' @references \cite{Bauer, J.O. (202Xa). High-dimensional block diagonal covariance structure detection using singular vectors.}
#' @references \cite{Wang, H., B. Li, and C. Leng (2009). Shrinkage tuning parameter selection with a diverging number of parameters, J. R. Stat. Soc. B 71 (3), 671–683.}
#' @references \cite{Wang, L., Y. Kim, and R. Li (2013). Calibrating nonconvex penalized regression in ultra-high dimension, Ann. Stat. 41 (5), 2505–2536.}
#'
#' @examples
#' #Replicate the simulation study (c) from Bauer (202Xa).
#'
#' \dontrun{
#' p <- 500 #Number of variables
#' n <- 250 #Number of observations
#' b <- 10  #Number of blocks
#' design <- "c" #Simulation design "a", "b", "c", or "d".
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
#' colnames(X) <- seq_len(p)
#'
#' bdsvd(X, standardize = FALSE)
#' }
#'
#' @importFrom irlba ssvd
#'
#' @export
bdsvd <- function(X,
                  dof.lim,
                  anp = "2",
                  standardize = TRUE,
                  max.iter,
                  trace = FALSE
) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA)")
    }

  p <- ncol(X)

  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(seq_len(p))
    }

  if (missing(dof.lim)) {
    dof.lim <- c(0, p - 1)
    }
  dof.lim <- sort(round(dof.lim))
  if (dof.lim[1] < 0) {
    dof.lim[1] <- 0
    }
  if (dof.lim[2] > p - 1) {
    dof.lim[2] <- p - 1
    }

  ANP <- c("1", "2", "3")
  if (!(anp %in% ANP))
    stop(paste(anp), " is an invalid option for anp")

  if (missing(max.iter)) {
    max.iter <- 500
    }

  sub.matrices <- list(colnames(X))
  results <- list()

  b <- 1
  while (TRUE) {
    if (length(sub.matrices) == 0) {
      break  #Stop when splitting is no longer supported
    }

    if (length(sub.matrices[[1]]) == 1) {
      if (trace) {
        cat("Block", b, ":", sub.matrices[[1]], "\n", sep = "\t")
      }
      b <- b + 1
      results <- c(results, sub.matrices[[1]])
      sub.matrices <- sub.matrices[-1]
      next
    }

    dof.split <- bdsvd.ht(X = X[, colnames(X) %in% sub.matrices[[1]]],
                          dof.lim = dof.lim,
                          standardize = standardize,
                          anp = anp,
                          max.iter = max.iter)$dof
    sub.results <- single.bdsvd(X = X[, colnames(X) %in% sub.matrices[[1]]], standardize = standardize, dof = dof.split)


    if (length(sub.matrices) == 1) {
      sub.matrices <- list()
    } else {
      sub.matrices <- sub.matrices[-1]
    }

    if (length(sub.results) == 1) {
      if (trace) {
        cat("Block", b, ":", sub.results[[1]], "\n", sep = "\t")
      }
      b <- b + 1
      results <- c(results, sub.results)
    } else {
      sub.matrices <- c(sub.matrices, sub.results)
    }

  }

  results <- results[order(sapply(results, length), decreasing = TRUE)]
  class(results) <- "bdsvd"
  return(results)

}



#' @title Covariance Matrix Simulation for BD-SVD
#'
#' @description This function generates covariance matrices based on the simulation studies described in Bauer (202Xa).
#'
#' @param p Number of variables.
#'
#' @param b Number of blocks. Only required for simulation design "c" and "d".
#'
#' @param design Simulation design "a", "b", "c", or "d".
#'
#' @return
#' A covariance matrix according to the chosen simulation design.
#'
#' @references \cite{Bauer, J.O. (202Xa). High-dimensional block diagonal covariance structure detection using singular vectors.}
#'
#' @examples
#' #The covariance matrix for simulation design (a) is given by
#' Sigma <- bdsvd.cov.sim(p = 500, b = 500, design = "a")
#'
#' @importFrom stats runif
#'
#' @export
bdsvd.cov.sim <- function(p = p,
                    b,
                    design = design
                    ) {

  DESIGN <- c("a", "b", "c", "d")
  if (!(design %in% DESIGN))
    stop(paste(design), " is an invalid design")

  if (design == "a") {
    Sigma <- diag(1, p, p)
    return(Sigma)
  }

  if (design == "b") {
    Sigma <- diag(stats::runif(p, 1, 5), p, p)
    return(Sigma)
  }

  if (design == "c") {
    if (missing(b))
      stop("Number of blocks b required for simulation design c.")

    if (!(p %% b == 0))
      stop("p must be divisible by b so that blocks of equal size can be created.")

    rho0 <- 0.2
    epsilon <- 0.1

    d <- p / b
    Sigma <- matrix(0, p, p)
    for (i in 1:b) {
      rho <- stats::runif(1, rho0 - epsilon, rho0 + epsilon)
      Sigma[(1 + d * (i - 1)) : (d + d * (i - 1)), (1 + d * (i - 1)) : (d + d * (i - 1))] <-
        ((1 - rho) * diag(1, d, d) + 2 * rho * rep(1, d) %*% t(rep(1, d)))
    }
    return(Sigma)
  }

  if (design == "d") {
    if (missing(b))
      stop("Number of blocks b required for simulation design d.")
    if (!(p %% b == 0))
      stop("p must be divisible by b so that blocks of equal size can be created.")

    rho0 <- 0.45
    epsilon <- 0.15
    omega <- 0.1

    d <- p / b
    Sigma <- matrix(0, p, p)
    for (B in 1:b) {
      rho <- stats::runif(1, rho0 - epsilon, rho0 + epsilon)
      Rii <- matrix(0, d, d)
      for (i in 1:d) {
        for (j in 1:d) {
          Rii[i, j] <- (-1)^(i + j) * rho^(abs(i - j)^(omega))
        }
      }
      diag(Rii) <- 1
      Sigma[(1 + d * (B - 1)) : (d + d * (B - 1)), (1 + d * (B - 1)):(d + d * (B - 1))] <- Rii
    }
    return(Sigma)
  }

}



#' @title Hyperparameter Tuning for BD-SVD
#'
#' @description Finds the number of non-zero elements of the sparse loading according to the high-dimensional
#' Bayesian information criterion (HBIC).
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param dof.lim Interval limits for the number of non-zero components in the sparse loading (degrees of freedom).
#' If \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom. Default is \code{dof.lim <- c(0, p-1)} which is highly recommended to check for
#' all levels of sparsity.
#'
#' @param anp Which regularization function should be used for the HBIC. \code{anp = "1"} implements \eqn{a_{np} = 1} which corresponds
#' to the BIC, \code{anp = "2"} implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (202Xa), and \code{anp = "3"}
#' implements \eqn{a_{np} = log(log(np))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package. The computation of the HBIC is outlined in Bauer (202Xa).
#'
#' @return
#' \item{dof}{
#'   The optimal number of nonzero components (degrees of freedom) according to the HBIC.
#' }
#' \item{BIC}{
#'   The HBIC for the different numbers of nonzero components.
#' }
#'
#' @seealso \link{bdsvd}, \link{single.bdsvd}
#'
#' @references \cite{Bauer, J.O. (202Xa). High-dimensional block diagonal covariance structure detection using singular vectors.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#' @references \cite{Wang, H., B. Li, and C. Leng (2009). Shrinkage tuning parameter selection with a diverging number of parameters, J. R. Stat. Soc. B 71 (3), 671–683.}
#' @references \cite{Wang, L., Y. Kim, and R. Li (2013). Calibrating nonconvex penalized regression in ultra-high dimension, Ann. Stat. 41 (5), 2505–2536.}
#'
#' @examples
#' #Replicate the illustrative example from Bauer (202Xa).
#'
#'
#' p <- 300 #Number of variables. In Bauer (202Xa), p = 3000
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
#' colnames(X) <- seq_len(p)
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' single.bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' @importFrom irlba ssvd
#'
#' @export
bdsvd.ht <- function(X,
                     dof.lim,
                     standardize = TRUE,
                     anp = "2",
                     max.iter
) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA)")
    }

  n <- nrow(X)
  p <- ncol(X)

  if (missing(dof.lim)) {
    dof.lim <- c(0, p - 1)
    }
  dof.lim <- sort(round(dof.lim))
  if (dof.lim[1] < 0) {
    dof.lim[1] <- 0
    }
  if (dof.lim[2] > p - 1) {
    dof.lim[2] <- p - 1
    }
  dof.grid <- dof.lim[1]:dof.lim[2]

  X <- scale(X, center = TRUE, scale = standardize)

  ANP <- c("1", "2", "3")
  if (!(anp %in% ANP))
    stop(paste(anp), " is an invalid option for anp")

  if (anp == "2") {
    a_np <- function(n, p) {
      1 / 2 * log(n * p)
      }
  } else if (anp == "1") {
    a_np <- function(n, p) {
      1
      }
  } else if (anp == "3") {
    a_np <- function(n, p) {
      log(log(n * p))
      }
  }

  if (missing(max.iter)) {
    max.iter <- 500
    }

  BIC <- vector(length = length(dof.grid))
  i <- 1
  for (dof in dof.grid) {
    v <- tryCatch(suppressWarnings(irlba::ssvd(x = X, k = 1, n = p - dof, maxit = max.iter)$v), error = function(e) e)
    if (inherits(v, "error")) {
      v <- matrix(0, nrow = p, ncol = 1)
    }

    u <- X %*% v
    BIC[i] <- log(norm(X - u %*% t(v), type = "F")^2 / n / p) + sum(v != 0) * log(n * p) / n / p *  a_np(n, p)
    i <- i + 1
  }

  dof <- dof.grid[order(BIC)[1]]
  BIC <- cbind.data.frame(BIC)
  rownames(BIC) <- dof.grid

  return(list(dof = dof, BIC = BIC))
}



#' @title Data Matrix Structure According to the Detected Block Structure.
#'
#' @description Either sorts the data matrix \eqn{X} according to the detected block structure \eqn{X_1 , ... , X_b}, ordered by the number
#' of variables that the blocks contain. Or returns the detected submatrices each individually in a list object.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param block.structure Output of \code{bdsvd()} or \code{single.bdsvd()} which identified the block structure.
#'
#' @param output Should the output be the data matrix ordered according to the blocks (\code{"matrix"}), or
#' a list containing the submatrices (\code{"submatrices"}). Default is \code{"matrix"}.
#'
#' @param block.order A vector that contains the order of the blocks detected by \code{bdsvd()} or \code{single.bdsvd()}.
#' The vector must contain the index of each blocks exactly once. Default is \code{1:b} where \code{b} is the total number of blocks.
#'
#' @return
#' Either the data matrix \code{X} with columns sorted according to the detected blocks, or a list containing the detected
#' submatrices.
#'
#' @seealso \link{bdsvd}, \link{single.bdsvd}
#'
#' @references \cite{Bauer, J.O. (202Xa). High-dimensional block diagonal covariance structure detection using singular vectors.}
#'
#' @examples
#' #Toying with the illustrative example from Bauer (202Xa).
#'
#'
#' p <- 150 #Number of variables. In Bauer (202Xa), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
#' colnames(X) <- seq_len(p)
#'
#' #Compute iterative BD-SVD
#' bdsvd.obj <- bdsvd(X, standardize = FALSE)
#'
#' #Obtain the data matrix X, sorted by the detected blocks
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "matrix") )
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "matrix", block.order = c(2,1,3)) )
#'
#' #Obtain the detected submatrices X_1, X_2, and X_3
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "submatrices")[[1]] )
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "submatrices")[[2]] )
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "submatrices")[[3]] )
#'
#' @export
bdsvd.structure <- function(X,
                            block.structure,
                            output = "matrix",
                            block.order
) {

  if (!inherits(block.structure, "bdsvd"))
    stop("block.structure must be the outcome of bdsvd() or single.bdsvd().")

  OUTPUT <- c("matrix", "submatrices")
  if (!(output %in% OUTPUT))
    stop(paste(output), " is an invalid argument for output")

  p <- ncol(X)
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(1:p)
    }

  b <- length(block.structure)
  if (b == p)
    return(result)

  ifelse(missing(block.order),
         block.order <- seq_along(block.structure),
         block.order <- as.integer(block.order))
  if (!identical(sort(block.order), 1:b)) {
    stop("block.order must contain the index of each block exactly once.")
    }

  if (output == "matrix") {
    result <- X[, colnames(X) %in% block.structure[[block.order[1]]], drop = FALSE]
    for (i in block.order[-1]) {
      result <- cbind.data.frame(result,  X[, colnames(X) %in% block.structure[[i]], drop = FALSE])
    }
  } else {
    result <- list()
    for (i in block.order) {
      result <- c(result, list(X[, colnames(X) %in% block.structure[[block.order[i]]], drop = FALSE]))
    }
  }

  return(result)
}



#' @title Block Detection
#'
#' @description This function returns the block structure of a matrix.
#'
#' @param V Numeric matrix which either contains the loadings or is a covariance matrix.
#'
#' @param threshold All absolute values of \code{V} below the threshold are set to zero.
#'
#' @return
#' An object of class \code{Block} containing the features and columns indices corresponding to each detected block.
#'
#' @seealso \link{bdsvd}, \link{single.bdsvd}
#'
#' @references \cite{Bauer, J.O. (202Xa). High-dimensional block diagonal covariance structure detection using singular vectors.}
#'
#' @examples
#' #In the first example, we replicate the simulation study for the ad hoc procedure
#' #Est_0.1 from Bauer (202Xa). In the second example, we manually compute the first step
#' #of BD-SVD, which can be done using the bdsvd() and/or single.bdsvd(), for constructed
#' #sparse loadings
#'
#' #Example 1: Replicate the simulation study (a) from Bauer (202Xa) for the ad hoc
#' #procedure Est_0.1.
#'
#' p <- 500 #Number of variables
#' n <- 125 #Number of observations
#' b <- 500 #Number of blocks
#' design <- "a"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' #Perform the ad hoc procedure
#' detect.blocks(cvCovEst::scadEst(dat = X, lambda = 0.2), threshold = 0)
#'
#' #Example 2: Manually compute the first step of BD-SVD
#' #for some loadings V that mirror the two blocks
#' #("A", "B") and c("C", "D").
#'
#' V <- matrix(c(1,0,
#'               1,0,
#'               0,1,
#'               0,1), 4, 2, byrow = TRUE)
#'
#' rownames(V) <- c("A", "B", "C", "D")
#' detected.blocks <- detect.blocks(V)
#'
#' #Variables in block one with corresponding column index:
#' detected.blocks[[1]]@features
#' detected.blocks[[1]]@block.columns
#'
#' #Variables in block two with corresponding column index:
#' detected.blocks[[2]]@features
#' detected.blocks[[2]]@block.columns
#'
#' @export
detect.blocks <- function(V,
                          threshold = 0
                          ) {

  if (missing(V)) {
    stop("V is required.")
    }

  if (length(rownames(V)) == 0) {
    rownames(V) == as.character(seq_len(nrow(V)))
    }

  threshold.matrx <- get.threshold.matrix(V, threshold)

  return(get.blocks(threshold.matrx, rownames(V)))
}



#' @title Single Iteration of Block Detection Using Singular Vectors (BD-SVD).
#'
#' @description Performs a single iteration of BD-SVD: splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (which is a sparse
#' approximation of the first right singular vector, i.e., a vector with many zero values) that mirrors the
#' shape of the covariance matrix.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param dof Number of non-zero components in the sparse loading (degrees of freedom). If
#' \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom.
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package.
#'
#' @return
#' A list containing the feature names of the submatrices of \code{X}. It is either of length one (no
#' split) or length two (split into two submatrices).
#'
#' @seealso \link{bdsvd}, \link{bdsvd.ht}
#'
#' @references \cite{Bauer, J.O. (202Xa). High-dimensional block diagonal covariance structure detection using singular vectors.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #Replicate the illustrative example from Bauer (202Xa).
#'
#' \dontrun{
#'
#' p <- 300 #Number of variables. In Bauer (202Xa), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' single.bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' }
#'
#' @importFrom irlba ssvd
#'
#' @export
single.bdsvd <- function(X,
                  dof,
                  standardize = TRUE,
                  max.iter
                  ) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA)")
    }

  p <- ncol(X)

  if (missing(dof)) {
    stop("Enter the degrees of freedom")
    }

  if (missing(max.iter)) {
    max.iter <- 500
    }

  X <- scale(X, center = TRUE, scale = standardize)
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(seq_len(p))
  }
  feature.names <- colnames(X)
  eigen <- list()

  v <- tryCatch(suppressWarnings(irlba::ssvd(x = X, k = 1, n = p - dof, maxit = max.iter)$v), error = function(e) e)
  if (inherits(v, "error")) {
    stop("dof = ", dof, " do not fit the structure of the singular vectors. You can use bdsvd.ht to find a suitable value for dof.")
  }

  if (length(which(v == 0)) == 0) {
    return(list(feature.names))
    }

  eigen$vectors <- cbind(v, 0)
  eigen$vectors[which(v == 0), 2] <- 1
  rownames(eigen$vectors) <- feature.names


  detected.blocks <- detect.blocks(eigen$vectors, 0)
  result <- list()
  result[[1]] <- detected.blocks[[1]]@features
  result[[2]] <- detected.blocks[[2]]@features
  class(result) <- "bdsvd"

  return(result)
}
