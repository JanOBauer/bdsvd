## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib bdsvd, .registration = TRUE
## usethis namespace: end


create.block <- function(feature.names, selected.features, block.columns) {
  if (length(feature.names) > 0) {
    selected.features <- feature.names[selected.features]
  }
  return(list(features = selected.features, block.columns = block.columns))
}

get.blocks <- function(threshold.matrix, feature.names) {
  p <- nrow(threshold.matrix)
  k <- ncol(threshold.matrix)
  columns <- seq_len(k)
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

  return(structure(blocks, class = "blocks"))
}

get.threshold.matrix <- function(loadings, threshold) {
  return((abs(loadings) > threshold) * 1)
}

find.blocks <- function(mat, column.idx) {
  if (length(column.idx) == 1) {
    rows <- which(mat[, column.idx] != 0)
  } else {
    rows <- which(rowSums(mat[, column.idx] != 0) > 0)
  }

  if (length(rows) == 1) {
    cols <- which(mat[rows, ] != 0)
  } else {
    cols <- which(colSums(mat[rows, ] != 0) > 0)
  }

  return(cols)
}



#' @title Block Detection Using Singular Vectors (BD-SVD).
#'
#' @description Performs BD-SVD iteratively to reveal the block structure. Splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (which is a sparse approximation of the
#' first right singular vector, i.e., a vector with many zero values) that mirrors the shape of the covariance matrix. This
#' procedure is continued iteratively until the block diagonal structure has been revealed.
#'
#' The data matrix ordered according to this revealed block diagonal structure can be obtained by \code{\link{bdsvd.structure}}.
#'
#' @param X Data matrix of dimension \eqn{n}x\eqn{p} with possibly \eqn{p >> n}.
#'
#' @param dof.lim Interval limits for the number of non-zero components in the sparse loading (degrees of freedom).
#' If \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom. Default is \code{dof.lim <- c(0, p-1)} which is highly recommended to check for
#' all levels of sparsity.
#'
#' @param anp Which regularization function should be used for the HBIC.
#' \itemize{
#'   \item \code{"1"}: implements \eqn{a_{np} = 1} which corresponds to the BIC.
#'   \item \code{"2"}: implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (2025).
#'   \item \code{"3"}: implements \eqn{a_{np} = log(log(np))}.
#'   \item \code{"4"}: implements \eqn{a_{np} = log(log(p))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#' }
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @param scores Compute scores?
#'
#' @param verbose Print out progress as iterations are performed. Default is \code{TRUE}.
#'
#' @details
#' The sparse loadings are computed using the method proposed by Shen & Huang (2008). The corresponding implementation is written in \code{Rcpp}/\code{RcppArmadillo}
#' for computational efficiency and is based on the \code{R} implementation by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} \pkg{irlba}.
#' However, the implementation has been adapted to better align with the scope of the \pkg{bdsvd} package.
#'
#' @return
#' A list containing the feature names of the submatrices of \code{X}. The length of the list equals
#' the number of submatrices.
#'
#' @seealso \code{\link{bdsvd.structure}}, \code{\link{bdsvd.ht}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#' @references \cite{Wang, H., B. Li, and C. Leng (2009). Shrinkage tuning parameter selection with a diverging number of parameters, J. R. Stat. Soc. B 71 (3), 671–683.}
#' @references \cite{Wang, L., Y. Kim, and R. Li (2013). Calibrating nonconvex penalized regression in ultra-high dimension, Ann. Stat. 41 (5), 2505–2536.}
#'
#' @examples
#' #Replicate the simulation study (c) from Bauer (2025).
#'
#' \dontrun{
#' p <- 500 #Number of variables
#' n <- 500 #Number of observations
#' b <- 10  #Number of blocks
#' design <- "c" #Simulation design "a", "b", "c", or "d".
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- seq_len(p)
#'
#' bdsvd(X, standardize = FALSE, anp = "4")
#' }
#'
#' @importFrom irlba irlba
#'
#' @export
bdsvd <- function(X,
                  dof.lim,
                  anp = "2",
                  standardize = TRUE,
                  max.iter,
                  scores = FALSE,
                  verbose = TRUE
) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA).")
  }

  X <- scale(X, scale = standardize)
  p <- ncol(X)

  if (is.null(colnames(X))) {
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

  ANP <- c("1", "2", "3", "4")
  if (!(anp %in% ANP))
    stop(sprintf("%s is an invalid option for anp.", anp))

  if (missing(max.iter)) {
    max.iter <- 500
  }

  sub.matrices <- list(colnames(X))
  results <- list()

  trace.blocks <- 0
  b <- 1
  while (TRUE) {
    if (length(sub.matrices) == 0) {
      break  #Stop when splitting is no longer supported
    }

    if (length(sub.matrices[[1]]) == 1) {
      if (verbose) {
        trace.blocks <- trace.blocks + 1
        trace.perc <- round(trace.blocks / p * 100, 2)
        cat("\rProgress:", trace.blocks, "variables (", trace.perc, "%) processed and split into", b, "blocks.    ")
      }
      b <- b + 1
      results <- c(results, sub.matrices[[1]])
      sub.matrices <- sub.matrices[-1]
      next
    }



    if(scores) {
      X.sub <- X[, colnames(X) %in% sub.matrices[[1]], drop = FALSE]
      v <- abs(irlba(X.sub, nv = 1)$v)
      main.signal <- colnames(X.sub)[which.max(v)]

      threshold <- 0.5 / sqrt(ncol(X.sub))
      X.sub.thresh <- X.sub[, v > threshold, drop = FALSE]


      if (ncol(X.sub.thresh) <= 1) {
        signal.block <- main.signal
      } else {
        dof.split <- bdsvd.ht(X = X.sub.thresh,
                              dof.lim = dof.lim,
                              standardize = FALSE,
                              anp = anp,
                              max.iter = max.iter)$dof
        sub.results <- single.bdsvd(X = X.sub.thresh, standardize = FALSE, dof = dof.split)
        signal.block <- Filter(function(x) main.signal %in% x, sub.results)[[1]]
      }

      other.block <- setdiff(colnames(X.sub), signal.block)

      if (length(other.block) == 0) {
        sub.results <- list(signal.block)
      } else {
        sub.results <- list(signal.block, other.block)
      }

    } else {
      X.sub <- X[, colnames(X) %in% sub.matrices[[1]], drop = FALSE]
      dof.split <- bdsvd.ht(X = X.sub,
                            dof.lim = dof.lim,
                            standardize = FALSE,
                            anp = anp,
                            max.iter = max.iter)$dof
      sub.results <- single.bdsvd(X = X.sub, standardize = FALSE, dof = dof.split)
    }



    if (length(sub.matrices) == 1) {
      sub.matrices <- list()
    } else {
      sub.matrices <- sub.matrices[-1]
    }

    if (dof.split == 0) { #previously: length(sub.results) == 1
      if (verbose) {
        trace.blocks <- trace.blocks + length(sub.results[[1]])
        trace.perc <- round(trace.blocks / p * 100, 2)
        cat("\rProgress:", trace.blocks, "variables (", trace.perc, "%) processed and split into", b, "blocks.    ")
      }
      b <- b + 1
      results <- c(results, sub.results)
    } else {
      sub.matrices <- c(sub.matrices, sub.results)
    }

  }

  cat("\n")
  out <- results[order(sapply(results, length), decreasing = TRUE)]
  return(structure(out, class = "bdsvd"))

}



#' @title Covariance Matrix Simulation for BD-SVD
#'
#' @description This function generates covariance matrices based on the simulation studies described in Bauer (2025).
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
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#'
#' @examples
#' #The covariance matrix for simulation design (a) is given by
#' Sigma <- bdsvd.cov.sim(p = 500, b = 500, design = "a")
#'
#' @importFrom stats runif
#'
#' @export
bdsvd.cov.sim <- function(p,
                          b,
                          design
) {

  DESIGN <- c("a", "b", "c", "d")
  if (!(design %in% DESIGN))
    stop(sprintf("%s is an invalid design.", design))

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
#' @param anp Which regularization function should be used for the HBIC.
#' \itemize{
#'   \item \code{"1"}: implements \eqn{a_{np} = 1} which corresponds to the BIC.
#'   \item \code{"2"}: implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (2025).
#'   \item \code{"3"}: implements \eqn{a_{np} = log(log(np))}.
#'   \item \code{"4"}: implements \eqn{a_{np} = log(log(p))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#' }
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @details
#' The sparse loadings are computed using the method proposed by Shen & Huang (2008). The corresponding implementation is written in \code{Rcpp}/\code{RcppArmadillo}
#' for computational efficiency and is based on the \code{R} implementation by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} \pkg{irlba}.
#' However, the implementation has been adapted to better align with the scope of the \pkg{bdsvd} package. The computation of the HBIC is outlined in Bauer (2025).
#'
#' @return
#' \item{dof}{
#'   The optimal number of nonzero components (degrees of freedom) according to the HBIC.
#' }
#' \item{BIC}{
#'   The HBIC for the different numbers of nonzero components.
#' }
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#' @references \cite{Wang, H., B. Li, and C. Leng (2009). Shrinkage tuning parameter selection with a diverging number of parameters, J. R. Stat. Soc. B 71 (3), 671–683.}
#' @references \cite{Wang, L., Y. Kim, and R. Li (2013). Calibrating nonconvex penalized regression in ultra-high dimension, Ann. Stat. 41 (5), 2505–2536.}
#'
#' @examples
#' #Replicate the illustrative example from Bauer (2025).
#'
#'
#' p <- 300 #Number of variables. In Bauer (2025), p = 3000
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- seq_len(p)
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' single.bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' @importFrom irlba irlba
#'
#' @export
bdsvd.ht <- function(X,
                     dof.lim,
                     standardize = TRUE,
                     anp = "2",
                     max.iter
) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA).")
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

  X <- scale(X, center = FALSE, scale = standardize)

  ANP <- c("1", "2", "3", "4")
  if (!(anp %in% ANP))
    stop(sprintf("%s is an invalid option for anp.", anp))

  if (missing(max.iter)) {
    max.iter <- 500
  }

  suppressWarnings(SVD <- irlba(X, nv = 1))
  BIC <- calc_BIC(X, SVD$u, SVD$v, dof.grid, max.iter, as.integer(anp), n, p)

  dof <- dof.grid[which.min(BIC)]
  BIC <- data.frame(BIC = BIC, row.names = dof.grid)

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
#' @seealso \code{\link{bdsvd}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#'
#' @examples
#' #Toying with the illustrative example from Bauer (2025).
#'
#'
#' p <- 150 #Number of variables. In Bauer (2025), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
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

  if (!inherits(block.structure, c("bdsvd", "blocks")))
    stop(sprintf(
      "block.structure must be a 'bdsvd' or 'blocks' object (got: '%s').\nE.g., pass the result of bdsvd(), single.bdsvd(), or detect.blocks().",
      paste(class(block.structure), collapse = "/")
    ), call. = FALSE)

  if (inherits(block.structure, "blocks")) {
    block.structure <- lapply(block.structure, `[[`, "features")
    class(block.structure) <- "bdsvd"
  }

  OUTPUT <- c("matrix", "submatrices")
  if (!(output %in% OUTPUT))
    stop(sprintf("%s is an invalid argument for output.", output))

  p <- ncol(X)
  if (is.null(colnames(X))) {
    colnames(X) <- as.character(1:p)
  }

  if (length(block.structure) == p) {
    warning("Number of blocks equals number of variables. X is returned unchanged.", call. = FALSE)
    return(X)
  }

  ifelse(missing(block.order),
         block.order <- seq_along(block.structure),
         block.order <- as.integer(block.order))

  if (output == "matrix") {
    result <- do.call(cbind.data.frame, lapply(block.order, function(i) {
      X[, colnames(X) %in% block.structure[[i]], drop = FALSE]
    }))
  } else {
    result <- lapply(block.order, function(i) X[, colnames(X) %in% block.structure[[i]], drop = FALSE])
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
#' An object of class \code{blocks} containing the features and columns indices corresponding to each detected block.
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#'
#' @examples
#' #In the first example, we replicate the simulation study for the ad hoc procedure
#' #Est_0.1 from Bauer (2025). In the second example, we manually compute the first step
#' #of BD-SVD, which can be done using the bdsvd() and/or single.bdsvd(), for constructed
#' #sparse loadings
#'
#' #Example 1: Replicate the simulation study (a) from Bauer (2025) for the ad hoc
#' #procedure Est_0.1.
#'
#'\dontrun{
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
#' }
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
#' detected.blocks[[1]]$features
#' detected.blocks[[1]]$block.columns
#'
#' #Variables in block two with corresponding column index:
#' detected.blocks[[2]]$features
#' detected.blocks[[2]]$block.columns
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
#' The sparse loadings are computed using the method proposed by Shen & Huang (2008). The corresponding implementation is written in \code{Rcpp}/\code{RcppArmadillo}
#' for computational efficiency and is based on the \code{R} implementation by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} \pkg{irlba}.
#' However, the implementation has been adapted to better align with the scope of the \pkg{bdsvd} package.
#'
#' @return
#' A list containing the feature names of the submatrices of \code{X}. It is either of length one (no
#' split) or length two (split into two submatrices).
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{bdsvd.ht}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #Replicate the illustrative example from Bauer (2025).
#'
#' \dontrun{
#'
#' p <- 300 #Number of variables. In Bauer (2025), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- 1:p
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' single.bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' }
#'
#' @importFrom irlba irlba
#'
#' @export
single.bdsvd <- function(X,
                         dof,
                         standardize = TRUE,
                         max.iter
) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA).")
  }

  p <- ncol(X)

  if (missing(dof)) {
    stop("Enter degrees of freedom (dof).")
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  X <- scale(X, center = TRUE, scale = standardize)
  if (is.null(colnames(X))) {
    colnames(X) <- as.character(seq_len(p))
  }

  if (length(unique(colnames(X))) != p) {
    stop("Variable names are not unique.")
  }

  feature.names <- colnames(X)
  eigen <- list()


  suppressWarnings(SVD <- irlba(X, nv = 1))
  v = calc_one_sparse_v_cpp(X, SVD$v, SVD$u, dof, max.iter)$v

  if (length(which(v == 0)) == 0) {
    return(list(feature.names))
  }

  eigen$vectors <- cbind(v, 0)
  eigen$vectors[which(v == 0), 2] <- 1
  rownames(eigen$vectors) <- feature.names


  detected.blocks <- detect.blocks(eigen$vectors, 0)
  result <- list()
  result[[1]] <- detected.blocks[[1]]$features
  result[[2]] <- detected.blocks[[2]]$features
  return(structure(result, class = "bdsvd"))
}














#' @export
print.blocks <- function(x, ...) {
  cat("Number of identified blocks:", length(x), "\n")
  invisible(x)
}



#' @export
print.bdsvd <- function(x, ...) {
  cat("Number of identified blocks:", length(x), "\n")
  cat("Access blocks with [[ ]], e.g., object[[1]] for the first block.\n")
  invisible(x)
}



