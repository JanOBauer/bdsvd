calc.reliability <- function(X, R, p, reliability) {
  if (reliability == "spectral") {
    reliability <- 1 - norm(R, type = "2") / p
    return(reliability)
  }
}



calc.distance <- function(X, R, idx.cluster, feature.names, p, linkage) {

  if (linkage == "average") {
    return(1 - mean(abs(R[feature.names %in% idx.cluster[[1]], feature.names %in% idx.cluster[[2]]])))
  }

  if (linkage == "single") {
    return(1 - max(abs(R[feature.names %in% idx.cluster[[1]], feature.names %in% idx.cluster[[2]]])))
  }

  if (linkage == "RV") {
    R1  <- R[feature.names %in% idx.cluster[[1]], feature.names %in% idx.cluster[[1]]]
    R2  <- R[feature.names %in% idx.cluster[[2]], feature.names %in% idx.cluster[[2]]]
    R12 <- R[feature.names %in% idx.cluster[[1]], feature.names %in% idx.cluster[[2]]]
    return(1 - RV.coef(R1, R2, R12))
  }


}




hcsvd.ht <- function(X,
                     k,
                     linkage,
                     reliability,
                     R,
                     feature.names,
                     max.iter,
                     trace = TRUE
                     ) {

  p <- ncol(X)

  if (k == "Kaiser") {
    k <- sum(eigen(R)$values >= 1)
  } else {
    k <- p - 1
  }

  if (p == 2) {
    cluster1 <- colnames(X)[1]
    cluster2 <- colnames(X)[2]

    if (reliability == "linkage") {
      idx.cluster <- list(feature.names[1], feature.names[2])
      cluster.reliability <- calc.distance(X, R, idx.cluster, feature.names, p, linkage)
    } else {
      cluster.reliability <- calc.reliability(X, R, p, reliability)
    }

    return(list(cluster1 = cluster1, cluster2 = cluster2, cluster.reliability = cluster.reliability, k.p = NA))
  }


  dof.grid <- 1 : (p - 1)
  distance <- -Inf
  i <- 1

  for (dof in dof.grid) {
    V <- tryCatch(suppressWarnings(irlba::ssvd(x = X, k = k, n = p - dof, maxit = max.iter)$v), error = function(e) e)
    if (inherits(V, "error")) next

    for (i in 1:k) {
      v <- V[, i]
      idx.cluster <- list(feature.names[v != 0], feature.names[v == 0])
      distance.ht <- calc.distance(X, R, idx.cluster, feature.names, p, linkage)

      if (distance.ht > distance) {
        distance <- distance.ht
        cluster1 <- idx.cluster[[1]]
        cluster2 <- idx.cluster[[2]]
      }
    }
  }

  if (reliability == "linkage") {
    cluster.reliability <- distance
  } else {
    cluster.reliability <- calc.reliability(X, R, p, reliability)
  }

  return(list(cluster1 = cluster1, cluster2 = cluster2, cluster.reliability = cluster.reliability, k.p = k / p))
}





RV.coef <- function(R1, R2, R12) {
  RV <- sum(diag(crossprod(R12, R12))) / sqrt(sum(diag(crossprod(R1))) * sum(diag(crossprod(R2))))
  return(RV)
}





#' @title Correlation Matrix Simulation for HC-SVD
#'
#' @description This function generates correlation matrices based on the simulation studies described in Bauer (202Xb).
#'
#' @param p Number of variables.
#'
#' @param b Number of blocks.
#'
#' @param design Simulation design "a" or "b".
#'
#' @return
#' A correlation matrix according to the chosen simulation design.
#'
#' @references \cite{Bauer, J.O. (202Xb). Hierarchical variable clustering using singular vectors.}
#'
#' @examples
#' #The correlation matrix for simulation design (a) is given by
#' #R <- hcsvd.cov.sim(p = 100, b = 5, design = "a")
#'
#' @export
hcsvd.cor.sim <- function(p = p,
                          b = b,
                          design = design
) {

  DESIGN <- c("a", "b")
  if (!(design %in% DESIGN))
    stop(paste(design), " is an invalid design")

  if (design == "a") {
    d <- p / b / 4

    if (!((p / b) %% 4 == 0))
      stop("p/b must be divisible by 4 for simulation design a.")

    R <- matrix(0, p, p)
    for (block in 1:b) {
      eps <- runif(1, -0.05, 0.05)
      Rii <- matrix(0.25 + eps, 4 * d, 4 * d)
      for (j in 3:2) {
        eps <- runif(1, -0.05, 0.05)
        Rii[1 : (j * d), 1 : (j * d)] <- Rii[1 : (j * d), 1 : (j * d)] + (0.25 + eps) * rep(1, (j * d)) %*% t(rep(1, (j * d)))
      }
      Rii[which(kronecker(diag(4), matrix(1, d, d))  == 1)] <- 0.95
      R[((block - 1) * 4 * d + 1) : (block * 4 * d), ((block - 1) * 4 * d + 1) : (block * 4 * d)] <- Rii

    }
    diag(R) <- 1
    return(R)
  }

  if (design == "b") {
    d <- 3

    if (!((p / b) %% 3 == 0))
      stop("p/b must be divisible by 3 for simulation design b.")

    R <- matrix(0, d * b, d * b)
    Rii <- matrix(0, d, d)
    for (i in 1 : b) {
      for (j in 1 : d) {
        omega <- runif(1, 0.75, 0.85)
        Rii[, j] <- (-1)^(j) * omega^((j - 1)^2)
      }
      R[(1 + d * (i - 1)) : (d + d * (i - 1)), (1 + d * (i - 1)) : (d + d * (i - 1))] <- Rii
    }
    R[lower.tri(R, diag = TRUE)] <- t(R)[lower.tri(R, diag = TRUE)]
    diag(R) <- 1
    R
    R <- kronecker(R, matrix(1, p / b / 3, p / b / 3))
    return(R)
  }

}





#' @title Hierarchical Variable Clustering Using Singular Vectors (HC-SVD).
#'
#' @description Performs HC-SVD to reveal the hierarchical variable structure as descried in Bauer (202Xb). For this divise approach, each cluster is split into two clusters iteratively. Potential splits
#' are identified by the first sparse loadings (which are sparse approximations of the first right singular vectors, i.e., vectors with many zero values) that
#' mirror the masked shape of the correlation matrix. This procedure is continued until each variable lies in a single cluster.
#'
#' @param X Data matrix of dimension \eqn{n x p}. The data matrix is standardized during the analysis by \code{hcsvd}.
#'
#' @param k Number of sparse loadings to be used. This should be \code{"all"} for all sparse loadings, or \code{"Kaiser"} for as many sparse loadings as
#' there are eigenvalues larger or equal to one (see Bauer (202Xb) for details). Selecting \code{"Kaiser"} reduces computation time.
#'
#' @param linkage The linkage function to be used. This should be one of \code{"average"}, \code{"single"}, or
#' \code{"RV"} (for RV-coefficient).
#'
#' @param reliability By default, the value of each cluster equals the distance calculated by the chosen linkage function.
#' If preferred, the value of each cluster can be assigned by its reliability. When \code{reliability = spectral}, the reliability is
#' calculated by the averaged spectral norm.
#'
#' @param R Sample correlation matrix of \code{X}. By default, \code{R <- cov(X)}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loadings.
#' Default is \code{200}.
#'
#' @param trace Print out progress as \eqn{p-1} iterations for divisive hierarchical clustering are performed.
#' Default is \code{TRUE}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package.
#'
#' @return
#' A list with two components:
#' \item{dist.matrix}{
#'  The ultrametric distance matrix (cophenetic matrix) of the HC-SVD structure as an object of class \code{dist}.
#' }
#' \item{u.cor}{
#'  The ultrametric correlation matrix of \eqn{X} obtained by HC-SVD as an object of class \code{matrix}.
#' }
#' \item{k.p}{
#'  A vector of length \eqn{p-1} containing the ratio \eqn{k_i/p_i} of the \eqn{k_i} sparse loadings used relative to all sparse
#'  loadings \eqn{p_i} for the split of each cluster. The ratio is set to \code{NA} if the cluster contains only two variables as the search
#'  for sparse loadings that reflect the split is not required in this case.
#' }
#'
#' @references \cite{Bauer, J.O. (202Xb). Hierarchical variable clustering using singular vectors.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #We replicate the simulation study in Bauer (202Xb)
#'
#' \dontrun{
#' p <- 100
#' n <- 300
#' b <- 5
#' design <- "a"
#'
#' Rho <- hcsvd.cor.sim(p = p, b = b, design = "a")
#' X <- scale(mvtnorm::rmvnorm(300, mean=rep(0,100), sigma=Rho, checkSymmetry = FALSE))
#' colnames(X) = 1:ncol(X)
#' hcsvd.obj <- hcsvd(X, k = "Kaiser")
#'
#' #The dendrogram can be obtained from the ultrametric distance matrix:
#' plot(hclust(hcsvd.obj$dist.matrix))
#' }
#'
#'
#' @importFrom irlba ssvd
#' @importFrom stats cov
#' @importFrom stats as.dist
#'
#' @export
hcsvd <- function(X,
                  k = "all",
                  linkage = "single",
                  reliability,
                  R,
                  max.iter,
                  trace = TRUE
                  ) {

  if (anyNA(X)) {
    stop("X contains missing value indicator (NA)")
  }

  K <- c("all", "Kaiser")
  if (!(k %in% K))
    stop(paste(k), " is an invalid argument for k")

  LINKAGE <- c("average", "single", "RV")
  if (!(linkage %in% LINKAGE))
    stop(paste(linkage), " is an invalid linkage function")

  if (missing(reliability)) {
    reliability <- "linkage"
  } else {
    RELIABILITY <- c("spectral")
    if (!(reliability %in% RELIABILITY))
      stop(paste(reliability), " is an invalid argument for internal consistency reliability")
    }

  if (missing(max.iter)) {
    max.iter <- 500
    }

  p <- ncol(X)
  X <- scale(X)
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(seq_len(p))
    }

  k.p <- c()

  if (missing(R)) {
    R <- stats::cov(X)
    }

  dist.matrix <- matrix(0, p, p)

  feature.names <- colnames(X)
  dimnames(R) <- list(feature.names, feature.names)
  dimnames(dist.matrix) <- list(feature.names, feature.names)

  sub.matrices <- list(colnames(X))
  iter <- 0
  while (TRUE) {
    while (length(sub.matrices[[1]]) == 1) {
      sub.matrices <- sub.matrices[-1]

      if (length(sub.matrices) == 0) {
        u.cor <- 1 - dist.matrix
        dist.matrix <- stats::as.dist(dist.matrix)
        attr(dist.matrix, "Size") <- p
        cat("\r======== FINISHED ========                    ")
        cat("\n")
        return(list(dist.matrix = dist.matrix, u.cor = u.cor, k.p = k.p))
      }
    }

    current.features <- feature.names[feature.names %in% sub.matrices[[1]]]
    hcsvd.ht <- hcsvd.ht(X = X[, feature.names %in% sub.matrices[[1]]],
                         k = k,
                         linkage = linkage,
                         reliability = reliability,
                         R = R[feature.names %in% sub.matrices[[1]], feature.names %in% sub.matrices[[1]]],
                         feature.names = current.features,
                         max.iter = max.iter,
                         trace = trace)

    dist.matrix[feature.names %in% hcsvd.ht$cluster1, feature.names %in% hcsvd.ht$cluster2] <- hcsvd.ht$cluster.reliability
    dist.matrix[feature.names %in% hcsvd.ht$cluster2, feature.names %in% hcsvd.ht$cluster1] <- hcsvd.ht$cluster.reliability

    sub.matrices <- sub.matrices[-1]
    sub.matrices <- c(sub.matrices, list(hcsvd.ht$cluster1, hcsvd.ht$cluster2))

    cat("\riteration", iter <- iter + 1, "of", p - 1, "(", round(iter / (p - 1) * 100, 2), "%)                    ")

    k.p <- c(k.p, hcsvd.ht$k.p)
  }
}
