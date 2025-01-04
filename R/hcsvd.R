calc.distance <- function(R, idx.cluster, feature.names, p, linkage) {

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




#' @importFrom utils head
calc.sparse.loadings <- function (X, SVD, k, n, maxit = 500, tol = 0.001) {
  n <- ncol(X) - n
  if (length(n) != k) {
    n <- rep(n, length.out = k)
  }

  soft <- function(X, u, p) {
    y <- crossprod(X, u)
    a <- abs(y)
    z <- apply(a, 2, sort)
    lambda <- vapply(seq(length(p)), function(j)  z[p[j], j], pi, USE.NAMES = FALSE)
    sign(y) * pmax(sweep(a, 2, lambda, `-`), 0)
  }
  SVD$v <- SVD$d * SVD$v
  iter <- 0
  delta.u <- Inf
  while (delta.u > tol && iter < maxit) {
    u <- SVD$u
    SVD$v <- soft(X, SVD$u, n)
    xsv <- X %*% SVD$v

    SVD$u <- qr.Q(qr(xsv))
    SVD$u <- sweep(SVD$u, 2, apply(xsv, 2, function(x) sign(head(x[x != 0], 1)))/
                     apply(SVD$u, 2, function(x) sign(head(x[x != 0], 1))), `*`)


    delta.u <- max(1 - diag(abs(crossprod(u, SVD$u))))
    iter <- iter + 1
  }
  if (iter >= maxit)
    warning("Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.")
  SVD$v <- SVD$v %*% diag(1/sqrt(apply(SVD$v, 2, crossprod)), k,
                          k)
  list(v = SVD$v, u = SVD$u, iter = iter)
}




hcsvd.ht <- function(R,
                     k,
                     linkage,
                     labels,
                     max.iter,
                     trace = TRUE
) {

  p <- ncol(R)
  if (p == 2) {
    idx.cluster <- list(labels[1], labels[2])
    cluster.distance <- calc.distance(R, idx.cluster, labels, p, linkage)

    return(list(cluster1 = idx.cluster[1], cluster2 = idx.cluster[2], cluster.distance = cluster.distance, k.p = NA))
  }

  lambda <- eigen(R)$values
  if (is.numeric(k)) {
    if (k == 1){
      k <- p - 1
    } else {
      k <- min(round(k * p), p - 1)
    }
  }
  if (k == "Kaiser") {
    k <- min(sum(lambda >= 1) + 2, p - 1)
  }


  # Check if all elements are identical
  all.identical <- length(unique(R[upper.tri(R, diag = FALSE)])) == 1
  if(all.identical){
    print(paste0("all identical: rownames(R) = ", rownames(R)))
    idx.cluster <- list(labels[1], labels[-1])

    cluster.distance <- calc.distance(R, idx.cluster, labels, p, linkage)

    return(list(cluster1 = idx.cluster[[1]], cluster2 = idx.cluster[[2]], cluster.distance = cluster.distance, k.p = k / p))
  }


  dof.grid <- 1 : (p - 1)
  distance <- -Inf
  i <- 1

  for (dof in dof.grid) {
    V <- tryCatch(suppressWarnings(calc.sparse.loadings(X = R, SVD = irlba::irlba(R, k), k = k, n = p - dof, maxit = max.iter)$v), error = function(e) e)
    if (inherits(V, "error")) next

    for (i in seq_len(k)) {
      v <- V[, i]
      idx.cluster <- list(labels[v != 0], labels[v == 0])
      distance.ht <- calc.distance(R, idx.cluster, labels, p, linkage)

      if (distance.ht > distance) {
        distance <- distance.ht
        cluster1 <- idx.cluster[[1]]
        cluster2 <- idx.cluster[[2]]

      }
    }



  }

  cluster.distance <- distance

  return(list(cluster1 = cluster1, cluster2 = cluster2, cluster.distance = cluster.distance, k.p = k / p))
}




is.corr.matrix <- function(R) {

  if (!is.matrix(R)) {
    return(FALSE)
  }

  if(max(abs(R)) > 1) {
    return(FALSE)
  }

  if (nrow(R) != ncol(R)) {
    return(FALSE)
  }

  if (!isSymmetric(R)) {
    return(FALSE)
  }

  if (!all(diag(R) == 1)) {
    return(FALSE)
  }

  eigenvalues <- eigen(R, only.values = TRUE)$values
  if(any(eigenvalues < 0)) {
    return(FALSE)
  }

  return(TRUE)
}




RV.coef <- function(R1, R2, R12) {
  RV <- sum(diag(crossprod(R12, R12))) / sqrt(sum(diag(crossprod(R1))) * sum(diag(crossprod(R2))))
  return(RV)
}





#' @title Correlation Matrix Simulation for HC-SVD
#'
#' @description This function generates correlation matrices based on the simulation studies described in Bauer (202X).
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
#' @references \cite{Bauer, J.O. (202X). Divisive hierarchical clustering identified by singular vectors.}
#'
#' @examples
#' #The correlation matrix for simulation design (a) is given by
#' #R <- hcsvd.cov.sim(p = 40, b = 5, design = "a")
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
        eps <- runif(1, -0.1, 0.1)
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
#' @description Performs HC-SVD to reveal the hierarchical variable structure as descried in Bauer (202X). For this divise approach, each cluster is split into two clusters iteratively. Potential splits
#' are identified by the first sparse loadings (which are sparse approximations of the first right eigenvectors, i.e., vectors with many zero values, of the correlation matrix) that
#' mirror the masked shape of the correlation matrix. This procedure is continued until each variable lies in a single cluster.
#'
#' @param R A correlation matrix of dimension \eqn{p x p} or a data matrix of dimension \eqn{n x p} an be provided. If a data matrix is supplied, it must be indicated by setting
#' \code{is.corr = FALSE}, and the correlation matrix will then be calculated as \code{cor(X)}.
#'
#' @param k Number of sparse loadings to be used. This should be either a numeric value between zero and one to indicate percentages, or \code{"Kaiser"} for as many sparse loadings as
#' there are eigenvalues larger or equal to one. For a numerical value between zero and one, the number of sparse loadings is determined as the corresponding share of the total number of loadings.
#' E.g., \code{k = 1} (100%) will use all sparse loadings (see Bauer (202X) for details).
#'
#' @param linkage The linkage function to be used. This should be one of \code{"average"}, \code{"single"}, or
#' \code{"RV"} (for RV-coefficient).
#'
#' @param is.corr Is the supplied object a correlation matrix. Default is \code{TRUE} and this parameter must be set to \code{FALSE} is
#' a data matrix instead of a correlation matrix is supplied.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loadings.
#' Default is \code{200}.
#'
#' @param trace Print out progress as \eqn{p-1} iterations for divisive hierarchical clustering are performed.
#' Default is \code{TRUE}.
#'
#' @details
#' The sparse loadings are computed using the method of Shen and Huang (2008), which is implemented based on the code
#' of Baglama, Reichel, and Lewis in \link[irlba]{ssvd} \{\link[=irlba]{irlba}\}, with slight modifications to suit our method.
#'
#' @return
#' A list with four components:
#' \item{hclust}{
#'  The clustering structure identified by HC-SVD as an object of type \code{hclust}.
#' }
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
#' @references \cite{Bauer, J.O. (202X). Divisive hierarchical clustering identified by singular vectors.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #We replicate the simulation study in Bauer (202X)
#'
#' \dontrun{
#' p <- 20
#' n <- 500
#' b <- 5
#' design <- "a"
#'
#' Rho <- hcsvd.cor.sim(p = p, b = b, design = "a")
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma = Rho, checkSymmetry = FALSE)
#' R <- cor(X)
#' hcsvd.obj <- hcsvd(R)
#'
#' #The object of hclust with corresponding dendrogram can be obtained
#' #directly from hcsvd.obj$hclust:
#' hc <- hcsvd.obj$hclust
#' plot(hc)
#'
#' #The dendrogram can also be obtained from the ultrametric distance matrix:
#' plot(hclust(hcsvd.obj$dist.matrix))
#' }
#'
#'
#' @importFrom stats cov
#' @importFrom stats as.dist
#'
#' @export
hcsvd <- function(R,
                  k = "Kaiser",
                  linkage = "average",
                  is.corr = TRUE,
                  max.iter,
                  trace = TRUE
) {

  if (is.corr && !is.corr.matrix(R)) {
    stop("R must be a correlation matrix. Set 'is.corr = FALSE' if you want to supply a data matrix")
  }

  if (!is.corr) {
    X <- R
    if (anyNA(X)) {
      stop("X contains missing value indicator (NA)")
    }
    R <- stats::cor(X)
  }

  if (!((k == "Kaiser") | (is.numeric(k) & (k > 0) & (k <= 1)))) {
    stop(paste(k), " is an invalid argument for k")
  }

  LINKAGE <- c("average", "single", "RV")
  if (!(linkage %in% LINKAGE)) {
    stop(paste(linkage), " is an invalid linkage function")
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  p <- ncol(R)
  if (length(colnames(R)) == 0 | length(rownames(R)) == 0) {
    dimnames(R) <- list(as.character(seq_len(p)), as.character(seq_len(p)))
  }

  if (length(unique(colnames(R))) != p) {
    stop("Variable names are not unique")
  }

  k.p <- c()

  dist.matrix <- matrix(0, p, p)
  labels <- colnames(R)
  dimnames(dist.matrix) <- list(labels, labels)

  merge <- matrix(0, p - 1, 2)
  height <- vector(length = p - 1)
  order <- labels

  sub.matrices <- list(colnames(R))
  cluster.count <- p - 2
  for (iter in 1:(p - 1)) {
    current.labels <- labels %in% sub.matrices[[1]]
    hcsvd.ht <- hcsvd.ht(R = R[current.labels, current.labels],
                         k = k,
                         linkage = linkage,
                         labels = labels[current.labels],
                         max.iter = max.iter,
                         trace = trace)

    k.p <- c(k.p, hcsvd.ht$k.p)

    cluster.rows <- labels %in% hcsvd.ht$cluster1
    cluster.cols <- labels %in% hcsvd.ht$cluster2
    dist.matrix[cluster.rows, cluster.cols] <- hcsvd.ht$cluster.distance
    dist.matrix[cluster.cols, cluster.rows] <- hcsvd.ht$cluster.distance

    sub.matrices <- sub.matrices[-1]

    height[p - iter] <- hcsvd.ht$cluster.distance
    for (i in 1:2) {
      cluster <- hcsvd.ht[[paste0("cluster", i)]]
      if (length(cluster) != 1) {
        merge[p - iter, i] <- cluster.count
        cluster.count <- cluster.count - 1
        sub.matrices <- append(sub.matrices, list(cluster))
      } else {
        merge[p - iter, i] <- -which(labels == cluster)
      }
    }

    order.cluster1 <- order %in% hcsvd.ht$cluster1
    order.cluster2 <- order %in% hcsvd.ht$cluster2
    order[order.cluster1 | order.cluster2] <- c(order[order.cluster2], order[order.cluster1])

    cat(sprintf("\rSplit %d out of %d (%.2f%%)           ", iter, p - 1, iter / (p - 1) * 100))
  }

  ordered.height <- order(height)

  merge <- merge[ordered.height, ]
  height <- height[ordered.height]

  not.changed <- matrix(TRUE, p - 1, 2)
  for (i in seq_len(p - 1)) {
    change.idx <- which(merge == i)
    merge[merge == i & not.changed] <- which(ordered.height == i)
    not.changed[change.idx] <- FALSE
  }

  hclust <- list(merge = merge, height = height, order = match(order, labels), labels = labels, method = linkage)
  class(hclust) <- "hclust"

  u.cor <- 1 - dist.matrix
  dist.matrix <- stats::as.dist(dist.matrix)
  attr(dist.matrix, "Size") <- p
  cat("\r======== FINISHED ========                    ")
  cat("\n")
  return(list(hclust = hclust, dist.matrix = dist.matrix, u.cor = u.cor, k.p = k.p))
}











