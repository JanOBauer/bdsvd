#' @importFrom irlba svdr
#' @importFrom stats var
#' @importFrom stats cov
calc.pca <- function(X, K, method, covariance, orthogonal) {

  n <- nrow(X)
  p <- ncol(X)

  if (p == 1) {
    v <- 1
    if (covariance) {
      l <- var(X)
    } else {
      l <- 1
    }

    return(list(l = l, v = v, K = 1))
  }

  if (method == "CDM") {
    n1 <- n %/% 2
    n2 <- n - n1

    if (K >= n1) {
      warning("K must be smaller than ", n1)
      K <- n1
    }

    if (K > p) {
      K <- p
    }

    X1 <- X[seq_len(n1), ]
    X2 <- X[seq(n1 + 1, n), ]
    S1 <- 1/sqrt(n1*n2) * X1 %*% t(X2)

    svd.S1 <- svdr(S1, k = K)

    lambda.tilde <- svd.S1$d
    u1.tilde <- svd.S1$u
    u2.tilde <- svd.S1$v

    h1.tilde <- sweep(t(X1) %*% u1.tilde, 2, sqrt(n * lambda.tilde), "/")
    h2.tilde <- sweep(t(X2) %*% u2.tilde, 2, sqrt(n * lambda.tilde), "/")

    signs <- sign(colSums(h1.tilde * h2.tilde))
    signs[signs == 0] <- 1
    h.tilde <- (h1.tilde + signs * h2.tilde) / 2
    h.tilde <- sweep(h.tilde, 2, sqrt(colSums(h.tilde^2)), "/")

    if (orthogonal) {
      svd.h.tilde <- svd(h.tilde)
      h.tilde <- svd.h.tilde$u %*% t(svd.h.tilde$v)
    }

    return(list(l = lambda.tilde, v = h.tilde, K = K))
  }

  if (method == "DM") {
    SD <- 1/n * X %*% t(X)

    svd.SD <- svdr(SD, k = K)

    lambda.tilde <- svd.SD$d
    u.tilde <- svd.SD$u

    h.tilde <- sweep(t(X) %*% u.tilde, 2, sqrt(n * lambda.tilde), "/")

    if (orthogonal) {
      svd.h.tilde <- svd(h.tilde)
      h.tilde <- svd.h.tilde$u %*% t(svd.h.tilde$v)
    }

    return(list(l = lambda.tilde, v = h.tilde, K = K))

  }

  if (method == "PCA") {

    if (n < p) {
      warning("The number of variables for a block exceeds the number of observations. It is recommended to use CDM or DM as the method")
    }

    S <- cov(X)
    svd.S <- eigen(S)
    l <- svd.S$values[1:K]
    v <- svd.S$vectors[, 1:K]

    return(list(l = l, v = v, K = K))
  }



}



#' @title High Dimensional Principal Component Analysis
#'
#' @description Performs a principal component analysis on the given data matrix using the methods of Yata and Aoshima (2009, 2010).
#'
#' @param X Data matrix of dimension \eqn{n}x\eqn{p} with possibly \eqn{p >> n}.
#'
#' @param K Number of principal components to be computed. If \code{K} is larger than the number of variables \eqn{p} contained in the data matrix, \eqn{K = p - 1} loadings are computed.
#'
#' @param scale Should the variables be scaled to have unit variance before the analysis takes place. Default is \code{TRUE}.
#'
#' @param method Which method should be used to calculate the eigenvectors (loadings) and eigenvalues. \code{method = "DM"} uses the method by Yata and Aoshima (2009) and
#' \code{method = "CDM"} uses the method by Yata and Aoshima (2010).
#'
#' @param orthogonal The estimated eigenvectors (loadings) computed using \code{method = "CDM"} (Yata and Aoshima, 2010) are orthogonal in the limit thus only approximately orthogonal
#' in the finite sample case. Should the loadings be orthogonalized. Default is \code{FALSE}.
#'
#' @details
#' This function performs principal component analysis using either the DM approach as described in Yata, K., Aoshima, M. (2009), or the CDM approach (Yata, K., Aoshima, M., 2010)
#' Note that there is also a code implementation of CDM available at 'Aoshima Lab' (\url{https://github.com/Aoshima-Lab/HDLSS-Tools/tree/main/CDM}) provided by Makoto Aoshima.
#'
#' @return
#' A list with the following components:
#' \item{v}{
#'  The first \code{K} estimated sparse singular vectors (loadings) if the data matrix \code{X}. The eigenvectors are
#'  orthogonalized if \code{orthogonal = TRUE}.
#' }
#' \item{l}{
#'  The corresponding first estimated eigenvalues of the identified block diagonal covariance matrix.
#' }
#' \item{K}{
#'  The number of sparse singular vectors (loadings) that have been computed.
#' }
#'
#' @references \cite{Yata, K., Aoshima, M. (2009). PCA consistency for non-Gaussian data in high dimension, low sample size contex, Commun. Stat. - Theory Methods 38, 2634–2652.}
#' @references \cite{Yata, K., Aoshima, M. (2010). Effective PCA for high-dimension, low-sample-size data with singular value decomposition of cross data matrix, J. Multivar. Anal. 101, 2060–2077.}
#'
#' @examples
#' #Example: run IS-PCA on a gene expression data set with two tissue types
#'
#' if (requireNamespace("dslabs", quietly = TRUE)) {
#' data("tissue_gene_expression", package = "dslabs")
#'
#' #We only select the two tissue types kidney (6) and liver (7)
#' Y <- as.numeric(tissue_gene_expression$y)
#' X <- scale(tissue_gene_expression$x[Y %in% c(6, 7), ], scale = FALSE)
#' Y <- Y[Y %in% c(6, 7)]
#'
#' # Run PCA
#' pca.obj <- cdm.pca(X, K = 2)
#' PC <- X %*% pca.obj$v
#'
#' # Plot the first two principal components
#' plot(PC, pch = Y-5, xlab = "PC1", ylab = "PC2")
#' }
#'
#' @export
cdm.pca <- function(X, K = 1, method = "CDM", scale = TRUE, orthogonal = FALSE) {

  n <- nrow(X)
  p <- ncol(X)

  if (p == 1) {
    v <- 1
    if (scale) {
      l <- 1
    } else {
      l <- var(X)
    }

    return(list(l = l, v = v, K = 1))
  }

  METHOD <- c("CDM", "DM")
  if (!(method %in% METHOD))
    stop(sprintf("%s is an invalid option for method.", method))

  X <- scale(X, scale = scale)

  if (K >= p) {
    warning("K must be smaller than p = ", p)
    K <- p - 1
  }

  pca.obj <- calc.pca(X = X, K = K, method = method, covariance = !scale, orthogonal = orthogonal)
  return(list(l = pca.obj$l, v = pca.obj$v, K = K))

}




#' @title Inherently Sparse Principal Component Analysis (IS-PCA).
#'
#' @description Performs IS-PCA
#'
#' @param X Data matrix of dimension \eqn{n}x\eqn{p} with possibly \eqn{p >> n}.
#'
#' @param K Number of singular vectors (loadings) to be computed for each identified submatrix. This means that \eqn{b*K} components are calculated if \eqn{b} blocks are identified by BD-SVD
#' (see Bauer (2025) and \code{\link{bdsvd}} for details). If \code{K} is larger than the number of variables \eqn{p_i} contained in a submatrix, \eqn{K = p_i} loadings are computed.
#'
#' @param block.structure Underlying block structure. This parameter is optional as otherwise the functions runs \code{bdsvd()}
#' to identify the underlying block structure. When supplied, it must be a '\code{bdsvd}', '\code{blocks}', or '\code{ispca}' object. E.g., pass the result of \code{bdsvd()}, \code{single.bdsvd()},
#' \code{ispca()}, or \code{detect.blocks()}. An identified block structure by any other method can be supplied using \code{detect.blocks()} (see example below).
#'
#' @param anp Which regularization function should be used for the HBIC.
#' \itemize{
#'   \item \code{"1"}: implements \eqn{a_{np} = 1} which corresponds to the BIC.
#'   \item \code{"2"}: implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (2025).
#'   \item \code{"3"}: implements \eqn{a_{np} = log(log(np))}.
#'   \item \code{"4"}: implements \eqn{a_{np} = log(log(p))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#' }
#'
#' @param covariance Perform IS-PCA on the covariance (\code{TRUE}) or correlation matrix (\code{FALSE}). Default is \code{TRUE}.
#'
#' @param method Which method should be used to calculate the eigenvectors (loadings) and eigenvalues. \code{method = "DM"} uses the method by Yata and Aoshima (2009),
#' \code{method = "CDM"} uses the method by Yata and Aoshima (2010), and \code{method = "PCA"} uses \code{\link[base]{svd}}. The methods by Yata and Aoshima (2009, 2010) are consistent in
#' HDLSS settings. Default is \code{method = "CDM"}.
#'
#' @param orthogonal The estimated eigenvectors (loadings) computed using \code{method = "CDM"} (Yata and Aoshima, 2010) are orthogonal in the limit thus only approximately orthogonal
#' in the finite sample case. Should the loadings be orthogonalized? Default is \code{FALSE}.
#'
#' @param standardize Standardize the data for block detection using BD-SVD. Default is \code{TRUE}.
#' Note: this does not affect the parameter \code{covariance}. PCA can be computed on the covariance
#' matrix regardless.
#'
#' @param verbose Print out progress for \code{\link{bdsvd}} as iterations are performed. Default is \code{TRUE}.
#'
#' @details
#' This function performs inherently sparse principal component analysis (IS-PCA) as introduced in Bauer (2026).
#'
#' @return
#' A list with the following components:
#' \describe{
#' \item{v}{
#'  The first estimated eigenvectors of the identified block diagonal covariance matrix (i.e., the loadings) as an
#'  object of type \code{matrix}. The eigenvectors are orthogonalized if \code{orthogonal = TRUE}.
#' }
#' \item{l}{
#'  The corresponding first estimated eigenvalues of the identified block diagonal covariance matrix.
#' }
#' \item{exp.var}{
#'  The explained variance of the first estimated eigenvalues \code{l}.
#' }
#' \item{X.b}{
#'  The detected submatrices using \code{\link{bdsvd}} as a \code{list} object.
#' }
#' \item{block.structure}{
#'  Either the block structure detected by \code{bdsvd()} or the user-supplied \code{block.structure}, depending on the input.
#' }
#' }
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{bdsvd.structure}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#' @references \cite{Bauer, J.O. (2026). Beyond regularization: inherently sparse principal component analysis. Stat. Comp.}
#' @references \cite{Yata, K., Aoshima, M. (2009). PCA consistency for non-Gaussian data in high dimension, low sample size contex, Commun. Stat. - Theory Methods 38, 2634–2652.}
#' @references \cite{Yata, K., Aoshima, M. (2010). Effective PCA for high-dimension, low-sample-size data with singular value decomposition of cross data matrix, J. Multivar. Anal. 101, 2060–2077.}
#'
#' @examples
#' #Example 1: run IS-PCA on a gene expression data set with two tissue types
#'
#' if (requireNamespace("dslabs", quietly = TRUE)) {
#' data("tissue_gene_expression", package = "dslabs")
#'
#' #We only select the two tissue types kidney (6) and liver (7)
#' Y <- as.numeric(tissue_gene_expression$y)
#' X <- scale(tissue_gene_expression$x[Y %in% c(6, 7), ], scale = FALSE)
#' Y <- Y[Y %in% c(6, 7)]
#'
#' # Run IS-PCA
#' ispca.obj <- ispca(X = X, anp = "1")
#' vhat <- ispca.obj$v[, 1:2]
#' ispc <- X %*% vhat
#'
#' # Percentage of non-zero components in the first two loadings
#' round(colSums(vhat != 0)/ncol(X), 2)
#'
#' # Plot the first two principal components
#' plot(ispc, pch = Y-5, xlab = "PC1", ylab = "PC2", main = "IS-PCA")
#'
#' # Compare to CDM-PCA (see cdm.pca(...))
#' pca.obj <- cdm.pca(X = X, K = 2)
#' pc <- X %*% pca.obj$v
#'
#' par(mfrow = c(1, 2))
#' plot(ispc, pch = Y-5, xlab = "PC1", ylab = "PC2", main = "IS-PCA")
#' plot(pc, pch = Y-5, xlab = "PC1", ylab = "PC2", main = "PCA")
#' par(mfrow = c(1, 1))
#' }
#'
#'
#' #Example 2: submit a block structure which was identified by any other approach. This can be
#' #done by transforming the block structure to an object of type 'blocks' using detect.blocks():
#'
#' if (requireNamespace("glasso", quietly = TRUE)) {
#' #Simulate a data matrix X with a block diagonal population covariance matrix.
#' set.seed(1)
#' n <- 100
#' p <- 4
#' Sigma <- bdsvd.cov.sim(p = p, b = 2, design = "c")
#'
#' X <-  matrix(rnorm(n * p), n, p) %*% chol(Sigma)
#' S <- cov(X)
#'
#' #Identify the block structure using glasso()
#' S.block <- glasso::glasso(S, 0.2)$w
#'
#' #S.blocks is a block diagonal matrix:
#' print(S.block != 0)
#' #We know extract the block information to an object of class 'blocks' using detect.blocks()
#' block.structure <- detect.blocks(S.block)
#' class(block.structure)
#'
#' #The block.structure of class 'blocks' can now be supplied to ispca()
#' ispca(X, block.structure = block.structure, verbose = FALSE)
#' }
#'
#' @importFrom irlba svdr
#'
#' @export
ispca <- function(X,
                  K = 1,
                  block.structure,
                  anp = "2",
                  covariance = TRUE,
                  method = "CDM",
                  orthogonal = FALSE,
                  standardize = FALSE,
                  verbose = TRUE
) {

  p <- ncol(X)

  if (K >= p) {
    warning("K >= p. K is now assigned to p - 1 (K <- p - 1).")
    K <- p - 1
  }

  METHOD <- c("CDM", "DM", "PCA")
  if (!(method %in% METHOD))
    stop(sprintf("%s is an invalid option for method.", method))

  if (is.null(colnames(X))) {
    colnames(X) <- as.character(seq_len(p))
  }



  if (verbose) {
    cat("### Run IS-SPCA:\n")
  }
  if (missing(block.structure)) {
    if (verbose) {
      cat("### 1. Perform BD-SVD\n")
    }

    X <- scale(X, center = TRUE, scale = !covariance)
    block.structure <- bdsvd(X, anp = anp, standardize = standardize, verbose = verbose)
    X.bd <- bdsvd.structure(X, block.structure, output = "submatrices")

  } else {
    if (verbose) {
      cat("### 1. BD-SVD object supplied\n")
    }

    if (!inherits(block.structure, c("bdsvd", "blocks"))) {
      stop(sprintf(
        "block.structure must be a 'bdsvd' or 'blocks' object (got: '%s').\nE.g., pass the result of bdsvd(), single.bdsvd(), or detect.blocks().",
        paste(class(block.structure), collapse = "/")
      ), call. = FALSE)
    }

    if (inherits(block.structure, "blocks")) {
      block.structure <- lapply(block.structure, `[[`, "features")
      class(block.structure) <- "bdsvd"
    }

    block.structure.names <- unlist(block.structure, use.names = FALSE)
    if (!identical(sort(colnames(X)), sort(block.structure.names))) {
      warning("Column names of X do not match the variable names in 'block.structure'.", call. = FALSE)
    }

    block.order <- seq_along(block.structure)
    X.bd <- list()
    for (i in block.order) {
      X.bd <- c(X.bd, list(X[, colnames(X) %in% block.structure[[block.order[i]]], drop = FALSE]))
    }
  }

  b <- length(X.bd)
  if (b == 1) {
    stop("BD-SVD identified only a single block thus no underlying sparse structure.")
  }

  if (verbose | missing(block.structure)) {
    cat("\r                                                                            ")

  cat("\r###    -->", b,"blocks detected\n")
  cat("\r### 2. Perform PCA\n")
  }

  p.i <- cumsum(sapply(X.bd, ncol))
  v <- matrix(0, p, K * b, dimnames = list(unlist(lapply(X.bd, colnames)), NULL))
  l <- numeric(K * b)
  idx.col <- 1

  for (block in seq_len(b)) {
    pca.i <- calc.pca(X = X.bd[[block]], K = K, method = method, covariance = covariance, orthogonal = orthogonal)

    idx.col.next <- idx.col +  pca.i$K - 1
    idx.row <- seq(c(1, p.i[-length(p.i)] + 1)[block], p.i[block])

    l[idx.col:idx.col.next] <- pca.i$l
    v[idx.row, idx.col:idx.col.next] <- pca.i$v
    idx.col <- idx.col.next + 1
  }

  order <- order(l, decreasing = TRUE)
  v <- v[colnames(X), order]
  l <- l[order]

  not.empty <- l != 0
  l <- l[not.empty]
  v <- v[, not.empty]

  tot.var <- ifelse(covariance, sum(apply(X, 2, var)), p)
  exp.var <- l/tot.var

  out <- list(v = v, l = l, exp.var = exp.var, X.b = X.bd, block.structure = block.structure)
  return(structure(out, class = "ispca"))
}




#' @title Principal (Sub)Matrices
#'
#' @description Identifies the principal (sub)matrices
#'
#' @param X Data matrix of dimension \eqn{n}x\eqn{p} with possibly \eqn{p >> n}.
#'
#' @param block.structure Underlying block structure. Must be a '\code{bdsvd}', '\code{blocks}', or '\code{ispca}' object. E.g., pass the result of \code{bdsvd()}, \code{single.bdsvd()},
#' \code{ispca()}, or \code{detect.blocks()}. An identified block structure by any other method can be supplied using \code{detect.blocks()} (see example below).
#'
#' @param rule Which rule should be used to choose principal submatrices. \code{rule = "cumvar"} selects the smallest
#' number of principal submatrices (ordered by explained variance) whose cumulative share is at least \code{value} (a proportion in (0, 1]).
#' \code{rule = "enrich"} selects all principal submatrices that explain \eqn{\times}\code{value} more than they should on average (see \code{value}).
#'
#' @param value Numeric parameter used by \code{rule}.
#' If \code{rule = "cumvar"}, \code{value} is the target cumulative proportion of explained variance (must be in \code{(0, 1]}).
#' Default is \code{0.8}.
#' If \code{rule = "enrich"}, \code{value} is the factor necessary to be selected compared to the equal-share baseline.
#' E.g., if a submatrix should on average explain 10% of the total explained variance and if \code{value = 2}, this submatrix
#' is only selected if it explains at least 2\eqn{x}10% = 20% of the total explained variance. Default is \code{2}.
#'
#' @details
#' This function selects the principal (sub)matrices as described in Bauer (2026).
#'
#' @return
#' A named list with the following components:
#' \describe{
#'   \item{prmats}{List of submatrices ordered by explained variance (rule = 'cumvar') or by factor (rule = 'enrich').
#'                 Each element \code{prmats[[b]]} is a named list with:
#'      \describe{
#'       \item{expl.var}{
#'       Proportion of total variance explained by block \eqn{b}.
#'       }
#'       \item{avg.var}{
#'       Average variance of the variables in block \eqn{b}.
#'       }
#'       \item{factor}{
#'       Enrichment factor \code{expl.var / avg.var} (see \code{value} argument of the function).
#'       }
#'       \item{feature.names}{
#'       Column names (variables) that belong to block \eqn{b}.
#'       }
#'       \item{p.b}{
#'       Number of variables in block \eqn{b}.
#'       }
#'      }
#'   }
#'
#'   \item{X.pr}{
#'   The data matrix of the kept submatrices/variables.
#'   }
#' }
#'
#' @section Access:
#' Submatrices can be accessed with list indexing, e.g., \code{res$prmats[[1]]$feature.names} gives the variable names of the first submatrix.
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{ispca}}
#'
#' @references \cite{Bauer, J.O. (2025). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat., 34(3), 1005–1016}
#' @references \cite{Bauer, J.O. (2026). Beyond regularization: inherently sparse principal component analysis. Stat. Comp.}
#'
#' @examples
#' #Example: principal submatrices of a gene expression data set with two tissue types
#'
#' if (requireNamespace("dslabs", quietly = TRUE)) {
#' data("tissue_gene_expression", package = "dslabs")
#'
#' #We only select the two tissue types kidney (6) and liver (7)
#' Y <- as.numeric(tissue_gene_expression$y)
#' X <- scale(tissue_gene_expression$x[Y %in% c(6, 7), ], scale = FALSE)
#' Y <- Y[Y %in% c(6, 7)]
#'
#'
#' #First: run IS-PCA (or submit a identified block structure using bdsvd(...) or detect.blocks(...))
#'
#' ispca.obj <- ispca(X = X, anp = "1")
#'
#'
#' #Second: extract the submatrices that explain at least 80% (default value) of the total variance
#'
#' res <- prmats(X, block.structure = ispca.obj)
#' res
#'
#' #One submatix is selected which contains 236 variables (out of 500) and explains
#' #81.67% of the total variance
#' length(res$prmats)
#' res$prmats[[1]]$p.b
#' round(res$prmats[[1]]$expl.var * 100, 2)
#'
#'
#' #Alternatively: extract the submatrices that explain five times more of the total variance
#' #than they should on average ('factor')
#'
#' res <- prmats(X, block.structure = ispca.obj, rule = "enrich", value = 1.5)
#' res
#'
#' #The highest 'factor' is 1.73
#' res <- prmats(X, block.structure = ispca.obj, rule = "enrich", value = 2)
#'
#'
#' }
#'
#' @export
prmats <- function(X,
                  block.structure,
                  rule = "cumvar",
                  value
) {

  p <- ncol(X)

  RULE <- c("cumvar", "enrich")
  if (!(rule %in% RULE))
    stop(sprintf("%s is an invalid option for rule", rule))

  if (is.null(colnames(X))) {
    colnames(X) <- as.character(seq_len(p))
  }


  if (missing(value)) {
    value <- if (rule == "cumvar") 0.8 else 2
  } else {
    if (rule == "cumvar" && (value <= 0 || value > 1))
      stop("For rule = 'cumvar', `value` must be in (0, 1].", call. = FALSE)

    if (rule == "enrich" && value <= 0)
      stop("For rule = 'enrich', `value` must be positive.", call. = FALSE)
  }

  CLASS <- c("bdsvd", "blocks", "ispca")
  if (!inherits(block.structure, CLASS))
    stop(sprintf(
      "block.structure must be a 'bdsvd', 'blocks', or 'ispca' object (got: '%s').\nE.g., pass the result of bdsvd(), single.bdsvd(), ispca(), or detect.blocks().",
      paste(class(block.structure), collapse = "/")
    ), call. = FALSE)

  if (inherits(block.structure, "blocks")) {
    block.structure <- lapply(block.structure, `[[`, "features")
    class(block.structure) <- "bdsvd"
  }

  if (inherits(block.structure, "ispca")) {
    block.structure <- block.structure$block.structure
    class(block.structure) <- "bdsvd"
  }

  block.structure.names <- unlist(block.structure, use.names = FALSE)
  if (!setequal(colnames(X), block.structure.names)) {
    warning("Column names of X do not match the variable names in 'block.structure'. It is recommended to check that the submitted block.structure corresponds to X.", call. = FALSE)
  }

  block.vars <- list()
  vars <- apply(X, 2, var)
  total.variance <- sum(vars)

  block.vars <- lapply(seq_along(block.structure), function(i) {
    feature.names <- block.structure[[i]]
    p.b <- length(feature.names)
    expl   <- (sum(vars[feature.names]) / total.variance)
    avg    <- (length(feature.names) / p)
    list(expl.var = expl,
         avg.var  = avg,
         factor   = expl / avg,
         feature.names = feature.names,
         p.b = p.b)
  })

  if (rule == "enrich") {
    act.factors <- vapply(block.vars, `[[`, numeric(1), "factor")
    block.idx <- which(act.factors >= value)
    if (!any(block.idx)) {
      warning(sprintf("No submatrix has a factor of %s or higher. The largest factor is %s (rounded to two decimals).\n  You might lower the 'value' argument and re-run prmats().",
                      value, round(max(act.factors), 2)), call. = FALSE)
      block.vars <- list(expl.var = NULL, factor = NULL, feature.names = NULL)
      return(list(prmats = block.vars, X.pr = NULL))
    }
    ord <- block.idx[order(act.factors[block.idx], decreasing = TRUE)]
    block.vars <- block.vars[ord]
    X.pr <- bdsvd.structure(X, block.structure, output = "matrix", block.order = ord)

    out <- list(prmats = block.vars, X.pr = X.pr)
    return(structure(out, class = "prmats"))
  }

  if (rule == "cumvar") {
    value <- value
    act.expl.vars <- vapply(block.vars, `[[`, numeric(1), "expl.var")
    ord.vars <- order(act.expl.vars, decreasing = TRUE)
    if (value == 1) {
      ord <- ord.vars
    } else {
      ord <- ord.vars[order(act.expl.vars, decreasing = TRUE)[seq_len(min(which(cumsum(act.expl.vars) >= value)))]]
    }
    block.vars <- block.vars[ord]
    X.pr <- bdsvd.structure(X, block.structure, output = "matrix", block.order = ord)

    out <- list(prmats = block.vars, X.pr = X.pr)
    return(structure(out, class = "prmats"))
  }

}








#' @export
print.ispca <- function(x, ...) {
  K <- min(length(x$l), 10)
  tab <- data.frame(
    PC                = seq_len(K),
    "Eigenvalue"      = formatC(x$l[seq_len(K)], digits = 2, format = "f"),
    "Expl.Var."       = formatC(x$exp.var[seq_len(K)] * 100,    digits = 2, format = "f"),
    "Cumul.Expl.Var." = formatC(cumsum(x$exp.var[seq_len(K)]) * 100,   digits = 2, format = "f"),
    row.names = NULL
  )

  cat("Number of identified blocks:", length(x$block.structure), "\n")
  cat(sprintf("First %d principal components:\n", K))
  print(tab, row.names = FALSE)
  cat("\nAvailable components:\n")
  cat(paste0("$", names(x), collapse = "\n"), "\n")
  invisible(x)
}



#' @export
print.prmats <- function(x, ...) {
  kept.expl.var <- round(sum(vapply(x$prmats, `[[`, numeric(1), 1)) * 100, 2)
  cat(sprintf("Number of kept submatrices: %s\n", length(x$prmats)))
  cat(sprintf("Total variables kept: %s\n", ncol(x$X.pr)))
  cat(sprintf("Total explained variance of kept submatrices: %s%%\n\n", kept.expl.var))

  cat("Submatrices are ordered by explained variance (rule = 'cumvar') ",
      "or by factor (rule = 'enrich'). Access them via $prmats.\n", sep = "")
  cat("The data matrix of the kept submatrices/variables is available via $X.pr.\n")
  invisible(x)
}




