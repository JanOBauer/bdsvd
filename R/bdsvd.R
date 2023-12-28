#' @importFrom methods new
create.block <- function(feature.names, selected.features, block.columns) {
  if (length(feature.names) > 0) {
    selected.features <- feature.names[selected.features]
  }
  return(new("Block", features = selected.features, block.columns = block.columns))
}

get.blocks <- function(threshold.matrix, feature.names){
  p <- nrow(threshold.matrix)
  k <- ncol(threshold.matrix)
  columns <- 1:k
  blocks <- list()

  if (any(colSums(threshold.matrix == 0) == p)) {
    stop("Zero column. Reduce regularization and/or threshold value.\n")
  }

  if (any(rowSums(threshold.matrix == 0) == k)) {
    if(k < p){
      stop("Zero row. Add more loadings, or reduce regularization and/or threshold value.\n")
    }
    else{
      stop("Zero row. Reduce regularization and/or threshold value.\n")
    }

  }

  while(length(columns) != 0){
    block.columns <- columns[1]
    while(!identical(block.columns, find.blocks(threshold.matrix, block.columns))){
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

get.feature.names <- function(X) {
  feature.names <- colnames(X)
  if (length(feature.names) <= 0) {
    feature.names <- seq_len(ncol(X))
  }

  return(feature.names)
}

get.threshold.matrix <- function(loadings, threshold){

  loadings[which(abs(loadings) <= threshold)] = 0
  loadings[which(abs(loadings) != 0)] = 1

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
setClass("Block", slots = c(features = "vector", block.columns = "vector"))



#' @title Covariance matrix simulation
#'
#' @description This function generates covariance matrices based on the simulation studies described in Bauer (202X).
#'
#' @param p Number of variables.
#'
#' @param b Number of blocks.
#'
#' @param design Simulation design "a", "b", "c", or "d".
#'
#' @return
#' A list containing the feature names of all detected blocks.
#'
#' @references
#' Bauer, J.O. (202X). \emph{High-dimensional block diagonal covariance structure detection using singular vectors}\cr
#'
#' @examples
#' #The covariance matrix for simulation design (a) is given by
#' Sigma <- cov.sim(p = 500, b = 500, design = "a")
#'
#' @importFrom stats runif
#'
#' @export
cov.sim <- function(p = p,
                    b = b,
                    design = design
                    ){
  if(design == "a"){
    Sigma <- diag(1, p, p)
    return(Sigma)
  }

  if(design == "b"){
    Sigma <- diag(runif(p, 1, 5), p, p)
    return(Sigma)
  }

  if(design == "c"){
    rho0 <- 0.2
    epsilon <- 0.1

    d <- p/b
    Sigma <- matrix(0, p, p)
    for(i in 1:b){
      rho <- runif(1, rho0-epsilon, rho0+epsilon)
      Sigma[(1 + d*(i-1)):(d + d*(i-1)), (1 + d*(i-1)):(d + d*(i-1))] <-
        ( (1-rho) * diag(1,d, d) + 2* rho * rep(1,d)%*%t(rep(1,d))   )
    }
    return(Sigma)
  }

  if(design == "d"){
    rho0 <- 0.45
    epsilon <- 0.15
    omega <- 0.1

    d <- p/b
    Sigma <- matrix(0, p, p)
    for(B in 1:b){
      rho <- runif(1, rho0-epsilon, rho0+epsilon)
      Rii <- matrix(0, d, d)
      for (i in 1:d) {
        for (j in 1:d) {
          Rii[i, j] <- (-1)^(i + j) * rho^(abs(i - j)^(omega))
        }
      }
      diag(Rii) <- 1
      Sigma[(1 + d*(B-1)):(d + d*(B-1)), (1 + d*(B-1)):(d + d*(B-1))] <- Rii
    }
    return(Sigma)
  }

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
#' An object of class Block containing the features and columns indices corresponding to each detected block.
#'
#' @seealso \link{bdsvd}, \link{iterative.bdsvd}
#'
#' @references
#' Bauer, J.O. (202X). \emph{High-dimensional block diagonal covariance structure detection using singular vectors}\cr
#'
#' @examples
#' #In the first example, we replicate the simulation study for the ad hoc procedure
#' #Est_0.1 in Bauer (202X). In the second example, we manually compute the first step
#' #of BD-SVD, provided in the functions iterative.bdsvd() and/or bdsvd(), for constructed
#' #sparse loadings
#'
#' #Example 1: Replicate the simulation study (a) from Bauer (202X) for the ad hoc procedure Est_0.1.
#' require(mvtnorm)
#'
#' p <- 500 #Number of variables
#' n <- 125 #Number of observations
#' b <- 500 #Number of blocks
#' design <- "a"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- cov.sim(p = p, b = b, design = design)
#' X <- rmvnorm(n, mean=rep(0,p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' #Perform the ad hoc procedure
#' require(cvCovEst)
#' detect.blocks(scadEst(dat = X, lambda = 0.2), threshold = 0)
#'
#' #Example 2: Manually compute the first step of BD-SVD
#' #for some loadings V that mirror the two blocks
#' #("A", "B") and c("C", "D").
#' require(mvtnorm)
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

  if(missing(V)){stop("V is required.")}

  if(length(rownames(V)) == 0){rownames(V) == as.character(1:nrow(V))}

  threshold.matrx <- get.threshold.matrix(V, threshold)

  return(get.blocks(threshold.matrx, rownames(V)))
}



#' @title Block detection using singular vectors (BD-SVD).
#'
#' @description Performs a single iteration of BD-SVD: splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (sparse
#' approximation of the first right singular vector, i.e., has many zero values) that mirrors the
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
#' @seealso \link{bdsvd.ht}, \link{iterative.bdsvd}
#'
#' @references
#' Bauer, J.O. (202X). \emph{High-dimensional block diagonal covariance structure detection using singular vectors}\cr
#' Shen, H. and Huang, J.Z. (2008). \emph{Sparse principal component analysis via regularized low rank matrix approximation}, \emph{J. Multivar. Anal. 99, 1015–1034.}, <\doi{10.1016/j.jmva.2007.06.007}>\cr
#'
#' @examples
#' #Replicate the illustrative example from Bauer (202X).
#'
#' require(mvtnorm)
#'
#' p <- 300 #Number of variables. In Bauer (202X), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- cov.sim(p = p, b = b, design = design)
#' X <- rmvnorm(n, mean=rep(0,p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' @importFrom irlba ssvd
#'
#' @export
bdsvd <- function(X,
                  dof, #number of nonzero components
                  standardize = TRUE,
                  max.iter
                  ) {

  if(anyNA(X)){stop("X contains missing value indicator (NA)")}

  p <- ncol(X)

  if(missing(dof)){stop("Enter the degrees of freedom")}

  if(missing(max.iter)){max.iter <- 500}

  X <- scale(X, center = TRUE, scale = standardize)
  feature.names <- get.feature.names(X = X)
  eigen <- list()

  v <- tryCatch(suppressWarnings(ssvd(x = X, k = 1, n = p - dof, maxit = max.iter)$v), error = function(e) e)
  if (inherits(v, "error")) {
    stop("dof = ", dof, " do not fit the structure of the singular vectors. You can use bdsvd.ht to find a suitable value for dof.")
  }

  if(length(which(v == 0)) == 0){return(list(feature.names))}

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
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package. The computation of the HBIC is outlined in Bauer (202X).
#'
#' @return
#' \item{dof}{
#'   The optimal number of nonzero components (degrees of freedom) according to the HBIC.
#' }
#' \item{BIC}{
#'   The HBIC for the different numbers of nonzero components.
#' }
#'
#' @seealso \link{bdsvd}, \link{iterative.bdsvd}
#'
#' @references
#' Bauer, J.O. (202X). \emph{High-dimensional block diagonal covariance structure detection using singular vectors}\cr
#' Shen, H. and Huang, J.Z. (2008). \emph{Sparse principal component analysis via regularized low rank matrix approximation}, \emph{J. Multivar. Anal. 99, 1015–1034}, <\doi{10.1016/j.jmva.2007.06.007}>\cr
#'
#' @examples
#' #Replicate the illustrative example from Bauer (202X).
#'
#' require(mvtnorm)
#'
#' p <- 300 #Number of variables. In Bauer (202X), p = 3000
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- cov.sim(p = p, b = b, design = design)
#' X <- rmvnorm(n, mean=rep(0,p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' @importFrom irlba ssvd
#'
#' @export
bdsvd.ht <- function(X,
                     dof.lim,
                     standardize = TRUE,
                     max.iter
                     ) {

  if(anyNA(X)){stop("X contains missing value indicator (NA)")}

  n <- nrow(X)
  p <- ncol(X)

  if(missing(dof.lim)){ dof.lim <- c(0,p-1) }
  dof.lim <- sort(round(dof.lim))
  if(dof.lim[1] < 0){ dof.lim[1] <- 0 }
  if(dof.lim[2] > p-1){ dof.lim[2] <- p-1 }
  dof.grid <- dof.lim[1]:dof.lim[2]

  X <- scale(X, center = TRUE, scale = standardize)

  if(missing(max.iter)){max.iter <- 500}

  BIC <- vector(length = length(dof.grid))
  i <- 1
  for(dof in dof.grid){

    v <- tryCatch(suppressWarnings(ssvd(x = X, k = 1, n = p - dof, maxit = max.iter)$v), error = function(e) e)
    if (inherits(v, "error")) {
      v <- matrix(0, nrow = p, ncol = 1)
    }

    u <- X%*%v
    BIC[i] <- log(norm(X - u%*%t(v), type = "F")^2/n/p) + sum(v!=0) * log(n*p)^(2/3)* log(log(n*p))/n/p * log(p)
    i <- i + 1
  }

  dof <- dof.grid[order(BIC)[1]]
  BIC <- cbind.data.frame(BIC)
  rownames(BIC) <- dof.grid

  return(list(dof = dof, BIC = BIC))
}



#' @title Iterative Block Detection Using Singular Vectors (BD-SVD).
#'
#' @description Performs BD-SVD iteratively to reveal the block structure. splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (sparse
#' approximation of the first right singular vector, i.e., has many zero values) that mirrors the
#' shape of the covariance matrix. This procedure is continued iteratively until the block diagonal structure has been revealed.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param dof.lim Interval limits for the number of non-zero components in the sparse loading (degrees of freedom).
#' If \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom. Default is \code{dof.lim <- c(0, p-1)} which is highly recommended to check for
#' all levels of sparsity.
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
#' @seealso \link{bdsvd}, \link{bdsvd.ht}
#'
#' @references
#' Bauer, J.O. (202X). \emph{High-dimensional block diagonal covariance structure detection using singular vectors}\cr
#'
#' @examples
#' #Replicate the simulation study (c) from Bauer (202X).
#'
#' require(mvtnorm)
#'
#' p <- 500 #Number of variables
#' n <- 250 #Number of observations
#' b <- 10  #Number of blocks
#' design <- "c" #Simulation design "a", "b", "c", or "d".
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- cov.sim(p = p, b = b, design = design)
#' X <- rmvnorm(n, mean=rep(0,p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' iterative.bdsvd(X, standardize = FALSE)
#'
#' @importFrom irlba ssvd
#'
#' @export
iterative.bdsvd <- function(X,
                           dof.lim,
                           standardize = TRUE,
                           max.iter,
                           trace = TRUE
                           ) {

  if(anyNA(X)){stop("X contains missing value indicator (NA)")}

  n <- nrow(X)
  p <- ncol(X)

  if(length(colnames(X)) == 0){colnames(X) <- as.character(1:p)}

  if(missing(dof.lim)){ dof.lim <- c(0,p-1) }
  dof.lim <- sort(round(dof.lim))
  if(dof.lim[1] < 0){ dof.lim[1] <- 0 }
  if(dof.lim[2] > p-1){ dof.lim[2] <- p-1 }

  if(missing(max.iter)){max.iter <- 500}

  sub.matrices <- list(colnames(X))
  results <- list()

  b <- 1
  while (TRUE) {
    for (i in 1:length(sub.matrices)) {

      if (length(sub.matrices) == 0) {
        break  #Stop when splitting is no longer supported
      }

      if(length(sub.matrices[[1]]) == 1){
        if(trace){
          cat("Block", b,":", sub.matrices[[1]],"\n",sep="\t")
        }
        b <- b+1
        results <- c(results, sub.matrices[[1]])
        sub.matrices <- sub.matrices[-1]
        next
      }

      dof.split <- bdsvd.ht(X = X[, colnames(X) %in% sub.matrices[[1]] ],
                           dof.lim = dof.lim,
                           standardize = standardize,
                           max.iter = max.iter)$dof
      sub.results <- bdsvd(X = X[, colnames(X) %in% sub.matrices[[1]] ], standardize = standardize, dof = dof.split)


      if(length(sub.matrices) == 1){
        sub.matrices <- list()
      } else{
        sub.matrices <- sub.matrices[-1]
      }

      if(length(sub.results) == 1){
        if(trace){
          cat("Block", b,":", sub.results[[1]],"\n",sep="\t")
        }
        b <- b+1
        results <- c(results, sub.results)
      }else{
        sub.matrices <- c(sub.matrices, sub.results)
      }

    }
    if (length(sub.matrices) == 0) {
      break  #Stop when splitting is no longer supported
    }

  }

  class(results) <- "bdsvd"
  return(results)
}





#' @title Data Matrix Structure According to the Detected Block Structure.
#'
#' @description Sorts the data matrix \eqn{X} according to the detected block structure \eqn{X_1 , ... , X_b} and returns
#' the detected submatrices individually.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param block.structure A user supplied block structure based on \code{iterative.bdsvd()} or \code{bdsvd()}.
#'
#' @param result Should the result be the data matrix ordered according to the blocks (\code{result == "matrix"}), or
#' a list containing the submatrices (\code{result == "submatrices"}). Default is \code{matrix}.
#'
#' @param block.order A vector that contains the order of the blocks detected by \code{iterative.bdsvd()} or \code{.bdsvd()}.
#' The vector must contain the index of each blocks exactly once. Default is \code{1:b} where \code{b} is the total number of blocks.
#'
#' @return
#' Either the data matrix \code{X} with columns sorted according to the detected blocks, or a list containing the detected
#' submatrices.
#'
#' @seealso \link{iterative.bdsvd}, \link{bdsvd}
#'
#' @references
#' Bauer, J.O. (202X). \emph{High-dimensional block diagonal covariance structure detection using singular vectors}\cr
#'
#' @examples
#' #Toying with the illustrative example from Bauer (202X).
#'
#' require(mvtnorm)
#'
#' p <- 300 #Number of variables. In Bauer (202X), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- cov.sim(p = p, b = b, design = design)
#' X <- rmvnorm(n, mean=rep(0,p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' #Compute iterative BD-SVD
#' bdsvd.obj <- iterative.bdsvd(X, standardize = FALSE)
#'
#' #Obtain the data matrix X, sorted by the detected blocks
#' colnames(result.bdsvd(X, bdsvd.obj, result <- "matrix") )
#' colnames(result.bdsvd(X, bdsvd.obj, result <- "matrix", block.order = c(2,1,3)) )
#'
#' #Obtain the detected submatrices X_1, X_2, and X_3
#' colnames(result.bdsvd(X, bdsvd.obj, result <- "submatrices")[[1]] )
#' colnames(result.bdsvd(X, bdsvd.obj, result <- "submatrices")[[2]] )
#' colnames(result.bdsvd(X, bdsvd.obj, result <- "submatrices")[[3]] )
#'
#' @export
result.bdsvd <- function(X,
                         block.structure,
                         result = c("matrix", "submatrices"),
                         block.order
                         ){

  if(class(block.structure) != "bdsvd"){stop("block.structure must be the outcome of bdsvd() or iterative.bdsvd().")}

  b <- length(block.structure)
  ifelse(missing(block.order),
         block.order <- 1:length(block.structure),
         block.order <- as.integer(block.order) )
  if(!identical(sort(block.order), 1:b)){stop("block.order must contain the index of each blocks exactly once.")}

  result <- match.arg(result)
  switch(result,
         "matrix" = {
           output <- X[, colnames(X) %in% block.structure[[block.order[1]]] ]
           for(i in block.order[-1]){
             output <- cbind.data.frame(output,  X[, colnames(X) %in% block.structure[[i]] ])
           }
         },
         "submatrices" = {
           output <- list()
           for(i in block.order){
             output <- c(output, list(X[, colnames(X) %in% block.structure[[block.order[i]]] ]) )
           }
         }
  )

  return(output)
}

