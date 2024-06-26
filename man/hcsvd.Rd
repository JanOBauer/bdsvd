% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hcsvd.R
\name{hcsvd}
\alias{hcsvd}
\title{Hierarchical Variable Clustering Using Singular Vectors (HC-SVD).}
\usage{
hcsvd(X, k = "all", linkage = "single", reliability, R, max.iter, trace = TRUE)
}
\arguments{
\item{X}{Data matrix of dimension \eqn{n x p}. The data matrix is standardized during the analysis by \code{hcsvd}.}

\item{k}{Number of sparse loadings to be used. This should be \code{"all"} for all sparse loadings, or \code{"Kaiser"} for as many sparse loadings as
there are eigenvalues larger or equal to one (see Bauer (202Xb) for details). Selecting \code{"Kaiser"} reduces computation time.}

\item{linkage}{The linkage function to be used. This should be one of \code{"average"}, \code{"single"}, or
\code{"RV"} (for RV-coefficient).}

\item{reliability}{By default, the value of each cluster equals the distance calculated by the chosen linkage function.
If preferred, the value of each cluster can be assigned by its reliability. When \code{reliability = spectral}, the reliability is
calculated by the averaged spectral norm.}

\item{R}{Sample correlation matrix of \code{X}. By default, \code{R <- cov(X)}.}

\item{max.iter}{How many iterations should be performed for computing the sparse loadings.
Default is \code{200}.}

\item{trace}{Print out progress as \eqn{p-1} iterations for divisive hierarchical clustering are performed.
Default is \code{TRUE}.}
}
\value{
A list with two components:
\item{dist.matrix}{
The ultrametric distance matrix (cophenetic matrix) of the HC-SVD structure as an object of class \code{dist}.
}
\item{u.cor}{
The ultrametric correlation matrix of \eqn{X} obtained by HC-SVD as an object of class \code{matrix}.
}
}
\description{
Performs HC-SVD to reveal the hierarchical variable structure as descried in Bauer (202Xb). For this divise approach, each cluster is split into two clusters iteratively. Potential splits
are identified by the first sparse loadings (which are sparse approximations of the first right singular vectors, i.e., vectors with many zero values) that
mirror the masked shape of the correlation matrix. This procedure is continued until each variable lies in a single cluster.
}
\details{
The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
the \code{irlba} package.
}
\examples{
#We replicate the simulation study in Bauer (202Xb)
p <- 100
n <- 300
b <- 5
design <- "a"

Rho <- hcsvd.cor.sim(p = 100, b = 5, design = "a")
X <- scale(mvtnorm::rmvnorm(300, mean=rep(0,100), sigma=Rho, checkSymmetry = FALSE))
colnames(X) = 1:ncol(X)
hcsvd.obj <- hcsvd(X, k = "Kaiser")

#The dendrogram can be obtained from the ultrametric distance matrix:
plot(hclust(hcsvd.obj$dist.matrix))


}
\references{
\cite{Bauer, J.O. (202Xb). Hierarchical variable clustering using singular vectors.}

\cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
}
