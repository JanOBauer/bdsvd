// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// Soft-thresholding function
// [[Rcpp::export]]
arma::mat soft_cpp(const arma::mat& X, const arma::mat& u, int q, int n) {
  arma::mat y = X.t() * u;
  arma::mat a = abs(y);
  arma::mat z = sort(a, "ascend", 0);

  arma::rowvec lambda(q);
  for (int j = 0; j < q; j++) {
    lambda(j) = z(n - 1, j);
  }

  a.each_row() -= lambda;
  a.transform([](double val) { return std::max(val, 0.0); });

  return sign(y) % a;
}




// Calculate First Sparse Loading
// [[Rcpp::export]]
Rcpp::List calc_one_sparse_v_cpp(const arma::mat& X, arma::vec v, arma::vec u, int n, int maxit = 500, double tol = 0.001) {

  if (n == 0) {
    //warning("n = 0: no sparsity constrains specified");
    return Rcpp::List::create(
      Rcpp::Named("v") = v,
      Rcpp::Named("u") = u
    );
  }

  int iter = 0;
  double delta = std::numeric_limits<double>::infinity();

  while (delta > tol && iter < maxit) {
    arma::vec u_check = u;
    v = soft_cpp(X, u, 1, n);
    u = X * v;
    u = u / arma::norm(u, 2);

    delta = abs(1 - abs(arma::dot(u_check, u)));

    ++iter;
  }

  if (iter >= maxit) {
    Rcpp::Rcerr << "Warning: Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.\n";
  }

  v = v / arma::norm(v, 2);
  //u = X * v;

  return Rcpp::List::create(
    Rcpp::Named("v") = v,
    Rcpp::Named("u") = u
  );
}




// Calculate the BIC for BD-SVD
// [[Rcpp::export]]
Rcpp::NumericVector calc_BIC(const arma::mat& X, arma::vec u, arma::vec v, const arma::vec& dof_grid,
                             int maxit, int anp, int n, int p) {

  int i = 0;
  int grid_size = dof_grid.n_elem;
  Rcpp::NumericVector BIC(grid_size);

  double n_BIC = n * p;

  double a_np = 1.0;
  switch (anp) {
  case 2:
    a_np = 0.5 * std::log(n_BIC);
    break;
  case 3:
    a_np = std::log(std::log(n_BIC));
    break;
  case 4:
    a_np = std::log(std::log(p));
    break;
  }
  double log_np = log(n_BIC) / (n_BIC);

  for (i = 0; i < grid_size; ++i) {
    double dof = dof_grid(i);

    Rcpp::List sparse_uv = calc_one_sparse_v_cpp(X, v, u, dof, maxit);
    u = Rcpp::as<arma::vec>(sparse_uv["u"]);
    v = Rcpp::as<arma::vec>(sparse_uv["v"]);

    u = X * v;
    double res_norm = norm(X - u * v.t(), "fro");
    double BIC_value = log(res_norm * res_norm / (n_BIC)) + (p - dof) * log_np * a_np;

    BIC(i) = BIC_value;
  }

  return BIC;
}

