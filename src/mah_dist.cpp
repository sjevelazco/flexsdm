#include <Rcpp.h>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mah_dist(const NumericMatrix& x, const NumericMatrix& y, const NumericMatrix& cov) {
  int nx = x.nrow();
  int ny = y.nrow();
  int ncol = x.ncol();

  if (ncol != y.ncol()) {
    stop("Matrices 'x' and 'y' must have the same number of columns.");
  }
  if (ncol != cov.nrow() || ncol != cov.ncol()) {
    stop("The covariance matrix 'cov' must be square and match the column dimension of 'x' and 'y'.");
  }

  NumericMatrix result(nx, ny);

  if (nx == 0 || ny == 0 || ncol == 0) {
    return result;
  }

  // Get the inverse of the covariance matrix by calling R's solve() once.
  // This is computationally efficient and handles matrix inversion robustly.
  Function solve("solve");
  NumericMatrix inv_cov = solve(cov);

  const double* ptr_x = x.begin();
  const double* ptr_y = y.begin();
  const double* ptr_inv = inv_cov.begin();
  double* ptr_res = result.begin();

  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double d2 = 0.0;
      bool has_na = false;

      // Compute the quadratic form: (x_i - y_j)^T * inv_cov * (x_i - y_j)
      for (int c = 0; c < ncol; c++) {
        double diff_c = ptr_x[i + (size_t)c * nx] - ptr_y[j + (size_t)c * ny];
        if (NumericVector::is_na(diff_c)) {
          has_na = true;
          break;
        }

        double row_sum = 0.0;
        for (int r = 0; r < ncol; r++) {
          double diff_r = ptr_x[i + (size_t)r * nx] - ptr_y[j + (size_t)r * ny];
          row_sum += diff_r * ptr_inv[r + (size_t)c * ncol];
        }
        d2 += row_sum * diff_c;
      }

      // Store sqrt of the squared Mahalanobis distance
      ptr_res[i + (size_t)j * nx] = has_na ? NA_REAL : std::sqrt(std::max(0.0, d2));
    }
  }

  // Handle dimnames: result rows get names from x, result columns get names from rows of y
  RObject x_dn = x.attr("dimnames");
  RObject y_dn = y.attr("dimnames");
  SEXP rn = (!Rf_isNull(x_dn) && as<List>(x_dn).size() > 0) ? as<List>(x_dn)[0] : R_NilValue;
  SEXP cn = (!Rf_isNull(y_dn) && as<List>(y_dn).size() > 0) ? as<List>(y_dn)[0] : R_NilValue;
  result.attr("dimnames") = List::create(rn, cn);

  return result;
}
