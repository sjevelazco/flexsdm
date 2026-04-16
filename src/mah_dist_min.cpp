#include <Rcpp.h>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mah_dist_min(const NumericMatrix& x, const NumericMatrix& y, const NumericMatrix& cov) {
  int nx = x.nrow();
  int ny = y.nrow();
  int ncol = x.ncol();

  if (ncol != y.ncol()) {
    stop("Matrices 'x' and 'y' must have the same number of columns.");
  }
  if (ncol != cov.nrow() || ncol != cov.ncol()) {
    stop("The covariance matrix 'cov' must be square and match the column dimension of 'x' and 'y'.");
  }

  if (nx == 0) return NumericVector(0);
  if (ny == 0 || ncol == 0) {
    return NumericVector(nx, NA_REAL);
  }

  // Get the inverse of the covariance matrix by calling R's solve() once.
  // This is computationally efficient and handles matrix inversion robustly.
  Function solve("solve");
  NumericMatrix inv_cov = solve(cov);

  const double* ptr_x = x.begin();
  const double* ptr_y = y.begin();
  const double* ptr_inv = inv_cov.begin();

  // Result vector: minimum distance for each row in x
  NumericVector min_dist(nx);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nx; i++) {
    // Check if row i of x has NA
    bool x_has_na = false;
    for (int k = 0; k < ncol; k++) {
      if (NumericVector::is_na(ptr_x[i + (size_t)k * nx])) {
        x_has_na = true;
        break;
      }
    }

    if (x_has_na) {
      min_dist[i] = NA_REAL;
      continue;
    }

    double best_d2 = R_PosInf;
    // Temporary buffer to store (x_i - y_j) to speed up quadratic form calculation
    std::vector<double> delta(ncol);

    for (int j = 0; j < ny; j++) {
      bool y_has_na = false;
      for (int k = 0; k < ncol; k++) {
        double val_y = ptr_y[j + (size_t)k * ny];
        if (NumericVector::is_na(val_y)) {
          y_has_na = true;
          break;
        }
        delta[k] = ptr_x[i + (size_t)k * nx] - val_y;
      }
      
      if (y_has_na) continue;

      double d2 = 0.0;
      // Compute the quadratic form: delta^T * inv_cov * delta
      for (int c = 0; c < ncol; c++) {
        double row_sum = 0.0;
        for (int r = 0; r < ncol; r++) {
          row_sum += delta[r] * ptr_inv[r + (size_t)c * ncol];
        }
        d2 += row_sum * delta[c];
      }

      if (d2 < best_d2) {
        best_d2 = d2;
      }
    }
    
    // Store square root of the minimum squared distance
    min_dist[i] = (best_d2 == R_PosInf) ? NA_REAL : std::sqrt(std::max(0.0, best_d2));
  }

  return min_dist;
}