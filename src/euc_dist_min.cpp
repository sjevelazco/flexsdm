#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector euc_dist_min(const NumericMatrix& x, const NumericMatrix& y) {
  int nx = x.nrow();
  int ny = y.nrow();
  int ncol = x.ncol();

  // Basic validation
  if (ncol != y.ncol()) {
    stop("Matrices 'x' and 'y' must have the same number of columns.");
  }

  // Check for empty inputs
  if (nx == 0) return NumericVector(0);
  if (ny == 0 || ncol == 0) {
    return NumericVector(nx, NA_REAL);
  }

  // Result vector: minimum distance for each row in x
  NumericVector min_dist(nx);

  // Use raw pointers for thread safety and faster column-major access
  const double* ptr_x = x.begin();
  const double* ptr_y = y.begin();

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int i = 0; i < nx; i++) {
    double best_sq = R_PosInf;
    
    for (int j = 0; j < ny; j++) {
      double sum_sq = 0.0;
      bool has_na = false;
      
      for (int k = 0; k < ncol; k++) {
        // Column-major: index = row + col * nrows
        double val_x = ptr_x[i + (size_t)k * nx];
        double val_y = ptr_y[j + (size_t)k * ny];
        
        if (NumericVector::is_na(val_x) || NumericVector::is_na(val_y)) {
          has_na = true;
          break;
        }
        
        double diff = val_x - val_y;
        sum_sq += diff * diff;
        
        // Early exit: if current sum_sq is already worse than best_sq
        if (sum_sq >= best_sq) break;
      }
      
      if (!has_na && sum_sq < best_sq) {
        best_sq = sum_sq;
      }
    }
    
    min_dist[i] = (best_sq == R_PosInf) ? NA_REAL : std::sqrt(best_sq);
  }

  return min_dist;
}