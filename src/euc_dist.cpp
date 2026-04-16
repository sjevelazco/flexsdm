#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix euc_dist(const NumericMatrix& x, const NumericMatrix& y) {
  int nx = x.nrow();
  int ny = y.nrow();
  int ncol = x.ncol();

  if (ncol != y.ncol()) {
    stop("Matrices 'x' and 'y' must have the same number of columns.");
  }

  NumericMatrix result(nx, ny);

  if (nx == 0 || ny == 0 || ncol == 0) {
    return result;
  }

  const double* ptr_x = x.begin();
  const double* ptr_y = y.begin();
  double* ptr_res = result.begin();

  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static)
  #endif
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double sum_sq = 0.0;
      bool has_na = false;
      
      for (int k = 0; k < ncol; k++) {
        double val_x = ptr_x[i + (size_t)k * nx];
        double val_y = ptr_y[j + (size_t)k * ny];
        
        if (NumericVector::is_na(val_x) || NumericVector::is_na(val_y)) {
          has_na = true;
          break;
        }
        
        double diff = val_x - val_y;
        sum_sq += diff * diff;
      }
      
      ptr_res[i + (size_t)j * nx] = has_na ? NA_REAL : std::sqrt(sum_sq);
    }
  }

  // Handle dimnames according to the R implementation logic
  RObject x_dn = x.attr("dimnames");
  RObject y_dn = y.attr("dimnames");
  SEXP rn = (!Rf_isNull(x_dn) && as<List>(x_dn).size() > 0) ? as<List>(x_dn)[0] : R_NilValue;
  SEXP cn = (!Rf_isNull(y_dn) && as<List>(y_dn).size() > 0) ? as<List>(y_dn)[0] : R_NilValue;
  result.attr("dimnames") = List::create(rn, cn);

  return result;
}
