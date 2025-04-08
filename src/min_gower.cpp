#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>

// Helper function to calculate the range of a numeric vector, with safety check
double range(const NumericVector& vec) {
  if (vec.size() == 0 || vec.size() == 1) {
    return 0.0;
  }
  double min_val = *std::min_element(vec.begin(), vec.end());
  double max_val = *std::max_element(vec.begin(), vec.end());
  return max_val - min_val;
}

// Gower distance function for mixed data types
// [[Rcpp::export]]
NumericVector min_gower_rcpp(DataFrame data1_r, DataFrame data2_r) {
  int n1 = data1_r.nrows();
  int n2 = data2_r.nrows();
  int p = data1_r.size(); // Number of columns

  NumericMatrix distance_matrix(n1, n2);

  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      double sum_sij_wij = 0.0;
      double sum_wij = 0.0;

      for (int k = 0; k < p; ++k) {
        RObject col1 = data1_r[k];
        RObject col2 = data2_r[k];

        if (col1.isNULL() || col2.isNULL()) continue;

        double sij = 0.0;
        double wij = 1.0; // Default weight is 1

        if (Rf_isNumeric(col1)) { // Numeric column
          NumericVector vec1 = as<NumericVector>(col1);
          NumericVector vec2 = as<NumericVector>(col2);

          if (NumericVector::is_na(vec1[i]) || NumericVector::is_na(vec2[j])) continue;

          double Rk = range(vec1);
          if (Rk != 0.0) {
            sij = std::abs(vec1[i] - vec2[j]) / Rk;
          }
        } else { // Categorical column
          CharacterVector vec1 = as<CharacterVector>(col1);
          CharacterVector vec2 = as<CharacterVector>(col2);

          if (vec1[i] != vec2[j]) sij = 1.0;
        }

        sum_sij_wij += sij * wij;
        sum_wij += wij;
      }

      distance_matrix(i, j) = (sum_wij == 0.0) ? NA_REAL : sum_sij_wij / sum_wij;
    }
  }

  NumericVector min_distances(n2, R_PosInf);
  for (int j = 0; j < n2; ++j) {
    for (int i = 0; i < n1; ++i) {
      if (!NumericVector::is_na(distance_matrix(i, j))) {
        min_distances[j] = std::min(min_distances[j], distance_matrix(i, j));
      }
    }
  }

  // Invert values 1 - min_distances
  for (double& min_distances : min_distances) { // Use '&' for reference to modify in-place
    min_distances = std::max(0.0, 1.0 - min_distances); // Calculate 1 - dist, then take max(0, result)
  }
  return min_distances;
}
