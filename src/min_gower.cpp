#include <Rcpp.h>
#include <unordered_map>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

double range_calc(const NumericVector& vec) {
  if (vec.size() == 0 || vec.size() == 1) {
    return 0.0;
  }
  double min_val = *std::min_element(vec.begin(), vec.end());
  double max_val = *std::max_element(vec.begin(), vec.end());
  return max_val - min_val;
}

// [[Rcpp::export]]
NumericVector min_gower_rcpp(DataFrame data1_r, DataFrame data2_r, int n_threads = 0) {
  int n1 = data1_r.nrows();
  int n2 = data2_r.nrows();
  int p = data1_r.size();
  
#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#endif
  
  // Separate numeric and categorical columns
  std::vector<NumericVector> num_cols1, num_cols2;
  std::vector<double> num_ranges;
  
  // For categorical, we'll store as integer codes for faster comparison
  std::vector<std::vector<int>> cat_cols1, cat_cols2;
  std::vector<int> cat_nlevels; // Number of unique categories across both datasets
  
  for (int k = 0; k < p; ++k) {
    RObject col1 = data1_r[k];
    RObject col2 = data2_r[k];
    
    if (Rf_isNumeric(col1)) {
      NumericVector vec1 = as<NumericVector>(col1);
      NumericVector vec2 = as<NumericVector>(col2);
      num_cols1.push_back(vec1);
      num_cols2.push_back(vec2);
      num_ranges.push_back(range_calc(vec1));
    } else {
      // Convert to character to get all unique levels
      CharacterVector char1, char2;
      
      if (Rf_inherits(col1, "factor")) {
        char1 = as<CharacterVector>(col1);
      } else {
        char1 = as<CharacterVector>(col1);
      }
      
      if (Rf_inherits(col2, "factor")) {
        char2 = as<CharacterVector>(col2);
      } else {
        char2 = as<CharacterVector>(col2);
      }
      
      // Create a mapping from string to integer code
      std::unordered_map<std::string, int> code_map;
      int next_code = 0;
      
      // First, add all unique values from data1
      for (int i = 0; i < n1; ++i) {
        if (char1[i] != NA_STRING) {
          std::string val = std::string(char1[i]);
          if (code_map.find(val) == code_map.end()) {
            code_map[val] = next_code++;
          }
        }
      }
      
      // Also add all unique values from data2 to ensure they have codes
      for (int j = 0; j < n2; ++j) {
        if (char2[j] != NA_STRING) {
          std::string val = std::string(char2[j]);
          if (code_map.find(val) == code_map.end()) {
            code_map[val] = next_code++;
          }
        }
      }
      
      // Convert both columns to integer codes
      std::vector<int> coded1(n1, NA_INTEGER);
      std::vector<int> coded2(n2, NA_INTEGER);
      
      for (int i = 0; i < n1; ++i) {
        if (char1[i] != NA_STRING) {
          std::string val = std::string(char1[i]);
          coded1[i] = code_map[val];
        }
      }
      
      for (int j = 0; j < n2; ++j) {
        if (char2[j] != NA_STRING) {
          std::string val = std::string(char2[j]);
          coded2[j] = code_map[val];
        }
      }
      
      cat_cols1.push_back(coded1);
      cat_cols2.push_back(coded2);
      cat_nlevels.push_back(next_code);
    }
  }
  
  int n_num = num_cols1.size();
  int n_cat = cat_cols1.size();
  
  NumericVector min_distances(n2);
  
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int j = 0; j < n2; ++j) {
    double best_similarity = -1.0;
    
    for (int i = 0; i < n1; ++i) {
      double sum_sij_wij = 0.0;
      double sum_wij = 0.0;
      
      // Process numeric columns
      for (int k = 0; k < n_num; ++k) {
        double val1 = num_cols1[k][i];
        double val2 = num_cols2[k][j];
        
        if (NumericVector::is_na(val1) || NumericVector::is_na(val2)) continue;
        
        double Rk = num_ranges[k];
        double sij = 0.0;
        if (Rk != 0.0) {
          sij = std::abs(val1 - val2) / Rk;
        }
        
        sum_sij_wij += sij;
        sum_wij += 1.0;
      }
      
      // Process categorical columns with integer codes
      for (int k = 0; k < n_cat; ++k) {
        int val1 = cat_cols1[k][i];
        int val2 = cat_cols2[k][j];
        
        if (val1 == NA_INTEGER || val2 == NA_INTEGER) continue;
        
        // For Gower distance, different categories = distance 1
        // This correctly handles categories from data2 not in data1
        // because they'll have different integer codes
        double sij = (val1 != val2) ? 1.0 : 0.0;
        
        sum_sij_wij += sij;
        sum_wij += 1.0;
      }
      
      if (sum_wij > 0.0) {
        double distance = sum_sij_wij / sum_wij;
        double similarity = 1.0 - distance;
        if (similarity > best_similarity) {
          best_similarity = similarity;
        }
      }
    }
    
    min_distances[j] = (best_similarity < 0.0) ? 0.0 : best_similarity;
  }
  
  return min_distances;
}