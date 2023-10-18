#include <cstdlib>
#include <iostream>
#include <string>
#include <random>
#include <RcppArmadillo.h>


using namespace arma;
using namespace Rcpp;

// Function to get indices of elements equal to 1
uvec get_indices_of_ones(uvec x) {
  int n = x.size();
  uvec indices;
  for (int i = 0; i < n; i++) {
    if (x[i] == 1) {
      // Create a new uvector one element longer and copy the elements
      arma::uvec newIndices(indices.n_elem + 1);
      newIndices.head(indices.n_elem) = indices;
      // Append the new element to the end
      newIndices(indices.n_elem) = i;
      // Replace the original indices with the newIndices if needed
      indices = newIndices;
    }
  }
  return indices;
}


mat get_woodbury_inv(mat U, mat D) {
  mat I = eye(U.n_rows, U.n_rows);
  mat sigma_inv = I - U*inv(inv(D) + U.t() * U) * U.t();
  return sigma_inv;
}

