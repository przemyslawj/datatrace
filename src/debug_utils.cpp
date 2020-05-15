#include "debug_utils.h"

using namespace Rcpp;


void printMat(arma::mat& m) {
  Debug("Matrix values:" << std::endl);
  for (int x = 0; x < m.n_rows; ++x) {
    for (int y = 0; y < m.n_cols; ++y) {
      Debug(m(x,y) << " ");  
    }
    Debug(std::endl);
  }
}

/*** R

*/
