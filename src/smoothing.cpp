#include "smoothing.h"
#include "cmath"

// [[Rcpp::export]]
arma::mat createGaussianKernel(const int kernelSize, const double var) {
  arma::mat kernel(kernelSize, kernelSize);
  for (int i = 0; i < kernelSize; ++i) {
    for (int j = 0; j < kernelSize; ++j) {
      double x = i + 0.5 - (double) kernelSize / 2;
      double y = j + 0.5 - (double) kernelSize / 2;
      kernel(i,j) = 1 / (2 * M_PI * var) * std::exp(-(x*x + y*y) / (2 * var)); 
    }
  }
  return kernel;
}

// [[Rcpp::export]]
arma::mat convolve(arma::mat M, arma::mat k) {
  return arma::conv2(M, k, "same");
}

double findSmoothMinOccupancy(arma::mat& occupancyMap, 
                              arma::mat& smoothOccupancyMap, 
                              double minOccupancy) {
  double minOccupancyVal = DBL_MAX;
  double smoothedMinOccupancy = DBL_MAX;
  for (int x = 0; x < smoothOccupancyMap.n_rows; ++x) {
    for (int y = 0; y < smoothOccupancyMap.n_cols; ++y) {
      int occupancyVal = occupancyMap(x,y);
      if (occupancyVal >= minOccupancy && occupancyVal <= minOccupancyVal) {
        if (smoothOccupancyMap(x, y) < smoothedMinOccupancy) {
          minOccupancyVal = occupancyVal;
          smoothedMinOccupancy = smoothOccupancyMap(x, y); 
        }
      }
    }
  }

  return smoothedMinOccupancy;
}

/*** R
createGaussianKernel(0, 2)

kernel = createGaussianKernel(3, 1)

M = matrix(0, nrow=5,ncol=5)
M[3,3] = 1
M.conv = convolve(M, kernel)
image(M.conv)
*/
