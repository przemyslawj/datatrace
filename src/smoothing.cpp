#include "smoothing.h"

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

/*** R
createGaussianKernel(0, 2)

kernel = createGaussianKernel(5, 2)
image(kernel)
*/
