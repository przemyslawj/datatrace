#ifndef SMOOTHING_INFO_H_
#define SMOOTHING_INFO_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cmath>

arma::mat createGaussianKernel(const int kernelSize, const double var);

#endif