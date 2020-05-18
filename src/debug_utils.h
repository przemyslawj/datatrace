#ifndef DEBUG_UTILS_H_
#define DEBUG_UTILS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

#define USEDEBUG

//#ifdef USEDEBUG
#define Debug(x) Rcout << x
#else
#define Debug(x)
#endif

void printMat(arma::mat& m);

#endif
