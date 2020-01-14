#ifndef MUTUAL_INFO_H_
#define MUTUAL_INFO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cmath>

using namespace Rcpp;


struct MI_Data {
  double MI;
  double MI_bias;
  int totalResponseBins;
};


class BinnedResponseModel {
public:
  NumericMatrix prob_r_given_s;
  NumericVector prob_response;
  NumericVector prob_stim;
  int nstim;
  int nresponse;
  int N;
};

BinnedResponseModel createResponseModel(NumericVector& response, 
                                        int nresponseBins,
                                        IntegerVector& stimulus, 
                                        int nstim);
MI_Data modelMutualInfo(BinnedResponseModel& m);

#endif