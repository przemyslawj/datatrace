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
  arma::mat prob_r_given_s;
  NumericVector prob_response;
  NumericVector prob_stim;
  int nresponse;
  int N;
};

class BinnedResponseModel2D : public BinnedResponseModel {
public:
  // Cube  (response x stim_x x stim_y) of probabilities
  arma::cube prob_r_given_xy;
  arma::cube total_r_given_xy;
  arma::mat count_stim_xy;
};

BinnedResponseModel2D create2DResponseModel(NumericVector& response, 
                                            int nresponseBins,
                                            IntegerVector& stimulus_x, 
                                            IntegerVector& stimulus_y, 
                                            int nstim_x,
                                            int nstim_y,
                                            int minOccurrence);

MI_Data modelMutualInfo(BinnedResponseModel2D& m);

BinnedResponseModel2D smooth2DResponseModel(BinnedResponseModel2D& m,
                                            arma::mat& kernel);

#endif