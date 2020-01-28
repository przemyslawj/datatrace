// Calculates place fields and two place cell related measures: Spatial information and Mutual information.
// The spatial information is tested for significance by shuffling the traces and comparing against the shuffled values
//
// The calculated mutual information is biased due to a finite and limited sampling of responses in the space.
// This bias is calculated using method from Panzeri and Treves, 1996.
//
#include <RcppArmadillo.h>
#include "mutual_info.h"
#include "trace_utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <cmath>
#include <float.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>


using namespace Rcpp;

class MfrModel {
public:
  NumericVector occupancyMap;
  NumericVector fr;
  double mfr;
  double min_trace;
  double sparsity;
  
  MfrModel(NumericVector& occupancyMap,
                    NumericVector& fr,
                    double mfr,
                    double min_trace,
                    double sparsity) {
    this->occupancyMap = occupancyMap;
    this->fr = fr;
    this->mfr = mfr;
    this->min_trace = min_trace;
    this->sparsity = sparsity;
  }
};

class SpatialInfoData {
public:
  NumericVector occupancyMap;
  NumericVector fr;
  double mfr = 0.0;
  double SI = 0.0;
  double MI = 0.0;
  double MI_bias = 0.0;
  double space_sampling_factor = 0.0;
  double sparsity = 0.0;

  SpatialInfoData() {
    occupancyMap = NumericVector(1, 0.0);
    fr = NumericVector(1, 0.0);
  }

  SpatialInfoData(NumericVector& occupancyMap,
                  NumericVector& fr,
                  double mfr,
                  double SI,
                  double MI,
                  double MI_bias,
                  double space_sampling_factor,
                  double sparsity) {
    this->occupancyMap = occupancyMap;
    this->fr = fr;
    this->mfr = mfr;
    this->SI = SI;
    this->MI = MI;
    this->MI_bias = MI_bias;
    this->space_sampling_factor = space_sampling_factor;
    this->sparsity = sparsity;
  }
};

MfrModel createMfrModel(IntegerVector& bin_xy,
                        int nstim,
                        NumericVector& trace,
                        double minOccupancy) {
  NumericVector totalActivityMap = NumericVector(nstim, 0.0);
  NumericVector occupancyMap = NumericVector(nstim, 0.0);
  NumericVector fr = NumericVector(nstim, 0.0);
  
  // Calculate occupancy and total activity maps
  double mfr = 0;
  for (int i = 0; i < trace.size(); ++i) {
    int xy = std::max(0, bin_xy[i] - 1);
    if (xy >= nstim) {
      stop("Stim value higher than the number of bins");
    }
    double trace_val = trace[i];
    occupancyMap[xy] += 1;
    totalActivityMap[xy] += trace_val;
    mfr += trace_val / trace.size();
  }  
  
  double sparsity = 0.0;
  for (int xy = 0; xy < nstim; ++xy) {
    if (occupancyMap[xy] >= minOccupancy) {
      fr[xy] = totalActivityMap[xy] / occupancyMap[xy];
      double p_s = ((double) occupancyMap[xy]) / trace.size();
      sparsity += p_s * std::pow(fr[xy], 2) / std::pow(mfr, 2);
    } else {
      fr[xy] = NAN;
    }
  }

  double min_trace = (double) *(std::min_element(trace.cbegin(), trace.cend()));
  return MfrModel(occupancyMap, fr, mfr, min_trace, sparsity);
}

// [[Rcpp::export]]
SEXP create_mfr_model(IntegerVector& bin_xy,
                      int nstim,
                      NumericVector& trace,
                      double minOccupancy) {
  MfrModel mfrModel = createMfrModel(bin_xy, nstim, trace, minOccupancy);
  List res;
  res["occupancyMap"] = mfrModel.occupancyMap;
  res["fr"] = mfrModel.fr;
  res["mfr"] = mfrModel.mfr;
  res["min_trace"] = mfrModel.min_trace;
  res["sparsity"] = mfrModel.sparsity;
  return res;
}

double calculateSI(MfrModel mfrModel, int minOccupancy, int N) {
  double SI = 0.0;
  
  // normalize the trace values to be positive and higher than 0.01
  double trace_offset = 0.0; 
  if (mfrModel.min_trace < 0.01) {
    trace_offset = mfrModel.min_trace - 0.01;
  }
  
  double mfr = mfrModel.mfr - trace_offset;
  for (int xy = 0; xy < mfrModel.occupancyMap.size(); ++xy) {
    if (mfrModel.occupancyMap[xy] >= minOccupancy) {
      Debug(" xy=" << xy << std::endl);
      double p_occupancy = (double) mfrModel.occupancyMap[xy] / N;
      double fr_xy = mfrModel.fr[xy] - trace_offset;
      double r_SI = p_occupancy * fr_xy / mfr * log2(fr_xy / mfr);
      Debug("partial SI=" << r_SI << std::endl);
      SI += r_SI;
    }
  }
  
  return SI;
}

SpatialInfoData calculateSpatialInformation(IntegerVector& bin_xy,
                                            int nstim,
                                            NumericVector& trace,
                                            NumericVector& binnedTrace,
                                            int nresponse,
                                            double minOccupancy) {

  const int N = trace.size();
  MfrModel mfrModel = createMfrModel(bin_xy, nstim, trace, minOccupancy);
  double SI = calculateSI(mfrModel, minOccupancy, N);
  
  BinnedResponseModel m = createResponseModel(binnedTrace, nresponse, bin_xy, nstim);
  MI_Data miData = modelMutualInfo(m);
  
  int occupiedBins = 0;
  for (int xy = 0; xy < nstim; ++xy) {
    if (mfrModel.occupancyMap[xy] >= minOccupancy) {
      ++occupiedBins;
    }
  }

  Debug("MFR=" << mfrModel.mfr);
  Debug(",N=" << N);
  Debug(", MI_bias=" << miData.MI_bias <<std::endl);
  double space_sampling_factor = ((double) occupiedBins) / nstim;
  return SpatialInfoData(mfrModel.occupancyMap, mfrModel.fr, mfrModel.mfr, 
                         SI, miData.MI, miData.MI_bias, space_sampling_factor, 
                         mfrModel.sparsity);
}

// [[Rcpp::export]]
SEXP calcPlaceField(IntegerVector& bin_xy,
                    int nstim,
                    NumericVector& trace,
                    NumericVector& binnedTrace,
                    IntegerVector& trialEnds,
                    int nshuffles,
                    int shuffleChunkLength,
                    double minOccupancy) {

  NumericVector shuffleSI = NumericVector(nshuffles, 0.0);
  NumericVector shuffleMI = NumericVector(nshuffles, 0.0);

  List result;
  result["field"] = NumericVector(nstim, 0.0);
  result["occupancy"] = NumericVector(nstim, 0.0);
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.size.50"] = 0;
  result["field.size.25"] = 0;
  result["shuffle.si"]= shuffleSI;
  result["shuffle.mi"]= shuffleMI;
  result["space.sampling.factor"] = 0.0;
  result["sparsity"] = 0.0;

  if (trace.size() == 0) {
    return(result);
  }
  int nresponse = (int) *(std::max_element(binnedTrace.cbegin(), binnedTrace.cend()));

  SpatialInfoData spatialInfoData = calculateSpatialInformation(bin_xy, nstim, 
                                                                trace, binnedTrace, 
                                                                nresponse, minOccupancy);
  NumericVector occupancyMap = spatialInfoData.occupancyMap;
  NumericVector fr = spatialInfoData.fr;

  double maxField = 0.0;
  for (int xy = 0; xy < nstim; ++xy) {
    if (occupancyMap[xy] >= minOccupancy && fr[xy] >= maxField) {
      maxField = fr[xy];
    }
  }
  
  // Calculate field size
  const double FIELD_BINARY_THRESH = 0.5;
  int nFieldBins50 = 0;
  int nFieldBins25 = 0;
  for (int xy = 0; xy < nstim; ++xy) {
    if (!std::isnan(fr[xy])) {
      if (fr[xy] >= FIELD_BINARY_THRESH * maxField) {
        ++nFieldBins50;
      }
      if (fr[xy] >= 0.25 * maxField) {
        ++nFieldBins25;
      }
    }
  }

  for (int i = 0; i < nshuffles; ++i) {
    //NumericVector shuffledTrace = chunkShuffle(trace, trialEnds, shuffleChunkLength);
    NumericVector shuffledTrace = randomShift(trace, trialEnds, shuffleChunkLength * 2);
    NumericVector shuffledBinnedTrace = randomShift(binnedTrace, trialEnds, shuffleChunkLength * 2);
    SpatialInfoData shuffleData = calculateSpatialInformation(bin_xy, nstim, 
                                                              shuffledTrace, shuffledBinnedTrace, 
                                                              nresponse, minOccupancy);
    shuffleSI[i] = shuffleData.SI;
    shuffleMI[i] = shuffleData.MI;
  }


  // Populate the result list
  result["field"] = fr;
  result["occupancy"] = occupancyMap;
  result["spatial.information"] = (float) spatialInfoData.SI;
  result["spatial.information.perspike"] = spatialInfoData.SI / spatialInfoData.mfr;
  result["mfr"] = spatialInfoData.mfr;
  result["field.size.50"] = ((double) nFieldBins50) / nstim * 100.0;
  result["field.size.25"] = ((double) nFieldBins25) / nstim * 100.0;
  result["shuffle.si"] = shuffleSI;
  result["mutual.info"] = spatialInfoData.MI;
  result["mutual.info.bias"] = spatialInfoData.MI_bias;
  result["shuffle.mi"] = shuffleMI;
  result["space.sampling.factor"] = spatialInfoData.space_sampling_factor;
  result["sparsity"] = spatialInfoData.sparsity;
  return(result);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
/*** R
xy =  c(1, 2, 3, 4, 5, 5, 4, 3)
binnedTrace=c(0, 0, 0, 1, 1, 0, 1, 0)
trace = binnedTrace / 2

pf = calcPlaceField(xy, 6, trace, binnedTrace + 1, c(7, length(trace)), 0, 2, 1)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace, traceQuantiles, 10,2), 4)
pf$spatial.information
pf$mutual.info
pf$mfr

calcPlaceField(xy, 3, trace, binnedTrace + 1, c(7, length(trace)), 0, 2, 1)
*/

