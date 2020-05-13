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

/**
 * Calculates MfrModel for the stimulus values and the corresponding firing signal.
 * 
 * Sparsity across the stimulus values, e.g. space, is calculated as:
 * sparsity = (\sum_x p_x * fr_x)^2 / (\sum_x (mfr)^2), where:
 * - x is a stimulus bin, e.g. spatial bin,
 * - p_x - probability of stimulus bin occupancy
 * - fr_x - firing rate for bin x
 * - mfr - mean firing rate
 */
MfrModel createMfrModel(IntegerVector& bin_xy,
                        int nstim,
                        NumericVector& trace,
                        int minOccupancy) {
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
                      int minOccupancy) {
  MfrModel mfrModel = createMfrModel(bin_xy, nstim, trace, minOccupancy);
  List res;
  res["occupancyMap"] = mfrModel.occupancyMap;
  res["fr"] = mfrModel.fr;
  res["mfr"] = mfrModel.mfr;
  res["min_trace"] = mfrModel.min_trace;
  res["sparsity"] = mfrModel.sparsity;
  return res;
}

/*
 * Spatial information per spike, calculated as: 
 * \sum_xy p_xy * fr_xy / mfr * log2(fr_xy / mfr), where:
 * - xy is a spatial bin
 * - p_xy - occupancy probability of that bin,
 * - fr_xy - firing rate or mean trace value in that bin
 * - mfr - mean firing rate or mean trace value
 */
double calculateSI(MfrModel mfrModel, int minOccupancy, int N) {
  double SI = 0.0;
  
  // normalize the firing rates to be positive and higher than 0.01
  // this is so the log values can be calculated
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


double getSpaceSamplingFactor(MfrModel& mfrModel, int nstim, int minOccupancy) {
  int occupiedBins = 0;
  for (int xy = 0; xy < nstim; ++xy) {
    if (mfrModel.occupancyMap[xy] >= minOccupancy) {
      ++occupiedBins;
    }
  }
  
  return ((double) occupiedBins) / nstim;
}


// [[Rcpp::export]]
SEXP calcPlaceField(IntegerVector& bin_xy,
                    int nstim,
                    NumericVector& trace,
                    NumericVector& binnedTrace,
                    int minOccupancy) {

  List result;
  result["field"] = NumericVector(nstim, 0.0);
  result["occupancy"] = NumericVector(nstim, 0.0);
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.size.50"] = 0;
  result["field.size.25"] = 0;
  result["space.sampling.factor"] = 0.0;
  result["sparsity"] = 0.0;

  if (trace.size() == 0) {
    return(result);
  }
  int nresponse = (int) *(std::max_element(binnedTrace.cbegin(), binnedTrace.cend()));

  MfrModel mfrModel = createMfrModel(bin_xy, nstim, trace, minOccupancy);
  double SI = calculateSI(mfrModel, minOccupancy, trace.size());
  double space_sampling_factor = getSpaceSamplingFactor(mfrModel, nstim, minOccupancy);
  Debug("MFR=" << mfrModel.mfr);
  NumericVector occupancyMap = mfrModel.occupancyMap;
  NumericVector fr = mfrModel.fr;

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

  BinnedResponseModel m = createResponseModel(binnedTrace, nresponse, bin_xy, nstim);
  MI_Data miData = modelMutualInfo(m);
  Debug(", MI_bias=" << miData.MI_bias <<std::endl);

  // Populate the result list
  result["field"] = fr;
  result["occupancy"] = occupancyMap;
  result["spatial.information"] = (float) SI;
  result["spatial.information.perspike"] = SI / mfrModel.mfr;
  result["mfr"] = mfrModel.mfr;
  result["field.size.50"] = ((double) nFieldBins50) / nstim * 100.0;
  result["field.size.25"] = ((double) nFieldBins25) / nstim * 100.0;
  result["mutual.info"] = miData.MI;
  result["mutual.info.bias"] = miData.MI_bias;
  result["space.sampling.factor"] = space_sampling_factor;
  result["sparsity"] = mfrModel.sparsity;
  return(result);
}


// [[Rcpp::export]]
SEXP placeFieldStatsForShuffled(IntegerVector& bin_xy,
                                int nstim,
                                NumericVector& trace,
                                NumericVector& binnedTrace,
                                IntegerVector& trialEnds,
                                int nshuffles,
                                int minShift,
                                int minOccupancy) {
  NumericVector shuffleSI = NumericVector(nshuffles, 0.0);
  NumericVector shuffleMI = NumericVector(nshuffles, 0.0);
  
  List result;
  result["shuffle.si"]= shuffleSI;
  result["shuffle.mi"]= shuffleMI;
  int nresponse = (int) *(std::max_element(binnedTrace.cbegin(), binnedTrace.cend()));
  for (int i = 0; i < nshuffles; ++i) {
    //NumericVector shuffledTrace = chunkShuffle(trace, trialEnds, shuffleChunkLength);
    NumericVector shuffledTrace = randomShift(trace, trialEnds, minShift);
    MfrModel mfrModel = createMfrModel(bin_xy, nstim, shuffledTrace, minOccupancy);
    double SI = calculateSI(mfrModel, minOccupancy, shuffledTrace.size());

    NumericVector shuffledBinnedTrace = randomShift(binnedTrace, trialEnds, minShift);
    BinnedResponseModel m = createResponseModel(shuffledBinnedTrace, nresponse, bin_xy, nstim);
    MI_Data miData = modelMutualInfo(m);

    shuffleSI[i] = SI;
    shuffleMI[i] = miData.MI;
  }
  
  result["shuffle.si"] = shuffleSI;
  result["shuffle.mi"] = shuffleMI;
  return result;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
/*** R
xy =        c(1, 2, 3, 4, 5, 5, 4, 3)
binnedTrace=c(0, 0, 0, 1, 0, 0, 1, 0)
trace = binnedTrace / 2

pf = calcPlaceField(xy, 6, trace, binnedTrace + 1, 1)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace, traceQuantiles, 10,2), 4)
pf$spatial.information
pf$mutual.info
pf$mfr

shuffle.pf = placeFieldStatsForShuffled(xy, 6, trace, binnedTrace + 1, c(4, length(trace)), 3, 1, 1)
shuffle.pf
*/

