// Calculates place fields and two place cell related measures: Spatial information and Mutual information.
// The spatial information is tested for significance by shuffling the traces and comparing against the shuffled values
//
// The calculated mutual information is biased due to a finite and limited sampling of responses in the space.
// This bias is calculated using method from Panzeri and Treves, 1996.
//
#include <RcppArmadillo.h>
#include "debug_utils.h"
#include "mutual_info.h"
#include "trace_utils.h"
#include "smoothing.h"
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
  arma::mat occupancyMap;
  arma::mat totalActivityMap;
  arma::mat fr;
  double mfr;
  double min_trace;
  double sparsity;
  double N;
  
  MfrModel(arma::mat& occupancyMap,
           arma::mat& totalActivityMap,
           arma::mat& fr,
           double mfr,
           double min_trace,
           double sparsity) {
    this->occupancyMap = occupancyMap;
    this->totalActivityMap = totalActivityMap;
    this->fr = fr;
    this->mfr = mfr;
    this->min_trace = min_trace;
    this->sparsity = sparsity;
    this->N = arma::accu(occupancyMap);
  }
};

arma::mat calcFrMap(arma::mat& occupancyMap, arma::mat& totalActivityMap, double minOccupancy) {
  arma::mat fr(occupancyMap.n_rows, occupancyMap.n_cols, arma::fill::zeros);
  for (int x = 0; x < occupancyMap.n_rows; ++x) {
    for (int y = 0; y < occupancyMap.n_cols; ++y) {
      if (occupancyMap(x,y) >= minOccupancy) {
        fr(x,y) = totalActivityMap(x,y) / occupancyMap(x,y);
      } else {
        fr(x,y) = NA_REAL;
      }
    }
  }
  return fr;
}

/**
 * Calculates MfrModel for the stimulus values and the corresponding firing signal.
 * 
 * Sparsity across the stimulus values, e.g. space, is calculated as:
 * sparsity = mfr^2 / (\sum_x (p_x * fr_x^2), where:
 * - x is a stimulus bin, e.g. spatial bin,
 * - p_x - probability of stimulus bin occupancy
 * - fr_x - firing rate for bin x
 * - mfr - mean firing rate
 * Definition as in Boccara, Science 2019
 */
MfrModel createMfrModel(IntegerVector& bin_x,
                        IntegerVector& bin_y,
                        int nbins_x,
                        int nbins_y,
                        NumericVector& trace,
                        double minOccupancy) {
  arma::mat totalActivityMap(nbins_x, nbins_y, arma::fill::zeros);
  arma::mat occupancyMap(nbins_x, nbins_y, arma::fill::zeros);
  
  // Calculate occupancy and total activity maps
  double mfr = 0;
  for (int i = 0; i < trace.size(); ++i) {
    int x = std::max(0, bin_x[i] - 1);
    int y = std::max(0, bin_y[i] - 1);
    if (x >= nbins_x || y >= nbins_y) {
      stop("Stim value higher than the number of bins");
    }
    double trace_val = trace[i];
    occupancyMap(x,y) += 1;
    totalActivityMap(x,y) += trace_val;
    mfr += trace_val / trace.size();
  }  
  
  arma::mat fr = calcFrMap(occupancyMap, totalActivityMap, minOccupancy);
  double sparsity = 0.0;
  for (int x = 0; x < occupancyMap.n_rows; ++x) {
    for (int y = 0; y < occupancyMap.n_cols; ++y) {
      if (occupancyMap(x,y) >= minOccupancy && !std::isnan(fr(x,y))) {
        double p_s = (double) occupancyMap(x,y) / trace.size();
        sparsity += std::pow(mfr, 2) / (p_s * std::pow(fr(x,y), 2));
      }
    }
  }

  double min_trace = (double) *(std::min_element(trace.cbegin(), trace.cend()));
  return MfrModel(occupancyMap, totalActivityMap, fr, mfr, min_trace, sparsity);
}

// [[Rcpp::export]]
SEXP create_mfr_model(IntegerVector& bin_x,
                      IntegerVector& bin_y,
                      int nbins_x,
                      int nbins_y,
                      NumericVector& trace,
                      double minOccupancy) {
  MfrModel mfrModel = createMfrModel(bin_x, bin_y, nbins_x, nbins_y, trace, minOccupancy);
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
double calculateSI(MfrModel mfrModel, double N) {
  double SI = 0.0;
  
  // normalize the firing rates to be positive and higher than 0.01
  // this is so the log values can be calculated
  double trace_offset = 0.0; 
  if (mfrModel.min_trace < 0.01) {
    trace_offset = mfrModel.min_trace - 0.01;
  }
  
  double mfr = mfrModel.mfr - trace_offset;
  for (int x = 0; x < mfrModel.occupancyMap.n_rows; ++x) {
    for (int y = 0; y < mfrModel.occupancyMap.n_cols; ++y) {
      if (!std::isnan(mfrModel.fr(x,y)) && mfrModel.fr(x,y) > 0) {
        Debug(" x=" << x << ", y=" << y << std::endl);
        double p_occupancy = (double) mfrModel.occupancyMap(x,y) / N;
        double fr_xy = mfrModel.fr(x,y) - trace_offset;
        double r_SI = p_occupancy * fr_xy / mfr * log2(fr_xy / mfr);
        Debug("partial SI=" << r_SI << std::endl);
        SI += r_SI;
      }
    }
  }
  
  return SI;
}

double getSpaceSamplingFactor(MfrModel& mfrModel) {
  int occupiedBins = 0;
  for (int x = 0; x < mfrModel.occupancyMap.n_rows; ++x) {
    for (int y = 0; y < mfrModel.occupancyMap.n_cols; ++y) {
      if (!std::isnan(mfrModel.fr(x,y))) {
        ++occupiedBins;
      }
    }
  }
  
  return ((double) occupiedBins) / (mfrModel.occupancyMap.n_rows * mfrModel.occupancyMap.n_cols);
}

// Smooths the fr map and totalActivity map convolving the kernel.
MfrModel smoothMfr(MfrModel& mfrModel, arma::mat& kernel, double minOccupancy) {
  arma::mat smoothOccupancy = conv2(mfrModel.occupancyMap, kernel, "same");
  arma::mat smoothTotalActivity = conv2(mfrModel.totalActivityMap, kernel, "same");

  double smoothedMinOccupancy = findSmoothMinOccupancy(mfrModel.occupancyMap, smoothOccupancy, minOccupancy);
  Debug("Smoothed minoccupancy=" << smoothedMinOccupancy << std::endl);

  // Calculate smooth FR map but using original occupancy values for thresholdin min occupancy
  arma::mat smoothFr(smoothOccupancy.n_rows, smoothOccupancy.n_cols, arma::fill::zeros);
  for (int x = 0; x < smoothOccupancy.n_rows; ++x) {
    for (int y = 0; y < smoothOccupancy.n_cols; ++y) {
      if (smoothOccupancy(x,y) >= smoothedMinOccupancy) {
        smoothFr(x,y) = smoothTotalActivity(x,y) / smoothOccupancy(x,y);
      } else {
        smoothFr(x,y) = NA_REAL;
      }
    }
  }
  
  MfrModel smoothedModel = MfrModel(smoothOccupancy,
                                    smoothTotalActivity,
                                    smoothFr,
                                    mfrModel.mfr,
                                    mfrModel.min_trace,
                                    mfrModel.sparsity);
  return smoothedModel;
}

// Set kernel size to 0 not to apply gaussian smoothing.
// [[Rcpp::export]]
SEXP calcPlaceField(IntegerVector& bin_x,
                    IntegerVector& bin_y,
                    int nbins_x,
                    int nbins_y,
                    NumericVector& trace,
                    NumericVector& binnedTrace,
                    int minOccupancy,
                    int kernelSize,
                    double gaussianVar) {

  List result;
  result["field"] = arma::mat(nbins_x, nbins_y, arma::fill::zeros);
  result["occupancy"] = arma::mat(nbins_x, nbins_y, arma::fill::zeros);
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

  MfrModel mfrModel = createMfrModel(bin_x, bin_y, nbins_x, nbins_y, trace, minOccupancy);
  Debug("Model before smoothing");
  #ifdef USE_DEBUG
  printMat(mfrModel.fr);
  #endif
  arma::mat kernel = createGaussianKernel(kernelSize, gaussianVar);
  if (kernelSize > 0) {
    mfrModel = smoothMfr(mfrModel, kernel, minOccupancy);
  }
  Debug("Model after smoothing");
  #ifdef USE_DEBUG
  printMat(mfrModel.fr);
  #endif
  double SI = calculateSI(mfrModel, mfrModel.N);
  double space_sampling_factor = getSpaceSamplingFactor(mfrModel);
  Debug("MFR=" << mfrModel.mfr);

  double maxField = 0.0;
  for (int x = 0; x < nbins_x; ++x) {
    for (int y = 0; y < nbins_y; ++y) {
      if (mfrModel.fr(x,y) >= maxField) {
        maxField = mfrModel.fr(x,y);
      }
    }
  }
  
  // Calculate field size
  const double FIELD_BINARY_THRESH = 0.5;
  int nFieldBins50 = 0;
  int nFieldBins25 = 0;
  for (int x = 0; x < nbins_x; ++x) {
    for (int y = 0; y < nbins_y; ++y) {
      if (!std::isnan(mfrModel.fr(x,y))) {
        if (mfrModel.fr(x,y) >= FIELD_BINARY_THRESH * maxField) {
          ++nFieldBins50;
        }
        if (mfrModel.fr(x,y) >= 0.25 * maxField) {
          ++nFieldBins25;
        }
      }
    }
  }

  BinnedResponseModel2D m = create2DResponseModel(binnedTrace, nresponse, bin_x, bin_y, nbins_x, nbins_y, minOccupancy);
  if (kernelSize > 0) {
    m = smooth2DResponseModel(m, kernel, minOccupancy);
  }
  MI_Data miData = modelMutualInfo(m);
  Debug(", MI_bias=" << miData.MI_bias <<std::endl);
  result["mutual.info"] = miData.MI;
  result["mutual.info.bias"] = miData.MI_bias;

  // Populate the result list
  result["field"] = mfrModel.fr;
  result["occupancy"] = mfrModel.occupancyMap;
  result["spatial.information"] = (float) SI;
  result["spatial.information.perspike"] = SI / mfrModel.mfr;
  result["mfr"] = mfrModel.mfr;
  result["field.size.50"] = ((double) nFieldBins50) / (nbins_x * nbins_y) * 100.0;
  result["field.size.25"] = ((double) nFieldBins25) / (nbins_x * nbins_y) * 100.0;
  result["space.sampling.factor"] = space_sampling_factor;
  result["sparsity"] = mfrModel.sparsity;
  return(result);
}

// [[Rcpp::export]]
SEXP placeFieldStatsForShuffled(IntegerVector& bin_x,
                                IntegerVector& bin_y,
                                int nbins_x,
                                int nbins_y,
                                NumericVector& trace,
                                NumericVector& binnedTrace,
                                IntegerVector& trialEnds,
                                int nshuffles,
                                int minShift,
                                double minOccupancy,
                                int kernelSize,
                                double gaussianVar) {
  NumericVector shuffleSI = NumericVector(nshuffles, 0.0);
  NumericVector shuffleMI = NumericVector(nshuffles, 0.0);
  arma::mat kernel = createGaussianKernel(kernelSize, gaussianVar);
  
  List result;
  result["shuffle.si"] = shuffleSI;
  result["shuffle.mi"] = shuffleMI;
  int nresponse = (int) *(std::max_element(binnedTrace.cbegin(), binnedTrace.cend()));
  for (int i = 0; i < nshuffles; ++i) {
    //NumericVector shuffledTrace = chunkShuffle(trace, trialEnds, shuffleChunkLength);
    NumericVector shuffledTrace = randomShift(trace, trialEnds, minShift);
    MfrModel mfrModel = createMfrModel(bin_x, bin_y, nbins_x, nbins_y, shuffledTrace, minOccupancy);
    if (kernelSize > 0) {
      Debug("FR before smoothing");
      #ifdef USE_DEBUG
      printMat(mfrModel.fr);
      #endif
      mfrModel = smoothMfr(mfrModel, kernel, minOccupancy);
      Debug("FR after smoothing");
      #ifdef USE_DEBUG
      printMat(mfrModel.fr);
      #endif
    }
    shuffleSI[i] = calculateSI(mfrModel, mfrModel.N);

    NumericVector shuffledBinnedTrace = randomShift(binnedTrace, trialEnds, minShift);
    BinnedResponseModel2D m = create2DResponseModel(shuffledBinnedTrace, nresponse, bin_x, bin_y, nbins_x, nbins_y, minOccupancy);
    if (kernelSize > 0) {
      m = smooth2DResponseModel(m, kernel, minOccupancy);
    }
    MI_Data miData = modelMutualInfo(m);
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
x =          c(1, 1, 1, 1, 1, 1, 1, 1)
y =          c(1, 2, 3, 4, 5, 5, 4, 3)
binnedTrace=c(0, 0, 0, 1, 0, 0, 1, 0)
trace = binnedTrace / 2

pf = calcPlaceField(x, y, 2, 5, trace, binnedTrace + 1, 1, 0, 0)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace, traceQuantiles, 10,2), 4)
pf$spatial.information
pf$mutual.info
pf$mfr

pf.smooth = calcPlaceField(x, y, 2, 7, trace, binnedTrace + 1, 1, 3, 0.5)
pf.smooth$field

shuffle.pf = placeFieldStatsForShuffled(bin_x=x, bin_y=y, nbins_x=2, nbins_y=5, 
                                        trace=trace, binnedTrace=binnedTrace + 1, 
                                        trialEnds=c(4, length(trace)), 
                                        nshuffles=5, minShift=1, minOccupancy=1, 
                                        kernelSize=3, gaussianVar=0.5)
shuffle.pf
*/

