#ifndef TRACE_UTILS_H_
#define TRACE_UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

//#define USEDEBUG

#ifdef USEDEBUG
#define Debug(x) Rcout << x
#else
#define Debug(x)
#endif


NumericVector chunkShuffle(NumericVector& trace,
                           IntegerVector& trialEnds, 
                           int shuffleChunkLength);

NumericVector randomShift(NumericVector& trace,
                          IntegerVector& trialEnds, 
                          int minShift);
#endif