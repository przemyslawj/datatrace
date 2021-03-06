// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bayesmax
SEXP bayesmax(NumericVector& prior, NumericVector& likelihood, NumericVector& pv);
RcppExport SEXP _datatrace_bayesmax(SEXP priorSEXP, SEXP likelihoodSEXP, SEXP pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type likelihood(likelihoodSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type pv(pvSEXP);
    rcpp_result_gen = Rcpp::wrap(bayesmax(prior, likelihood, pv));
    return rcpp_result_gen;
END_RCPP
}
// poisson_prob
double poisson_prob(double lambda, double k);
RcppExport SEXP _datatrace_poisson_prob(SEXP lambdaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(poisson_prob(lambda, k));
    return rcpp_result_gen;
END_RCPP
}
// gauss_prob
double gauss_prob(double m, double sig, double p);
RcppExport SEXP _datatrace_gauss_prob(SEXP mSEXP, SEXP sigSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_prob(m, sig, p));
    return rcpp_result_gen;
END_RCPP
}
// bayesmax_mfr
SEXP bayesmax_mfr(NumericVector& prior, NumericVector& mfr, NumericVector& pv);
RcppExport SEXP _datatrace_bayesmax_mfr(SEXP priorSEXP, SEXP mfrSEXP, SEXP pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mfr(mfrSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type pv(pvSEXP);
    rcpp_result_gen = Rcpp::wrap(bayesmax_mfr(prior, mfr, pv));
    return rcpp_result_gen;
END_RCPP
}
// mutual_info
SEXP mutual_info(NumericVector& response, int nresponseBins, IntegerVector& stimulus, int nstim, int minStimOccurrence);
RcppExport SEXP _datatrace_mutual_info(SEXP responseSEXP, SEXP nresponseBinsSEXP, SEXP stimulusSEXP, SEXP nstimSEXP, SEXP minStimOccurrenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< int >::type nresponseBins(nresponseBinsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type stimulus(stimulusSEXP);
    Rcpp::traits::input_parameter< int >::type nstim(nstimSEXP);
    Rcpp::traits::input_parameter< int >::type minStimOccurrence(minStimOccurrenceSEXP);
    rcpp_result_gen = Rcpp::wrap(mutual_info(response, nresponseBins, stimulus, nstim, minStimOccurrence));
    return rcpp_result_gen;
END_RCPP
}
// mutual_info2D
SEXP mutual_info2D(NumericVector& response, int nresponseBins, IntegerVector& stimulus_x, IntegerVector& stimulus_y, int nstim_x, int nstim_y, int minStimOccurrence, int kernelSize, double gaussianVar);
RcppExport SEXP _datatrace_mutual_info2D(SEXP responseSEXP, SEXP nresponseBinsSEXP, SEXP stimulus_xSEXP, SEXP stimulus_ySEXP, SEXP nstim_xSEXP, SEXP nstim_ySEXP, SEXP minStimOccurrenceSEXP, SEXP kernelSizeSEXP, SEXP gaussianVarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< int >::type nresponseBins(nresponseBinsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type stimulus_x(stimulus_xSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type stimulus_y(stimulus_ySEXP);
    Rcpp::traits::input_parameter< int >::type nstim_x(nstim_xSEXP);
    Rcpp::traits::input_parameter< int >::type nstim_y(nstim_ySEXP);
    Rcpp::traits::input_parameter< int >::type minStimOccurrence(minStimOccurrenceSEXP);
    Rcpp::traits::input_parameter< int >::type kernelSize(kernelSizeSEXP);
    Rcpp::traits::input_parameter< double >::type gaussianVar(gaussianVarSEXP);
    rcpp_result_gen = Rcpp::wrap(mutual_info2D(response, nresponseBins, stimulus_x, stimulus_y, nstim_x, nstim_y, minStimOccurrence, kernelSize, gaussianVar));
    return rcpp_result_gen;
END_RCPP
}
// mutual_info_with_shuffles
SEXP mutual_info_with_shuffles(NumericVector& response, int nresponseBins, IntegerVector& stimulus, int nstim, IntegerVector& trialEnds, int nshuffles, int shuffleChunkLength, int minStimOccurrence);
RcppExport SEXP _datatrace_mutual_info_with_shuffles(SEXP responseSEXP, SEXP nresponseBinsSEXP, SEXP stimulusSEXP, SEXP nstimSEXP, SEXP trialEndsSEXP, SEXP nshufflesSEXP, SEXP shuffleChunkLengthSEXP, SEXP minStimOccurrenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< int >::type nresponseBins(nresponseBinsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type stimulus(stimulusSEXP);
    Rcpp::traits::input_parameter< int >::type nstim(nstimSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type trialEnds(trialEndsSEXP);
    Rcpp::traits::input_parameter< int >::type nshuffles(nshufflesSEXP);
    Rcpp::traits::input_parameter< int >::type shuffleChunkLength(shuffleChunkLengthSEXP);
    Rcpp::traits::input_parameter< int >::type minStimOccurrence(minStimOccurrenceSEXP);
    rcpp_result_gen = Rcpp::wrap(mutual_info_with_shuffles(response, nresponseBins, stimulus, nstim, trialEnds, nshuffles, shuffleChunkLength, minStimOccurrence));
    return rcpp_result_gen;
END_RCPP
}
// create_mfr_model
SEXP create_mfr_model(IntegerVector& bin_x, IntegerVector& bin_y, int nbins_x, int nbins_y, NumericVector& trace, double minOccupancy);
RcppExport SEXP _datatrace_create_mfr_model(SEXP bin_xSEXP, SEXP bin_ySEXP, SEXP nbins_xSEXP, SEXP nbins_ySEXP, SEXP traceSEXP, SEXP minOccupancySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type bin_x(bin_xSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type bin_y(bin_ySEXP);
    Rcpp::traits::input_parameter< int >::type nbins_x(nbins_xSEXP);
    Rcpp::traits::input_parameter< int >::type nbins_y(nbins_ySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< double >::type minOccupancy(minOccupancySEXP);
    rcpp_result_gen = Rcpp::wrap(create_mfr_model(bin_x, bin_y, nbins_x, nbins_y, trace, minOccupancy));
    return rcpp_result_gen;
END_RCPP
}
// calcPlaceField
SEXP calcPlaceField(IntegerVector& bin_x, IntegerVector& bin_y, int nbins_x, int nbins_y, NumericVector& trace, NumericVector& binnedTrace, int minOccupancy, int kernelSize, double gaussianVar);
RcppExport SEXP _datatrace_calcPlaceField(SEXP bin_xSEXP, SEXP bin_ySEXP, SEXP nbins_xSEXP, SEXP nbins_ySEXP, SEXP traceSEXP, SEXP binnedTraceSEXP, SEXP minOccupancySEXP, SEXP kernelSizeSEXP, SEXP gaussianVarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type bin_x(bin_xSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type bin_y(bin_ySEXP);
    Rcpp::traits::input_parameter< int >::type nbins_x(nbins_xSEXP);
    Rcpp::traits::input_parameter< int >::type nbins_y(nbins_ySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type binnedTrace(binnedTraceSEXP);
    Rcpp::traits::input_parameter< int >::type minOccupancy(minOccupancySEXP);
    Rcpp::traits::input_parameter< int >::type kernelSize(kernelSizeSEXP);
    Rcpp::traits::input_parameter< double >::type gaussianVar(gaussianVarSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPlaceField(bin_x, bin_y, nbins_x, nbins_y, trace, binnedTrace, minOccupancy, kernelSize, gaussianVar));
    return rcpp_result_gen;
END_RCPP
}
// placeFieldStatsForShuffled
SEXP placeFieldStatsForShuffled(IntegerVector& bin_x, IntegerVector& bin_y, int nbins_x, int nbins_y, NumericVector& trace, NumericVector& binnedTrace, IntegerVector& trialEnds, int nshuffles, int minShift, double minOccupancy, int kernelSize, double gaussianVar);
RcppExport SEXP _datatrace_placeFieldStatsForShuffled(SEXP bin_xSEXP, SEXP bin_ySEXP, SEXP nbins_xSEXP, SEXP nbins_ySEXP, SEXP traceSEXP, SEXP binnedTraceSEXP, SEXP trialEndsSEXP, SEXP nshufflesSEXP, SEXP minShiftSEXP, SEXP minOccupancySEXP, SEXP kernelSizeSEXP, SEXP gaussianVarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type bin_x(bin_xSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type bin_y(bin_ySEXP);
    Rcpp::traits::input_parameter< int >::type nbins_x(nbins_xSEXP);
    Rcpp::traits::input_parameter< int >::type nbins_y(nbins_ySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type binnedTrace(binnedTraceSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type trialEnds(trialEndsSEXP);
    Rcpp::traits::input_parameter< int >::type nshuffles(nshufflesSEXP);
    Rcpp::traits::input_parameter< int >::type minShift(minShiftSEXP);
    Rcpp::traits::input_parameter< double >::type minOccupancy(minOccupancySEXP);
    Rcpp::traits::input_parameter< int >::type kernelSize(kernelSizeSEXP);
    Rcpp::traits::input_parameter< double >::type gaussianVar(gaussianVarSEXP);
    rcpp_result_gen = Rcpp::wrap(placeFieldStatsForShuffled(bin_x, bin_y, nbins_x, nbins_y, trace, binnedTrace, trialEnds, nshuffles, minShift, minOccupancy, kernelSize, gaussianVar));
    return rcpp_result_gen;
END_RCPP
}
// createGaussianKernel
arma::mat createGaussianKernel(const int kernelSize, const double var);
RcppExport SEXP _datatrace_createGaussianKernel(SEXP kernelSizeSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type kernelSize(kernelSizeSEXP);
    Rcpp::traits::input_parameter< const double >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(createGaussianKernel(kernelSize, var));
    return rcpp_result_gen;
END_RCPP
}
// convolve
arma::mat convolve(arma::mat M, arma::mat k);
RcppExport SEXP _datatrace_convolve(SEXP MSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve(M, k));
    return rcpp_result_gen;
END_RCPP
}
// isRunning
LogicalVector isRunning(DataFrame& df, double min_run_velocity, double mean_run_velocity, double window_dur_ms);
RcppExport SEXP _datatrace_isRunning(SEXP dfSEXP, SEXP min_run_velocitySEXP, SEXP mean_run_velocitySEXP, SEXP window_dur_msSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type min_run_velocity(min_run_velocitySEXP);
    Rcpp::traits::input_parameter< double >::type mean_run_velocity(mean_run_velocitySEXP);
    Rcpp::traits::input_parameter< double >::type window_dur_ms(window_dur_msSEXP);
    rcpp_result_gen = Rcpp::wrap(isRunning(df, min_run_velocity, mean_run_velocity, window_dur_ms));
    return rcpp_result_gen;
END_RCPP
}
// chunkShuffle
NumericVector chunkShuffle(NumericVector& trace, IntegerVector& trialEnds, int shuffleChunkLength);
RcppExport SEXP _datatrace_chunkShuffle(SEXP traceSEXP, SEXP trialEndsSEXP, SEXP shuffleChunkLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type trialEnds(trialEndsSEXP);
    Rcpp::traits::input_parameter< int >::type shuffleChunkLength(shuffleChunkLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(chunkShuffle(trace, trialEnds, shuffleChunkLength));
    return rcpp_result_gen;
END_RCPP
}
// randomShift
NumericVector randomShift(NumericVector& trace, IntegerVector& trialEnds, int minShift);
RcppExport SEXP _datatrace_randomShift(SEXP traceSEXP, SEXP trialEndsSEXP, SEXP minShiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type trialEnds(trialEndsSEXP);
    Rcpp::traits::input_parameter< int >::type minShift(minShiftSEXP);
    rcpp_result_gen = Rcpp::wrap(randomShift(trace, trialEnds, minShift));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_datatrace_bayesmax", (DL_FUNC) &_datatrace_bayesmax, 3},
    {"_datatrace_poisson_prob", (DL_FUNC) &_datatrace_poisson_prob, 2},
    {"_datatrace_gauss_prob", (DL_FUNC) &_datatrace_gauss_prob, 3},
    {"_datatrace_bayesmax_mfr", (DL_FUNC) &_datatrace_bayesmax_mfr, 3},
    {"_datatrace_mutual_info", (DL_FUNC) &_datatrace_mutual_info, 5},
    {"_datatrace_mutual_info2D", (DL_FUNC) &_datatrace_mutual_info2D, 9},
    {"_datatrace_mutual_info_with_shuffles", (DL_FUNC) &_datatrace_mutual_info_with_shuffles, 8},
    {"_datatrace_create_mfr_model", (DL_FUNC) &_datatrace_create_mfr_model, 6},
    {"_datatrace_calcPlaceField", (DL_FUNC) &_datatrace_calcPlaceField, 9},
    {"_datatrace_placeFieldStatsForShuffled", (DL_FUNC) &_datatrace_placeFieldStatsForShuffled, 12},
    {"_datatrace_createGaussianKernel", (DL_FUNC) &_datatrace_createGaussianKernel, 2},
    {"_datatrace_convolve", (DL_FUNC) &_datatrace_convolve, 2},
    {"_datatrace_isRunning", (DL_FUNC) &_datatrace_isRunning, 4},
    {"_datatrace_chunkShuffle", (DL_FUNC) &_datatrace_chunkShuffle, 3},
    {"_datatrace_randomShift", (DL_FUNC) &_datatrace_randomShift, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_datatrace(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
