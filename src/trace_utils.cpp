#include <Rcpp.h>
using namespace Rcpp;

#include "trace_utils.h"

NumericVector chunkShuffle(NumericVector& trace,
                           IntegerVector& trialEnds, 
                          int shuffleChunkLength) {
  
  // Shuffle the trace keeping the order within small chunks
  int nchunks = std::ceil(((double) trace.size() / shuffleChunkLength));
  std::vector<int> chunkShuffle(nchunks);
  for (int j = 0; j < nchunks; ++j) {
    chunkShuffle[j] = j;
  }
  
  NumericVector shuffledTrace(trace.size(), 0.0);
  
  // Shuffle chunks only within the same trial
  int trialStart = 0;
  for (int trial_i = 0; trial_i < trialEnds.size(); ++trial_i) {
    int trialEnd = std::min((int) trialEnds[trial_i] / shuffleChunkLength, nchunks);
    std::random_shuffle(chunkShuffle.begin() + trialStart, chunkShuffle.begin() + trialEnd);
    trialStart = trialEnd;
  }
  
  // Offset from the start to randomize chunk boundaries
  int offset = std::rand() % ((int) shuffleChunkLength/2);
  for (int chunk = 0; chunk < nchunks; ++chunk) {
    int src_chunk_start = chunk * shuffleChunkLength + offset;
    int target_chunk_start = chunkShuffle[chunk] * shuffleChunkLength + offset;
    for (int j = 0; j < shuffleChunkLength; ++j) {
      int src_index = (src_chunk_start + j) % trace.size();
      int target_index = (target_chunk_start + j) % trace.size();
      shuffledTrace[target_index] = trace[src_index];
    }
  }
  
  return shuffledTrace;
}

NumericVector randomShift(NumericVector& trace,
                          IntegerVector& trialEnds, 
                          int minShift) {
  NumericVector shiftedTrace(trace.size(), 0.0);
  
  int trialStart = 0;
  for (int trial_i = 0; trial_i < trialEnds.size(); ++trial_i) {
    trialEnds[trial_i] = std::min((int) trialEnds[trial_i], (int) trace.size());
    int ntrial = trialEnds[trial_i] - trialStart;
    if (ntrial <= 2 * minShift) {
      minShift = std::max(0, ntrial / 2 - 1);
    }
    int shift = std::rand() % (ntrial - 2 * minShift) + minShift;
    for (int i = 0; i < ntrial; ++i) {
      int within_trial_i = (i + shift) % ntrial;
      shiftedTrace[trialStart + i] = trace[trialStart + within_trial_i];
    }
    
    trialStart = trialEnds[trial_i];
  }
  
  return shiftedTrace;
}



/***
 * Requires DF to be ordered by trial_id, cell_id and timestamp, so the ajacent rows have
 * sequential timestamps.
 */
// [[Rcpp::export]]
LogicalVector isRunning(DataFrame& df,
                        double min_run_velocity,
                        double mean_run_velocity,
                        double window_dur_ms) {

  LogicalVector result(df.nrows());
  
  NumericVector velocity = df["velocity"];
  NumericVector pos_x = df["x"];
  NumericVector pos_y = df["y"];
  IntegerVector timestamp = df["timestamp"];
  CharacterVector trial_id = df["trial_id"];
  
  int i = 0;
  while (i < df.nrows()) {
    int j = i + 1;
    while(j < df.nrows() && 
          velocity[i] >= min_run_velocity && 
          trial_id[i] == trial_id[j] &&
          velocity[j] >= min_run_velocity) {
      j = j + 1;
    }
    j = j - 1;
    
    double distx = pos_x[j] - pos_x[i];
    double disty = pos_y[j] - pos_y[i];
    double dist = sqrt(distx * distx + disty * disty);
    double dur_ms = std::max(timestamp[j] - timestamp[i], 1);
    double vel = dist / dur_ms * 1000;
    bool is_running = (vel >= mean_run_velocity) && (dur_ms >= window_dur_ms);
    for (int k = i; k <= j; ++k) {
      result[k] = is_running;
    }
    i = j + 1;
  }
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
