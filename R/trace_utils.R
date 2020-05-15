xybins = 20

zscore = function(trace) {
  x = trace - mean(trace)
  sd_est = IQR(x) / 1.349
  return(x / sd_est)
}

sem = function(x) sqrt( var(x, na.rm=TRUE) / length(x))

read.data.trace = function(caimg_result_dir, filter_exp_title = NA) {
  print(paste('Processing dir: ', caimg_result_dir))
  cellmapping_file = file.path(caimg_result_dir, 'cell_mapping.csv')
  cellmapping.df = read.csv(cellmapping_file)
  traces_file = file.path(caimg_result_dir, 'traces_and_positions.csv')
  data = fread(traces_file)

  if (!is.na(filter_exp_title)) {
    data = data[exp_title == filter_exp_title,]
  }

  if (nrow(data) == 0) {
    return(data)
  }

  data.traces = melt.traces(data)
  data.traces = left_join(data.traces, cellmapping.df, by=c('cell'='cell_no'))
  dir_parts = str_split(caimg_result_dir, '/')
  animal = dir_parts[[1]][length(dir_parts[[1]]) - 2]
  if (nrow(data.traces) > 0) {
    data.traces$animal = animal
  }
  data.traces = data.table(data.traces)
  data.traces = data.traces[smooth_trans_x >= 0 & smooth_trans_y >= 0, ]

  setnames(data.traces, c('smooth_trans_x', 'smooth_trans_y', 'smooth_heading_angle'), c('x', 'y', 'angle'))
  return(data.traces)
}

.calc.event.vec = function(deconv_trace, deconv.threshold=0.1, minpeakdistance=1) {
  deconv.thr.val = deconv.threshold * max(deconv_trace)
  peak.mat = findpeaks(deconv_trace,
                       minpeakheight=deconv.thr.val,
                       minpeakdistance=minpeakdistance,
                       zero='+')
  res = rep(FALSE, length(deconv_trace))
  if (!is.null(peak.mat)) {
    res[peak.mat[,2]] = TRUE
  }
  res
}

detect.events = function(data.traces, deconv.threshold=0.1, minpeakdistance=1) {
  data.traces[, is.event := .calc.event.vec(.SD$deconv_trace, deconv.threshold, minpeakdistance) , by=c('animal', 'date', 'cell_id')]
  return(data.traces)
}

overlapping.timebin.traces = function(data.traces, timebin.dur.msec=500) {
  overlap.dur.msec = timebin.dur.msec / 2
  t1 = timebin.traces(data.traces, timebin.dur.msec)
  t1[, time_bin := time_bin * 2]

  t2 = timebin.traces(data.traces[, timestamp := timestamp + overlap.dur.msec], timebin.dur.msec)
  t2[, time_bin := time_bin * 2 + 1]

  merged.t = as.data.table(rbind(t1, t2))
  setorder(merged.t, time_bin, cell_id)
  return(merged.t)
}

timebin.traces = function(data.traces, timebin.dur.msec=200) {
  if (nrow(data.traces) == 0) {
    return(data.frame())
  }

  max.timestamp = max(data.traces$timestamp)

  data.traces[, abs_timestamp := (trial-1)*(max.timestamp+timebin.dur.msec) + timestamp]
  data.traces[, time_bin := floor(abs_timestamp/timebin.dur.msec) %>% as.integer]
  timebinned.traces = data.traces[,
                                  lapply(.SD, mean),
                                  by=.(animal, date, trial_id, trial, exp_title, cell_id, time_bin),
                                  .SDcols = !c('is.event')]
  setorder(timebinned.traces, time_bin, cell_id)
  nevents.traces = data.traces[, .(nevents=sum(is.event)),
                                 by=.(animal, date, trial_id, trial, exp_title, cell_id, time_bin)]
  setorder(nevents.traces, time_bin, cell_id)
  timebinned.traces$nevents = nevents.traces$nevents
  timebinned.traces
}

# Bins data.traces[[stim.var]] values, by dividing the into bins of bin.width.
# Bins are indexed from 1.
stimbin.traces = function(data.traces, stim.var, nbins, max.width=100) {
  if (nrow(data.traces) == 0) {
    return(data.traces)
  }
  stim.var = enquo(stim.var)
  bin.width = max.width / nbins

  bin.var.name = paste0('bin.', quo_name(stim.var))
  data.traces %>%
    dplyr::mutate(!!bin.var.name := as.integer(ceiling(pmin(!!stim.var, max.width) / bin.width)))
}


to_1dim = function(x, y, nbins.y) {
  # Vals from 1 to xybins * xybins
  as.integer(nbins.y * (x-1) + y)
}

from_1dim = function(z, nbins.y) {
  list(x=as.integer(ceiling(z/nbins.y)), y=(z-1)%%nbins.y + 1)
}

get.quantiles.fun = function(quantile.fractions) {
  function(vals) {  quantile(vals, quantile.fractions) + 0.001 }
}

bin.responses = function(df, get.bin.thresholds.fun, binned.var='trace') {
  df = data.table(df)
  if (nrow(df) == 0) {
    return(df)
  }

  get.response.bin = function(vals, get.bin.thresholds.fun) {
    bin.thresholds = get.bin.thresholds.fun(vals)
    map_int(vals, ~ dplyr::first(which(.x <= bin.thresholds)), default=length(bin.thresholds))
  }

  binned.df = df[, response_bin := get.response.bin(.SD[[binned.var]], get.bin.thresholds.fun),
                 by=c('animal', 'date', 'cell_id')]
  binned.df
}

bin.nevents = function(df, nevents.thr=c(0.5, 2.0, 6.0)) {
  nevents.binning.fun = function(nevents) {
    nevents.thr
  }
  bin.responses(df, nevents.binning.fun, binned.var='ntrace')
}

bin.time.space = function(data.traces,
                          nbins.x,
                          nbins.y,
                          get.bin.thresholds.fun=NULL,
                          binned.var='trace',
                          timebin.dur.msec=200) {
  data.traces = data.table(data.traces)
  if (nrow(data.traces) == 0) {
    return(data.traces)
  }

  timebinned.traces = timebin.traces(data.traces[x >= 0 & y >= 0, ],
                                     timebin.dur.msec = timebin.dur.msec) %>%
    stimbin.traces(x, nbins.x) %>%
    stimbin.traces(y, nbins.y) %>%
    data.table()
  if (!is.null(get.bin.thresholds.fun)) {
    binned.traces = bin.responses(timebinned.traces, get.bin.thresholds.fun, binned.var=binned.var)
  } else {
    binned.traces = timebinned.traces
  }
  binned.traces[, bin.xy := to_1dim(bin.x, bin.y, nbins.y)]

  binned.traces
}

get.trial.ends = function(timestamps) {
  c(which(diff(timestamps) < 0), length(timestamps))
}

# Fast implementation of melting with data.tables
melt.traces = function(data) {
  .remove.cols.expr = function(data, cols.exp) {
    cols = grep(cols.exp, colnames(data))
    if (length(cols) > 0) {
      data[, grep(cols.exp, colnames(data)):=NULL]
    }
  }
  .remove.cols.expr(data, "^events_")
  .remove.cols.expr(data, "^dist")

  trace.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'trace_')]
  #events.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'events_')]
  deconv.measure.vars = colnames(data)[stringr::str_starts(colnames(data), 'deconvTrace_')]
  melted.df = melt(data,
       measure = list(trace.measure.vars,
                      #events.measure.vars,
                      deconv.measure.vars),
       value.name = c('trace',
                      #'nevents',
                      'deconv_trace'))

  melted.df = melted.df[, ('cell') := variable %>% trimws %>% as.integer][order(date,trial_id,cell,timestamp), !'variable']
  return(melted.df)
}

gather.traces = function(data) {
  data.traces = data %>%
    select(-starts_with('events_'), -starts_with('deconvTrace_')) %>%
    gather('cell', 'trace', starts_with('trace_')) %>%
    mutate(cell = str_replace(cell, 'trace_[.]?[.]?','')) %>%
    arrange(date, trial_id, cell, timestamp)

  data.deconv_traces = data %>%
    select('date', 'trial_id', 'timestamp', starts_with('deconvTrace_')) %>%
    gather('cell', 'deconv_trace', starts_with('deconvTrace_')) %>%
    mutate(cell = str_replace(cell, 'deconvTrace_[.]?[.]?','')) %>%
    arrange(date, trial_id, cell, timestamp)

  #data.events = data %>%
  #  select('date', 'trial_id', 'timestamp', starts_with('events_')) %>%
  #  gather('cell', 'nevents', starts_with('events_')) %>%
  #  mutate(cell = str_replace(cell, 'events_[.]?[.]?','')) %>%
  #  arrange(date, trial_id, cell, timestamp)

  data = cbind(data.traces,
               data.deconv_traces$deconv_trace)
               #data.events$nevents)
  #colnames(data)[(ncol(data)-1):ncol(data)] = c('deconv_trace', 'nevents')
  colnames(data)[(ncol(data))] = 'deconv_trace'
  data$cell = as.integer(data$cell)

  return(data)
}


# Filter based on velocity calculated for the time window. Filtering is done on epochs
# Slow code, use isRunning for better performance.
filter.running = function(df, min.run.velocity=2, mean.run.velocity=3, window.dur.ms=500) {
  df = data.table(df)
  setorder(df, trial_id, exp_title, cell_id, timestamp)
  is.running=rep(FALSE, nrow(df))
  i = 1
  while (i < nrow(df)) {
    j = i + 1
    while(j < nrow(df) &&
          df$velocity[i] >= min.run.velocity &&
          df$trial_id[i] == df$trial_id[j] &&
          df$velocity[j] >= min.run.velocity) {
      j = j + 1
    }
    j = j - 1

    dist = norm2(df$smooth_trans_x[j] - df$smooth_trans_x[i],
                 df$smooth_trans_y[j] - df$smooth_trans_y[i])
    dur.ms = max(df$timestamp[j] - df$timestamp[i], 1)
    vel = dist / dur.ms * 1000
    is.running[i:j] = (vel >= mean.run.velocity) && (dur.ms >= window.dur.ms)
    i = j + 1
  }

  df[which(is.running),]
}

# Map traces to principal components
pca.binned.traces = function(binned.traces) {
  response.matrix = reshape2::acast(binned.traces, time_bin ~ cell_id, value.var='response_bin')
  pca.res = prcomp(response.matrix, center = TRUE, scale. = TRUE)
  cumvar = cumsum(pca.res$sdev*pca.res$sdev)
  cumvar.portion = cumvar / cumvar[length(cumvar)]
  npc = which(cumvar.portion >= 0.9) %>% first
  response.df = as.data.frame(pca.res$x[,1:npc])
  response.df$time_bin = as.integer(rownames(response.matrix))
  pc.measure.vars = colnames(response.df)[stringr::str_starts(colnames(response.df), 'PC')]
  melt.pc.df = melt(response.df, measure=pc.measure.vars, value.name='mean.trace') %>%
    data.table()
  melt.pc.df[,cell_id:= as.integer(stringr::str_replace(variable, 'PC', ''))]
  metadata.df = dplyr::select(binned.traces, animal, date, trial_id, trial, time_bin, bin.xy, bin.x, bin.y) %>%
    distinct() %>% data.table()
  melt.pc.df = melt.pc.df[metadata.df, on='time_bin']

  quantile.fractions = c(0.2, 0.4, 0.6, 0.8, 1.0)
  binned.pc.df = bin.responses(melt.pc.df, get.quantiles.fun(quantile.fractions))
  setkey(binned.pc.df, time_bin, cell_id)
  setorder(binned.pc.df, time_bin, cell_id)

  return(binned.pc.df)
}
