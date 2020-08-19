.norm2 = function(x, y) {
  sqrt(x**2 + y**2)
}

# Create df with smoothed values from matrix representation
create.pf.df = function(M, occupancyM, min.occupancy.sec=1, frame.rate=20, sigma=1.4, 
                        smooth=FALSE, filter.circular=TRUE) {
  if (smooth) {
    M = gauss2dsmooth(M, lambda=sigma, nx=11, ny=11)
  }
  df1 = reshape2::melt(M) 
  
  min.occupancy = min.occupancy.sec * frame.rate
  
  smoothedOccupancy = occupancyM
  min.smoothed.occupancy = 0
  if (smooth) {
    min.smoothed.occupancy = min.occupancy * 1/(2*pi*sigma^2)
    smoothedOccupancy = gauss2dsmooth(occupancyM, lambda=sigma, nx=11, ny=11)
  }
  
  df_org = reshape2::melt(smoothedOccupancy) %>%
    dplyr::filter(value >= min.smoothed.occupancy) 
  if (filter.circular) {
    mid.pt = (dim(M)[1] + 1) / 2 
    df_org = dplyr::filter(df_org, .norm2(Var1 - mid.pt, Var2 - mid.pt) <= (mid.pt + 0.25))
  }
  df2 = left_join(df_org, df1, by=c('Var1'='Var1', 'Var2'='Var2'), suffix=c('.occupancy', '.field'))
  
  return(df2)
}

geom_pf = function(df, max.x, max.y) {
  jet.colours = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  list(
    geom_raster(data=filter(df, !is.na(value.field)),
                mapping=aes(x=Var1, y=max.y-Var2, fill=value.field), interpolate=FALSE),
    scale_fill_gradientn(colours=jet.colours(7)),
    xlim(c(0, max.x)),
    ylim(c(0, max.y)),
    theme_void()
  )
}

jet.colours = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

plot.pf = function(df, max.x=NA, max.y=NA) {
  g = df %>%
    filter(!is.na(value.field)) %>%
    ggplot(aes(x=Var1, y=max.y-Var2)) +
    geom_raster(aes(fill=value.field), interpolate=FALSE) +
    scale_fill_gradientn(colours=jet.colours(7)) +
    theme_void() 
  
  if (!is.na(max.x)) {
    g = g + xlim(c(0, max.x))
  }
  if (!is.na(max.y)) {
    g = g + ylim(c(0, max.y))
  }
  
  return(g)
}

plot.trace.events = function(df, event.var=nevents, event.thr=0.5) {
  event.var = enquo(event.var)
  df.events = filter(df, !!event.var >= event.thr)
  mid.pt = c(50,50)
  df %>%
    arrange(timestamp) %>%
    filter( sqrt((mid.pt[1] - x)^2 + (mid.pt[1] - y)^2) < mid.pt[1]) %>%
    ggplot(aes(x=x, y=100-y)) +
    geom_path(aes(group=trial_id), alpha=0.6) +
    geom_point(data=df.events, color='red') +
    xlim(c(0, 100)) + ylim(c(0, 100)) +
    theme_void()
}

plot.event.vectors = function(df, event.var=nevents, event.thr=0.5) {
  event.var = enquo(event.var)
  arrow.len = 5
  df$is_running = as.logical(df$is_running)
  df.events = dplyr::filter(df, !!event.var >= event.thr) %>%
    dplyr::mutate(arrow.x = x + cos(angle / 180 * pi) * arrow.len,
                  arrow.y = y + sin(angle / 180 * pi) * arrow.len)
  ggplot(df.events, aes(x=x, y=100-y)) +
    geom_segment(aes(x=x, y=100-y, xend=arrow.x, yend=100-arrow.y, color=is_running),
                 arrow = ggplot2::arrow(length = unit(0.015, "npc"))) +
    xlim(c(0, 100)) + ylim(c(0, 100)) +
    theme_void()
}


field.cor = function(field1, field2, max.xy, make.cor.plot=FALSE) {
  result = list()
  joined.fields = inner_join(field1, field2, by=c('Var1', 'Var2')) %>%
    filter(!is.na(value.field.x), !is.na(value.field.y)) %>%
    # Avoid calculating correlation on the edges: the blue highly correlated and increases the overall correlation
    filter(Var1 < max.xy - 1, Var2 < max.xy - 1, Var1 > 2, Var2 > 2)  
  result$cor = cor(joined.fields$value.field.x, joined.fields$value.field.y)
  
  if (make.cor.plot) {
    m.field.x = mean(joined.fields$value.field.x)
    sd.field.x = sd(joined.fields$value.field.x)
    m.field.y = mean(joined.fields$value.field.y)
    sd.field.y = sd(joined.fields$value.field.y)
    
    joined.fields = mutate(
      joined.fields,
      value.field = (value.field.x - m.field.x) * (value.field.y - m.field.y) /
        (nrow(joined.fields) - 1) / sd.field.x / sd.field.y)
    
    g = ggplot(joined.fields) +
      geom_raster(aes(x=Var1, y=max.xy-Var2, fill=value.field), interpolate=FALSE) +
      scale_fill_gradient2(low = 'blue', mid = 'white', high='red', midpoint = 0.0) +
      xlim(c(0, max.xy)) + ylim(c(0, max.xy)) +
      theme_void()
    result$g = g
  }
  
  return(result)
}


cell.spatial.info = function(cell.df, 
                             nstim.x,
                             nstim.y,
                             generate.plots=FALSE, 
                             nshuffles=0,
                             bin.hz=5,
                             shuffle.shift.sec=10,
                             trace.var='trace',
                             binned.trace.var='response_bin',
                             min.occupancy.sec=1,
                             kernel.size=9,
                             gaussian.var=2) {
  cell.df = data.table(cell.df)
  
  cell.nevents = nrow(cell.df[nevents > 0,])
  trace.vals = cell.df[[trace.var]]

  if (length(trace.vals) == 0) {
    return(list(cell_info=list(),
                field=matrix(),
                occupancy=matrix(),
                g=NULL))
  }
  cell_name = cell.df$cell_id[1]
  trial_ends = get.trial.ends(cell.df$timestamp)
  trace.quantiles = quantile(trace.vals, c(0.2, 0.99, 1.0), na.rm=TRUE)

  nstim = nstim.x * nstim.y
  pf = calcPlaceField(bin_x=cell.df$bin.x, 
                      bin_y=cell.df$bin.y, 
                      nbins_x=nstim.x, 
                      nbins_y=nstim.y, 
                      trace=trace.vals,
                      binnedTrace=cell.df[[binned.trace.var]],
                      minOccupancy=min.occupancy.sec * bin.hz,
                      kernelSize=kernel.size,
                      gaussianVar=gaussian.var)
  
  shuffle.pf = placeFieldStatsForShuffled(bin_x=cell.df$bin.x, 
                                          bin_y=cell.df$bin.y, 
                                          nbins_x=nstim.x, 
                                          nbins_y=nstim.y, 
                                          trace=trace.vals,
                                          binnedTrace=cell.df[[binned.trace.var]],
                                          trialEnds=as.integer(trial_ends), 
                                          nshuffles=nshuffles,
                                          minShift=shuffle.shift.sec * bin.hz,
                                          minOccupancy=min.occupancy.sec * bin.hz, # min occupancy
                                          kernelSize=kernel.size,
                                          gaussianVar=gaussian.var)
  si.signif.thresh = quantile(shuffle.pf$shuffle.si, 0.95, na.rm=TRUE)[[1]]
  si.signif = pf$spatial.information >= si.signif.thresh
  mi.signif.thresh = quantile(shuffle.pf$shuffle.mi, 0.95, na.rm=TRUE)[[1]]
  mi.signif = pf$mutual.info >= mi.signif.thresh

  cell_info = list(cell_id=cell_name,
                   spatial.information=pf$spatial.information,
                   si.signif.thresh=si.signif.thresh,
                   signif.si=si.signif,
                   mi.signif.thresh=mi.signif.thresh,
                   signif.mi=mi.signif,
                   field.size.50=pf$field.size.50,
                   field.size.25=pf$field.size.25,
                   spatial.information.perspike=pf$spatial.information.perspike,
                   trace.mean=pf$mfr,
                   mutual.info=pf$mutual.info,
                   mutual.info.bias=pf$mutual.info.bias,
                   space.sampling.factor = pf$space.sampling.factor,
                   sparsity = pf$sparsity,
                   quantile20=trace.quantiles[["20%"]],
                   quantile99=trace.quantiles[["99%"]],
                   nevents=cell.nevents)

  # find field max value and pos in the smoothed values
  pf.df = create.pf.df(pf$field, 
                       pf$occupancy, 
                       min.occupancy.sec=min.occupancy.sec, 
                       frame.rate=bin.hz, 
                       smooth=FALSE) %>%
    filter(!is.na(value.field))
  if (nrow(pf.df) > 0) {
    max.row = pf.df[which.max(pf.df$value.field),]
  } else { # no bin with high enough occupancy
    max.row=list(value.field=0, Var1=-1, Var2=-1)
  }
  cell_info$field.max = max.row$value.field
  cell_info$field.mean = mean(pf.df$value.field, na.rm = TRUE)
  cell_info$field.max.x = max.row$Var1 / nstim.x * 100
  cell_info$field.max.y = max.row$Var2/ nstim.y * 100
  
  g.placefield=NA
  if (generate.plots) {
    cell_event_rate = cell.nevents / nrow(cell.df) * bin.hz

    g.placefield = plot.pf(pf.df, nstim.x, nstim.y) +
      labs(title=paste0('Cell ', cell_name, ' MER = ', format(cell_event_rate, digits=2), ' Hz',
                        '\nSI = ', format(pf$spatial.information, digits=2),
                        ifelse(si.signif, '*', ''),
                        ' SI per spike = ', format(pf$spatial.information.perspike, digits=2),
                        '\nMI - bias = ', format(pf$mutual.info - pf$mutual.info.bias, digits=3),
                        ifelse(mi.signif, '*', ''),
                        '\nsparsity = ', format(pf$sparsity, digits=2)),
           fill=paste(trace.var, 'values')) +
      theme(text = element_text(size=4),
            plot.title = element_text(size=4))
  }

  return(list(cell_info=cell_info,
              field=pf$field,
              occupancy=pf$occupancy,
              g=g.placefield))
}
