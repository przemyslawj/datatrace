.norm2 = function(x, y) {
  sqrt(x**2 + y**2)
}

.to.matrix = function(M, nstim.x, nstim.y) {
  if (!is.matrix(M)) {
    M = matrix(M, nrow=nstim.x, ncol=nstim.y, byrow=TRUE) 
  }
  return(M)
}

# Create df with smoothed values from matrix representation
create.pf.df = function(M, occupancyM, min.occupancy.sec=1, frame.rate=20, sigma = 1.4) {
  M1 = gauss2dsmooth(M, lambda=sigma, nx=11, ny=11)
  df1 = reshape2::melt(M1) 
  
  min.occupancy = min.occupancy.sec * frame.rate
  min.smoothed.occupancy = min.occupancy * 1/(2*pi*sigma^2)
  smoothedOccupancy = gauss2dsmooth(occupancyM, lambda=sigma, nx=11, ny=11)
  #TODO: filtering below assumes a circular maze
  mid.pt = mean(1:dim(M)[1])
  df_org = reshape2::melt(smoothedOccupancy) %>%
    dplyr::filter(value >= min.smoothed.occupancy) %>%
    dplyr::filter(.norm2(Var1 - mid.pt, Var2 - mid.pt) <= mid.pt)
  df2 = left_join(df_org, df1, by=c('Var1'='Var1', 'Var2'='Var2'), suffix=c('.occupancy', '.conv'))
  
  return(df2)
}
  
plot.pf = function(df, max.x, max.y) {
  jet.colours = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  ggplot(df, aes(x=Var1, y=max.y-Var2)) +
    geom_raster(aes(fill=value.conv), interpolate=FALSE) +
    scale_fill_gradientn(colours=jet.colours(7)) +
    xlim(c(0, max.x)) + ylim(c(0, max.y)) +
    theme_void() 
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
    #facet_grid(trial_id ~ cell_id) +
    xlim(c(0, 100)) + ylim(c(0, 100)) +
    theme_void()
}


field.cor = function(field1, field2, max.xy, make.cor.plot=FALSE) {
  result = list()
  joined.fields = inner_join(field1, field2, by=c('Var1', 'Var2')) %>%
    # Avoid calculating correlation on the edges: the blue highly correlated and increases the overall correlation
    filter(Var1 < max.xy - 1, Var2 < max.xy - 1, Var1 > 2, Var2 > 2)  
  result$cor = cor(joined.fields$value.conv.x, joined.fields$value.conv.y)
  
  if (make.cor.plot) {
    m.conv.x = mean(joined.fields$value.conv.x)
    sd.conv.x = sd(joined.fields$value.conv.x)
    m.conv.y = mean(joined.fields$value.conv.y)
    sd.conv.y = sd(joined.fields$value.conv.y)
    
    joined.fields = mutate(
      joined.fields,
      value.conv = (value.conv.x - m.conv.x) * (value.conv.y - m.conv.y) /
        (nrow(joined.fields) - 1) / sd.conv.x / sd.conv.y)
    
    g = ggplot(joined.fields) +
      geom_raster(aes(x=Var1, y=max.xy-Var2, fill=value.conv), interpolate=FALSE) +
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
                             trace.var='trace',
                             binned.trace.var='response_bin',
                             min.occupancy.sec=1) {
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
  pf = calcPlaceField(as.integer(to_1dim(cell.df$bin.x, cell.df$bin.y, nstim.y)), 
                      nstim,
                      trace.vals,
                      cell.df[[binned.trace.var]],
                      as.integer(trial_ends), 
                      nshuffles,
                      2 * bin.hz,
                      min.occupancy.sec * bin.hz) # min occupancy
  pf.field = .to.matrix(pf$field, nstim.x, nstim.y)
  pf.occupancy = .to.matrix(pf$occupancy, nstim.x, nstim.y)

  si.signif.thresh = quantile(pf$shuffle.si, 0.95, na.rm=TRUE)[[1]]
  si.signif = pf$spatial.information >= si.signif.thresh
  mi.signif.thresh = quantile(pf$shuffle.mi, 0.95, na.rm=TRUE)[[1]]
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
                   mfr=pf$mfr,
                   mutual.info=pf$mutual.info,
                   mutual.info.bias=pf$mutual.info.bias,
                   space.sampling.factor = pf$space.sampling.factor,
                   sparsity = pf$sparsity,
                   quantile20=trace.quantiles[["20%"]],
                   quantile99=trace.quantiles[["99%"]],
                   nevents=cell.nevents)

  # find field max value and pos in the smoothed values
  pf.df = create.pf.df(pf.field, pf.occupancy, frame.rate=bin.hz)
  if (nrow(pf.df) > 0) {
    max.row = pf.df[which.max(pf.df$value.conv),]
  } else { # no bin with high enough occupancy
    max.row=list(value.conv=0, Var1=-1, Var2=-1)
  }
  cell_info$field.max = max.row$value.conv
  cell_info$field.mean = mean(pf.df$value.conv)
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
             field=pf.field,
             occupancy=pf.occupancy,
             g=g.placefield))
}
