source('./0-Setup.R')

# Cols / lty scheme (max: 21)
tcols = rep(cols, 3)
tlty  = c(rep(1, length(cols)),
          rep(2, length(cols)),
          rep(3, length(cols)))

# Loop over datasets ----
for (set in dataSets) {
  cat('\nData set : ', set, '\n')
  
  ## Get data ----
  data = read.csv(
    file = file.path('..', 'data', paste0(set, '_Data.csv')),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  systems <- data[, 1]
  rownames(data) = systems
  Ref <- data[, 2]
  Data <- data[, -c(1, 2), drop = FALSE]
  methList <- colnames(Data)
  Errors <- Ref - Data
  if (set %in% relSets)
    Errors <- 100 * Errors / Ref
  
  # Split large sets
  nMeth = length(methList)
  nMax  = 18
  nPlot = floor(nMeth / nMax)
  if (nMeth %% nMax != 0)
    nPlot = nPlot + 1
  
  ## Plot ----
  frow = 1
  fcol = 2
  for (iPlot in 1:nPlot) {
    first = (iPlot - 1) * nMax + 1
    last  = min(nMeth, iPlot * nMax)
    
    fname = paste0('Fig_SI_', set, '.png')
    if (iPlot > 1)
      fname = paste0('Fig_SI_', set, '_', iPlot, '.png')
    
    png(
      file = file.path('..', 'results', 'figs', fname),
      width = fcol * reso,
      height = frow * reso
    )
    par(
      mfrow = c(frow, fcol),
      pty = pty,
      mar = mar,
      mgp = mgp,
      tcl = tcl,
      lwd = lwd,
      cex = cex
    )
    
    ### (a) ----
    plot(
      ecdf(abs(Errors[, first])),
      do.points = FALSE,
      lwd = lwd,
      col = tcols[1],
      lty = tlty[1],
      xlim = c(0, hd(abs(unlist(
        Errors
      )), 0.99)),
      xlab = paste0('|Errors| (', units[set], ')'),
      ylim = c(0, 1),
      yaxs = 'i',
      ylab = 'Probability',
      main = 'ECDF'
    )
    grid(col = 'gray70')
    for (i in (first + 1):last)
      lines(
        ecdf(abs(Errors[, i])),
        do.points = FALSE,
        lty = tlty[i - first + 1],
        lwd = lwd,
        col = tcols[i - first + 1]
      )
    legend(
      'bottomright',
      bty = 'n',
      cex = 0.5,
      legend = methList[first:last],
      col = tcols,
      lty = tlty
    )
    box()
    mtext(
      text = '(a)',
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3
    )
    
    ### (b) ----
    plot(ineq::Lc(abs(Errors[, first])),
         lwd = lwd,
         col = tcols[1],
         lty = tlty[1])
    grid(col = 'gray70')
    X = rnorm(10000)
    lines(ineq::Lc(abs(X)),
          lty = 2,
          lwd = lwd,
          col = 'gray70')
    for (i in (first + 1):last)
      lines(ineq::Lc(abs(Errors[, i])),
            lwd = lwd,
            col = tcols[i - first + 1],
            lty = tlty[i - first + 1])
    box()
    mtext(
      text = '(b)',
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3
    )
    
    dev.off()
    
  }
}