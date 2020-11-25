library(ErrViewLib)

removeGlobalOutliers = FALSE
correctTrend = FALSE

resultsTab = file.path('..', 'results', 'tables', 'allStats.csv')
if (file.exists(resultsTab))
  file.remove(resultsTab)

dataSets = c('BOR2019',
             'PER2018')

for (removeGlobalOutliers in c(FALSE, TRUE))
  for (correctTrend in c(FALSE, TRUE))
    for (set in dataSets) {
      cat('\nData set : ', set, '\n')

      # Get data ####
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

      if (correctTrend) {
        # Correct linear trend ####
        for (i in 1:ncol(Errors)) {
          x = Data[, i]
          y = Errors[, i]
          y = residuals(lm(y ~ x))
          Errors[, i] = y
        }
        colnames(Errors) = paste0('lc-', colnames(Errors))
      }

      if (removeGlobalOutliers) {
        # GLobal outliers  (out of 95% CI)  ####
        qlim = t(apply(Errors, 2, quantile, probs = c(0.025, 0.975)))
        ql1  = matrix(
          qlim[, 1],
          ncol = ncol(Errors),
          nrow = nrow(Errors),
          byrow = TRUE
        )
        ql2  = matrix(
          qlim[, 2],
          ncol = ncol(Errors),
          nrow = nrow(Errors),
          byrow = TRUE
        )
        out  = rowSums(Errors < ql1 | Errors > ql2) == ncol(Errors)
        if (any(out)) {
          cat('- Outliers (CI95) : ', systems[out], '\n')
          Errors = Errors[!out,]
          Data   = Data[!out,]
          systems = systems[!out]
        }
      }

      # Estimate stats ####
      stats = c('mue',
                'mse',
                'rmsd',
                'skew',
                'kurt',
                'q95hd',
                'W',
                'gini',
                'lasym')
      bs = estBS1(
        Errors,
        props = stats,
        eps = 1,
        do.sip = FALSE,
        silent = TRUE
      )
      df = data.frame(Dataset = set, Methods = methList)
      df[, 'lc']  = correctTrend
      df[, 'out'] = removeGlobalOutliers
      for (stat in stats) {
        df[, stat] = signif(bs[[stat]]$val, 5)
        df[, paste0('u_', stat)] = signif(bs[[stat]]$unc, 2)
      }
      rownames(df) = NULL

      # Save data
      write.table(
        df,
        file = resultsTab,
        sep = ',',
        append = TRUE,
        col.names = set == dataSets[1],
        row.names = FALSE
      )
    }
