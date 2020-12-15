source('./0-Setup.R')

rm('dft')

stats = c(
  'mue',
  'mse',
  'rmsd',
  'q95hd',
  'hrmode',
  'skew',
  'skewgm',
  'skewgm_mcf',
  'kurt',
  'kurtcs',
  'gini',
  'gimc',
  'bmax',
  'gmax'
)

for (removeGlobalOutliers in c(FALSE, TRUE))
  for (correctTrendDegree in c(-1, 1))
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
      if(set %in% relSets)
        Errors <- 100 * Errors / Ref
      
      if (removeGlobalOutliers) {
        # Global outliers  (out of 95% CI)  ####
        qlim = t(apply(Errors, 2, quantile, 
                       probs = c(0.025, 0.975)))
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
      
      if (correctTrendDegree >= 0) {
        if(correctTrendDegree == 0) {
          # Center ####
          for (i in 1:ncol(Errors)) {
            x = Data[, i]
            y = Errors[, i]
            y = residuals(lm(y ~ 1))
            Errors[, i] = y
          }
          colnames(Errors) = paste0('c-', colnames(Errors))
          
        } else {
          # Correct linear trend ####
          for (i in 1:ncol(Errors)) {
            x = Data[, i]
            y = Errors[, i]
            y = residuals(lm(y ~ x))
            Errors[, i] = y
          }
          colnames(Errors) = paste0('lc-', colnames(Errors))
        }
      } else if (correctTrendDegree == -2) {
        # Mode centering
        for (i in 1:ncol(Errors)) {
          y = Errors[, i]
          Errors[, i] = y - hrmode(y) 
        }
        colnames(Errors) = paste0('mc-', colnames(Errors))
      }
      
      # Estimate stats ####
      bs = estBS1(
        Errors,
        props = stats,
        eps = 1,
        do.sip = FALSE,
        silent = TRUE
      )
      df = data.frame(Dataset = set, Methods = methList)
      df[, 'ctd'] = correctTrendDegree
      df[, 'out'] = removeGlobalOutliers
      for (stat in stats) {
        df[, stat] = signif(bs[[stat]]$val, 6)
        df[, paste0('u_', stat)] = signif(bs[[stat]]$unc, 3)
      }
      
      if(!exists('dft'))
        dft = df
      else
        dft = rbind(dft,df)
      
    }

# Save data
write.csv(dft,file = resultsTab, row.names = FALSE)
