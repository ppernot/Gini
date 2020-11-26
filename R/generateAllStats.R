library(ErrViewLib)

lcp = function(X, index = 1:length(X), p=0.95, ...) {
  # Lorenz curve
  X = sort(abs(X[index]))
  pr = (1:length(X)) / length(X)
  lc = cumsum(X)/sum(X)
  # Linear interpolation
  iup = which(pr >= p)[1]
  ilo = iup - 1
  Lcp = lc[ilo]+ (p-pr[ilo])*
    (lc[iup]-lc[ilo])/(pr[iup]-pr[ilo])

  return(1-Lcp)
}

# lcpmue = function(X, index = 1:length(X), ...) {
#   # Lorenz curve
#   X = sort(abs(X[index]))
#   pr = (1:length(X)) / length(X)
#   mue = mean(X)
#   p = pr[which(X >= mue)[1]]
#   lc = cumsum(X)/sum(X)
#   # Linear interpolation
#   iup = which(pr >= p)[1]
#   ilo = iup - 1
#   Lcp = lc[ilo]+ (p-pr[ilo])*
#     (lc[iup]-lc[ilo])/(pr[iup]-pr[ilo])
#
#   return(1-Lcp)
# }


pmue = function(X, index = 1:length(X), ...) {
  # P(abs(X) >= 2*MUE)
  X = X[index]
  sum(abs(X) >= 3*mue(X))/length(X)
}
p2med = function(X, index = 1:length(X), ...) {
  X = abs(X[index])
  sum(X >= 2*median(X))/length(X)
}

resultsTab = file.path('..', 'results', 'tables', 'allStatsXXX.csv')
# if (file.exists(resultsTab))
#   file.remove(resultsTab)
rm('dft')

dataSets = c(
  'Ref_NormBias',
  'Ref_Student',
  'Ref_StudentBias1',
  'Ref_StudentBias2',
  'Ref_GandH',
  'Ref_GandHBias1',
  'Ref_GandHBias2',
  'BOR2019',
  'HAI2018',
  'NAR2019',
  'PER2018',
  'SCH2018',
  'THA2015', # Need relative errors
  'WU2015', # Need relative errors
  'ZAS2019',
  'ZHA2018'
)
relSets = c('THA2015','WU2015') # Use relative errors

for (removeGlobalOutliers in c(FALSE, TRUE))
  for (correctTrendDegree in c(-1,0,1)[3])
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
        Errors <- Errors / Ref

      if (removeGlobalOutliers) {
        # Global outliers  (out of 95% CI)  ####
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
      }

      # Estimate stats ####
      stats = c(
        'mue',
        'mse',
        'rmsd',
        'skew',
        'kurt',
        'kurtcs',
        'q95hd',
        'gini',
        'lasym',
        'pietra'
      )
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
        df[, stat] = signif(bs[[stat]]$val, 5)
        df[, paste0('u_', stat)] = signif(bs[[stat]]$unc, 2)
      }
      # rownames(df) = NULL

      if(!exists('dft'))
        dft = df
      else
        dft = rbind(dft,df)

    }

# Save data
write.csv(dft,file = resultsTab, row.names = FALSE)

stop()

# COmbine with previous results
# resultsTab = file.path('..', 'results', 'tables', 'allStats.csv')
# D = read.csv(file= resultsTab)
# D = cbind(D,dft[,5:6])
# write.csv(D,file = resultsTab, row.names = FALSE)



for(dat in unique(dft$Dataset)) {
  sel = dft$Dataset == dat
  plot(dft[sel,'q95hd'],
       dft[sel,'lcp'],
       main = dat,
       col = cols[6])
  text(dft[sel,'q95hd'],
       dft[sel,'lcp'],
       dft[sel,'Methods'])
}



# Shape ####
png(file='../article/fig_corr_shapestats_cs_centered.png', width=1600, height=1600)
# par(
#   mfrow = c(1, 1),
#   pty = pty,
#   mar = mar,
#   mgp = mgp,
#   tcl = tcl,
#   lwd = lwd,
#   cex = cex
SAPlot = function(X,cex=1) {
  if (missing(X)) {
    print('Scatterplot and correlation pairs for sample X',quote=F)
    print(" ",quote=F)
    print("Call : SAPlot(X)",quote=F)
    print("-X   : (oblig) a MxN matrix of M values for N variables",quote=F)
    print("-cex : graphical parameter",quote=F)
    return(invisible())
  }
  sdX=apply(X,2,sd) # Identify fixed params to exclude from plot
  par(cex=cex,cex.axis=4, mai = c(1,1,1,1),
      mgp  = c(3, 2, 0),
      tcl  = -0.5)
  pairs(X[,sdX != 0], gap=0, cex.labels = 6,
        upper.panel=panel.cor,
        diag.panel =panel.hist,
        lower.panel=panel.smooth )
}
panel.hist <- function(x,...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks;
  nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  grid(col='brown',ny=0)
  rect(breaks[-nB],0,breaks[-1],y,col="orange",...)
}
panel.cor <- function(x,y,digits=2,prefix="",cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,method="spearman")
  ra = abs(r)
  txt <- format(c(r,0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)*ra
  text(0.5,0.5,txt,cex = cex.cor,col=ifelse(r>=0,4,2))
}
panel.smooth <- function (x, y, cex = 4, col.smooth = "red",
                          span = 2/3, iter = 3, ...) {
  maxPoints=500
  nP=min(maxPoints,length(x))
  iSamp = seq.int(1,length(x),length.out=nP)
  x1=x[iSamp]
  y1=y[iSamp]
  green_tr=rgb(unlist(t(col2rgb("darkgreen"))),
               alpha=150,maxColorValue = 255)
  grid(col='brown')
  points(x1, y1, pch = 19, col = green_tr, lwd=0, cex = cex)
}
sel = !dft$out & dft$ctd == -1
M = cbind(dft$kurtcs[sel], boot::logit(dft$W)[sel],
          dft$lcp[sel],dft$gini[sel])
colnames(M)=c('KurtCS','logit(W)','l95','G')
SAPlot(M,cex=1)
# sel = dft$kurt < 10
# SAPlot(M[sel,],cex=2)
dev.off()
stop()

# dft = read.csv(file=resultsTab)
sel1 = !dft$out & dft$ctd == -1
sel2 = dft$out & dft$ctd == -1
plot(dft$kurtcs[sel],dft$lcp[sel],type = 'n', cex=0.1)
segments(dft$kurtcs[sel1],dft$lcp[sel1],dft$kurtcs[sel2],dft$lcp[sel2],
         col=as.numeric(factor(dft$Dataset[sel1])))
abline(v=0);abline(h=0.15)
text(dft$kurtcs[sel1],dft$lcp[sel1],dft$Methods[sel1],cex=0.5,
     col=as.numeric(factor(dft$Dataset[sel1])))
