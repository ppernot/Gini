# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('genefilter',version = "3.10")

fdgi = function(X, index=1:nrow(X),...){
  v1 = ErrViewLib::gini(X[index,1])
  v2 = ErrViewLib::gini(X[index,2])
  diff = v2 - v1
  return(diff)
}

dataSets = c(
  'BOR2019',
  'NAR2019',
  'PER2018',
  'SCH2018',
  'THA2015', # Need relative errors
  'WU2015', # Need relative errors
  'ZAS2019',
  'ZHA2018'
)
relSets = c('THA2015','WU2015') # Use relative errors

rm('dfg')
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
  
  # Mode-Centered data ####
  cErrors = Errors
  mo = list()
  for (i in 1:ncol(Errors)) {
    tag = paste0(set,':',meth)
    X = Errors[,i]
    mode = genefilter::half.range.mode(X,B=1000)
    mo[[tag]] = mode
    cErrors[, i] = X-mode
  }
  
  # Significance of Gini change when centering ####
  # Bootstrap Gini differences between raw and centered data
  dgi = udgi = gi = gic = si = list()
  for (meth in methList) { 
    tag = paste0(set,':',meth)
    bs = boot::boot(
      cbind(Errors[[meth]],cErrors[[meth]]),
      statistic = fdgi,
      R=1000)
    dgi[[tag]]  = bs$t0
    udgi[[tag]] = sd(bs$t,na.rm = TRUE)
    gi[[tag]] = ErrViewLib::gini(Errors[[meth]])
    gic[[tag]] = ErrViewLib::gini(cErrors[[meth]])
    si[[tag]] = sd(Errors[[meth]])
  }
  # Compute proba of diff to be null
  peq0 = list()
  for (meth in methList) {
    tag = paste0(set,':',meth)
    d  = dgi[[tag]]
    ud = udgi[[tag]]
    p = min(pnorm(0,d,ud),1-pnorm(0,d,ud))
    peq0[[tag]]  = p
    # cat(meth,'\t',prettyUnc(d,ud,2),'\t P(dgi == 0) = ',signif(p,2),'\n')
  }
  df = data.frame(d = unlist(dgi), ud = unlist(udgi), p = unlist(peq0),
                  gi = unlist(gi), gic = unlist(gic), mo = unlist(mo),
                  si = unlist(si))
  if(!exists('dfg'))
    dfg = df
  else
    dfg = rbind(dfg,df)
}
# Estimate percentages 
pdneg = mean(dfg$p < 0.05 & dfg$d < 0)
pdnul = mean(dfg$p >= 0.05)
pdpos = mean(dfg$p < 0.05 & dfg$d > 0)

