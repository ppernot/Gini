
gh = function(N = 1000,
              g = 0,
              h = 0) {
  x = rnorm(N)
  if (g == 0)
    x = x * exp(h * x ^ 2 / 2)
  if (g > 0)
    x = (exp(g * x) - 1) * exp(h * x ^ 2 / 2) / g
  return(x)
}
hrmode = function(X, index = 1:length(X), ...) {
  genefilter::half.range.mode(X[index])
}

skewgm = function(X, index = 1:length(X), ...) {
  X = X[index]
  m = hd(X,0.5)
  s = (mean(X) - m) / mean(abs(X - m))
  return(s)
}

skewgmf = function(X, index = 1:length(X), ...) {
  X = X[index]
  X = abs(X - hrmode(X)) # Absolute mode-centered errors
  m = hd(X,0.5) # Median
  s = (mean(X) - m) / mean(abs(X - m))
  return(s)
}

gimc = function(X, index = 1:length(X), ...) {
  X = X[index]
  ErrViewLib::gini(X - hrmode(X))
}

cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1.5, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 4.5

png(file = file.path('..', 'article', 'fig_G_vs_Skew_mode.png'), 
    width=2400, height=2400)
par(
  mfrow = c(2, 2),
  mar = mar,
  pty = pty,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

# a ----

N = 1e5
gi = sk = c()
i = 0
dmax = 0
for (g in seq(0,1,by=0.1)) {
  # Asymmetry
  for (h in seq(0,0.5,by=0.1)) {
    # Tails
    i = i + 1
    X = gh(N, g, h)
    X = X - hrmode(X)
    gi[i] = ErrViewLib::gini(X)
    sk[i] = skewgm(abs(X))
  }
}

plot(
  gi, sk, pch=16, col = cols[2],
  xlim = c(0.4,0.7), xlab = 'G',
  ylim = c(0.2,0.8), ylab = expression(beta[GMF]),
  main = paste0('Ref / N=',N)
)
grid()

i=0
gi = sk = c()
for (nu in seq(2,20,by=1)) {
  i = i + 1
  X = rt(N, df=nu)
  X = X - hrmode(X)
  gi[i] = ErrViewLib::gini(X)
  sk[i] = skewgm(abs(X))
}

points(
  gi, sk, 
  pch=17,
  col= cols[6]
)
legend(
  'topleft', bty = 'n',
  legend = c('g-and-h','Student'),
  col = cols[c(2,6)],
  pch = c(16,17)
)
box() 
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

# b ----

N = 150
gi = sk = c()
i = 0
dmax = 0
for (g in seq(0,1,by=0.1)) {
  # Asymmetry
  for (h in seq(0,0.5,by=0.1)) {
    # Tails
    i = i + 1
    X = gh(N, g, h)
    X = X - hrmode(X)
    gi[i] = ErrViewLib::gini(X)
    sk[i] = skewgm(abs(X))
  }
}
plot(
  gi, sk, pch=16, col = cols[2],
  xlim = c(0.4,0.7), xlab = 'G',
  ylim = c(0.2,0.8), ylab = expression(beta[GMF]),
  main = paste0('Ref / N=',N)
)
grid()

i=0
gi = sk = c()
for (nu in seq(2,20,by=1)) {
  i = i + 1
  X = rt(N, df=nu)
  X = X - hrmode(X)
  gi[i] = ErrViewLib::gini(X)
  sk[i] = skewgm(abs(X))
}

points(
  gi, sk, 
  pch=17,
  col= cols[6]
)
box() 
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)


# c ----
dataSets = c(
  'BOR2019',
  'NAR2019',
  'PER2018',
  'SCH2018',
  'THA2015', # Need relative errors
  'WU2015',  # Need relative errors
  'ZAS2019',
  'ZHA2018'
)
relSets = c('THA2015','WU2015') # Use relative errors
units = c('eV','kcal/mol','kcal/mol',
          'eV','%','%','kcal/mol','eV/atom')
names(units) = dataSets

tcols = rep(cols,6)

rm('dft')
for (set in dataSets) {
  cat('\nData set : ', set, '\n')
  # Get data
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
  
  # Bootstrapping
  stats = c( 'gimc','skewgmf')
  bs = estBS1(
    Errors,
    props = stats,
    eps = 1,
    do.sip = FALSE,
    silent = TRUE
  )
  df = data.frame(Dataset = set, Methods = methList)
  for (stat in stats) {
    df[, stat] = bs[[stat]]$val
    df[, paste0('u_', stat)] = bs[[stat]]$unc
  }
  if(!exists('dft'))
    dft = df
  else
    dft = rbind(dft,df)
}

icol = as.numeric(factor(dft$Dataset)) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(dft$Dataset))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17

plot(
  dft$gimc, dft$skewgmf,
  col=cols[icol], 
  pch=pch,
  xlim = c(0.4,0.7), xlab = 'G',
  ylim = c(0.2,0.8), ylab = expression(beta[GMF]),
  main = 'Literature'
)
grid()

pch = unique(as.numeric(factor(dft$Dataset)))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
legend(
  'topleft', bty = 'o', box.col = 'white', ncol = 1,
  legend = paste0('mc-',unique(dft$Dataset)),
  col = unique(cols[icol]),
  pch = pch,
  cex=0.8
)
box()
mtext(
  text = '(c)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

# d ----
ux = dft$u_gimc
uy = dft$u_skewgmf
plot(
  ux,uy,
  pch=16,
  col = cols[icol],
  xlim = c(0.,0.03), xlab = 'u(G)',
  ylim = c(0.,0.12), ylab = expression(u(beta[GMF]))
)
grid()
# abline(a=0,b=1,lty=1,col='gray90')
abline(a=0,b=2,lty=2,col='gray50')
abline(a=0,b=5,lty=3,col='gray50')
box()
mtext(
  text = '(d)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()