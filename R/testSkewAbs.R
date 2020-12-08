cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1.5, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 4.5


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

skewgm_f = function(X, index = 1:length(X), ...) {
  return(skewgm(abs(X[index])))
}

skewgm_mcf = function(X, index = 1:length(X), ...) {
  X = X[index]
  X = abs(X - hrmode(X)) # Absolute mode-centered errors
  return(skewgm(X))
}

kurtcs_f = function(X, index = 1:length(X), ...) {
  return( ErrViewLib::kurtcs( abs(X[index]) ) )
}

kurtcs_mcf = function(X, index = 1:length(X), ...) {
  X = X[index]
  X = abs(X - hrmode(X)) # Absolute mode-centered errors
  return( ErrViewLib::kurtcs(X) )
}

gimc = function(X, index = 1:length(X), ...) {
  X = X[index]
  ErrViewLib::gini(abs(X - hrmode(X)))
}

gimc_pm = function(X, index = 1:length(X), ...) {
  X = X[index]
  X = X - hrmode(X)
  Gm = ErrViewLib::gini(X[X<=0]) 
  Gp = ErrViewLib::gini(X[X>=0]) 
  return(max(Gm,Gp))
}


# Pre-calculate stats for lit datasets ----

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
  stats = c( 
    'gini',
    'gimc',
    'gimc_pm',
    'skewgm',
    'skewgm_f',
    'skewgm_mcf',
    'kurtcs',
    'kurtcs_f',
    'kurtcs_mcf')

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

# Precalculate for G vs. Skew  ----
N = 1e6
gi = sk = c()
i = 0
for (g in seq(0,1,by=0.1)) {
  # Asymmetry
  for (h in seq(0,0.5,by=0.1)) {
    # Tails
    i = i + 1
    X = gh(N, g, h)
    gi[i] = gimc(X)
    sk[i] = skewgm_mcf(X)
  }
}
for (nu in seq(2,20,by=1)) {
  i = i + 1
  X = rt(N, df=nu)
  gi[i] = gimc(X)
  sk[i] = skewgm_mcf(X)
}

# GF vs GMCF & BMCF ----
icol = as.numeric(factor(dft$Dataset)) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(dft$Dataset))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
png(file = '../article/fig_G_vs_GMC.png', 
    width = 2400, height=2400)

par(
  mfrow = c(2, 2),
  mar = mar,
  pty = pty,
  mgp = mgp,
  tcl = tcl,
  lwd = 4,
  cex = 4.5,
  xaxs = 'i',
  yaxs = 'i'
)

x = dft$gini
y = dft$gimc
ux = dft$u_gini
uy = dft$u_gimc
plot(x, y, pch=pch,
     xlim = c(0.1,0.7), xlab = expression(G[F]),
     ylim = c(0.1,0.7), ylab = expression(G[MCF]),
     col = cols[icol]
)
grid()
# segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
abline(a=0,b=1)
pch1 = unique(as.numeric(factor(dft$Dataset)))
sel1 = pch1 <= length(cols)
pch1[sel1] = 16
pch1[!sel1]= 17
legend(
  'bottomright', bty = 'o', box.col = 'white', ncol = 1,
  legend = unique(dft$Dataset),
  col = unique(cols[icol]),
  pch = pch1,
  cex=0.8
)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

plot(ux, uy, pch=pch,
     xlim = c(0.0,0.035), xlab = expression(u(G[F])),
     ylim = c(0.0,0.035), ylab = expression(u(G[MCF])),
     col = cols[icol]
)
grid()
# segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
abline(a=0,b=1)
abline(a=0,b=2,lty=2)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)


x = dft$gimc
y = dft$skewgm_mcf
ux = dft$u_gimc
uy = dft$u_skewgm_mcf
plot(x, y, pch=pch, type ='n',
     xlim = c(0.4,0.7), xlab = expression(G[MCF]),
     ylim = c(0.0,0.8), ylab = expression(beta[MCF]),
     col = cols[icol]
)
grid()
# segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
io = order(gi)
lines(gi[io],sk[io],col='gray50',lty=3,lwd=4)
points(x, y, pch=pch, col = cols[icol])
box()
mtext(
  text = '(c)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

plot(ux, uy, pch=pch,
     xlim = c(0.0,0.03), xlab = expression(u(G[MCF])),
     ylim = c(0.0,0.12), ylab = expression(u(beta[MCF])),
     col = cols[icol]
)
grid()
abline(a=0,b=1)
abline(a=0,b=2,lty=2)
abline(a=0,b=5,lty=2)
box()
mtext(
  text = '(d)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()
  
png(file = '../article/fig_GMCF_vs_Kurt_Skew.png', 
    width = 2400, height=1200)

par(
  mfrow = c(1, 2),
  mar = mar,
  pty = pty,
  mgp = mgp,
  tcl = tcl,
  lwd = 4,
  cex = 4.5,
  xaxs = 'i',
  yaxs = 'i'
)

x = dft$gimc
y = abs(dft$skewgm)
ux = dft$u_gini
uy = dft$u_skewgm
plot(x, y, pch=pch,
     xlim = c(0.4,0.7), xlab = expression(G[MCF]),
     ylim = c(-0.02,0.6), ylab = expression(abs(beta[GM])),
     col = cols[icol]
)
grid()
# segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
pch1 = unique(as.numeric(factor(dft$Dataset)))
sel1 = pch1 <= length(cols)
pch1[sel1] = 16
pch1[!sel1]= 17
legend(
  'topleft', bty = 'o', box.col = 'white', ncol = 1,
  legend = unique(dft$Dataset),
  col = unique(cols[icol]),
  pch = pch1,
  cex=0.8
)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

y =  dft$kurtcs
uy = dft$u_kurtcs
plot(x, y, pch=pch,
     xlim = c(0.4,0.7), xlab = expression(G[MCF]),
     ylim = c(-0.5,6.0), ylab = expression(kappa[CS]),
     col = cols[icol]
)
grid()
# segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()

resultsTab = file.path('..', 'results', 'tables', 'allStats_modeNew.csv')
D = read.csv(resultsTab)

# Ranking -----
png(file='../article/fig_GMCF_vs_Rank.png', 
    width=2400, height=1200)
mar  = c(3, 5, 1.1, 1)
par(mfrow = c(1, 2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

eps2 = 0.5

sel = D$ctd == -1 & !D$out & substr(D$Dataset,1,4) != 'Ref_'
setsNames = unique(D$Dataset[sel])
icol = as.numeric(factor(D$Dataset[sel]))
d = D[sel,]
d[['rmue']] = d$mue
for(s in setsNames) {
  s1 = which(d$Dataset == s)
  d$rmue[s1] = rank(d$mue[s1])
}

jcol= rep('gray80',length(sel))
s2 = pnorm(eps2,d$gimc,d$u_gimc) < 0.95 #d$gini >= 0.5
jcol[s2] = cols[2]

plot(d$rmue, icol, pch=19, col=jcol,
     xaxt = 'n', xlim = c(0.8, 10.2), xlab = 'rank(MUE)',
     yaxt = 'n', ylab ='', ylim = c(0.5,8.5))
axis(side = 1, at =1:10, gap.axis = 1/4) # Last arg enforces plot of "10"
mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
legend(
  6.2, 8.8,
  bty = 'n', cex = 0.75,
  legend = c(
    paste0('G < ',eps2), 
    paste0('G > ',eps2)
  ),
  pch = 19,
  col = c('gray80',cols[2]),
  y.intersp = 0.9
)
box()

# Remove global outliers ####
sel = D$ctd == -1 & D$out & substr(D$Dataset,1,4) != 'Ref_'
icol = as.numeric(factor(D$Dataset[sel]))
setsNames = unique(D$Dataset[sel])
d = D[sel,]
d[['rmue']] = d$mue
for(s in setsNames) {
  s1 = which(d$Dataset == s)
  d$rmue[s1] = rank(d$mue[s1])
}

jcol= rep('gray80',length(sel))
s2 = pnorm(eps2,d$gimc,d$u_gimc) < 0.95 #d$gini >= 0.5
jcol[s2] = cols[2]

plot(d$rmue, icol, pch=19, col=jcol,
     xaxt = 'n', xlim = c(0.8, 10.2), xlab = 'rank(MUE)',
     yaxt = 'n', ylab ='', ylim = c(0.5,8.5))
axis(side = 1, at =1:10, gap.axis = 1/4)
mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
legend(
  5.2, 8.5,
  bty = 'n', cex = 0.9,
  title = 'Removed outliers',
  legend = '',
  pch = -1
)
box()


dev.off()

#####################################


x = dft$gini
y = dft$gimc_pm
plot(x, y,pch=16,
     xlim = c(0.1,0.7),
     ylim = c(0.4,0.7), 
     col = cols[icol]
)
abline(a=0,b=1)
ux = dft$u_gini
uy = dft$u_gimc_pm
segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])

x = dft$gimc
y = dft$gimc_pm
plot(x, y,pch=16,
     xlim = c(0.45,0.55),
     ylim = c(0.4,0.7), 
     col = cols[icol]
)
abline(a=0,b=1)
ux = dft$u_gimc
uy = dft$u_gimc_pm
segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])

plot(ux, uy, pch=16,
     xlim = c(0,0.04),
     ylim = c(0,0.08), 
     col = cols[icol]
)
abline(a=0,b=1)
abline(a=0,b=3)

























# Plot ----

mc = TRUE


if(mc) {
  png(file = file.path('..', 'article', 'fig_G_vs_Skew_mcf.png'), 
      width=2400, height=2400)
} else {
  png(file = file.path('..', 'article', 'fig_G_vs_Skew_f.png'), 
      width=2400, height=2400)
}

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
    if(mc)
      X = X - hrmode(X)
    gi[i] = ErrViewLib::gini(X)
    sk[i] = skewgm_f(X)
  }
}

plot(
  gi, sk, pch=16, col = cols[2],
  xlim = c(0.4,0.7), xlab = expression(G[F]),
  ylim = c(0.2,0.8), ylab = expression(beta[GMF]),
  main = paste0('Ref / N=',N)
)
grid()

i=0
gi = sk = c()
for (nu in seq(2,20,by=1)) {
  i = i + 1
  X = rt(N, df=nu)
  if(mc)
    X = X - hrmode(X)
  gi[i] = ErrViewLib::gini(X)
  sk[i] = skewgm_f(X)
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
    if(mc)
      X = X - hrmode(X)
    gi[i] = ErrViewLib::gini(X)
    sk[i] = skewgm_f(X)
  }
}
plot(
  gi, sk, pch=16, col = cols[2],
  xlim = c(0.4,0.7), xlab = expression(G[F]),
  ylim = c(0.2,0.8), ylab = expression(beta[GMF]),
  main = paste0('Ref / N=',N)
)
grid()

i=0
gi = sk = c()
for (nu in seq(2,20,by=1)) {
  i = i + 1
  X = rt(N, df=nu)
  if(mc)
    X = X - hrmode(X)
  gi[i] = ErrViewLib::gini(X)
  sk[i] = skewgm_f(X)
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
icol = as.numeric(factor(dft$Dataset)) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(dft$Dataset))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17

if(mc) {
  x = dft$gimc
  y = dft$skewgm_mcf
  ux = dft$u_gimc
  uy = dft$u_skewgm_mcf
} else {
  x = dft$gini
  y = dft$skewgm_f
  ux = dft$u_gini
  uy = dft$u_skewgm_f
}
plot(
  x, y,
  col=cols[icol], 
  pch=pch,
  xlim = c(0.4,0.7), xlab = expression(G[F]),
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

plot(
  ux,uy,
  pch=16,
  col = cols[icol],
  xlim = c(0.,0.03), xlab = expression(u(G[F])),
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