gimc = function(X, index = 1:length(X), ...) {
  X = X[index]
  ErrViewLib::gini(abs(X - hrmode(X)))
}
hrmode = function(X, index = 1:length(X), ...) {
  genefilter::half.range.mode(X[index])
}
bmax = function(X, index = 1:length(X), ...) {
  X = X[index]
  m = hrmode(X)
  gfun = function (b,X)
    ErrViewLib::gini(abs(X - b))
  o = optim(
    par   = m,
    fn    = gfun,
    # method = "Brent",
    lower = min(X),
    upper = max(X),
    control = list(fnscale=-1),
    X     = X
  )
  return(o$par)
}
gmax = function(X, index = 1:length(X), ...) {
  X = X[index]
  m = hrmode(X)
  gfun = function (b,X)
    ErrViewLib::gini(abs(X - b))
  o = optim(
    par   = m,
    fn    = gfun,
    # method = "Brent",
    lower = min(X),
    upper = max(X),
    control = list(fnscale=-1),
    X     = X
  )
  return(ErrViewLib::gini(abs(X - o$par)))
}

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

rm('dftest')
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
    'gimc',
    'gmax',
    # 'gimc',
    'hrmode',
    'bmax'
  )
  
  bs = ErrViewLib::estBS1(
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
  if(!exists('dftest'))
    dftest = df
  else
    dftest = rbind(dftest,df)
}

# Plots ----

# cols = rev(inlmisc::GetColors(8))[1:7]
# pty  = 's'
# mar  = c(3, 3, 1.1, 1)
# mgp  = c(2, .75, 0)
# tcl  = -0.5
# lwd  = 1 #4
# cex  = 1 #4.5
# 
# # ColorScheme
# icol = as.numeric(factor(dftest$Dataset)) %% length(cols)
# icol[icol==0] = length(cols)
# pch = as.numeric(factor(dftest$Dataset))
# sel1 = pch<=length(cols)
# pch[sel1] = 16
# pch[!sel1]= 17
# 
# par(mfrow = c(2,2),
#     pty = pty,
#     mar = mar,
#     mgp = mgp,
#     tcl = tcl,
#     lwd = lwd,
#     cex = cex)
# 
# x  = dftest$gimc
# y  = dftest$gmax
# ux = dftest$u_gimc
# uy = dftest$u_gmax
# plot(x, y, pch=pch,
#      xlim = c(0.4,0.7), xlab = expression(G[MCF]),
#      ylim = c(0.4,0.7), ylab = expression(G[max]),
#      col = cols[icol]
# )
# grid()
# # segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# # segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
# abline(a=0,b=1)
# box()
# 
# plot(ux, uy, pch=pch,
#      xlim = c(0.0,0.035), xlab = expression(u(G[MCF])),
#      ylim = c(0.0,0.035), ylab = expression(u(G[max])),
#      col = cols[icol]
# )
# grid()
# abline(a=0,b=1)
# abline(a=0,b=1/2,lty=2)
# 
# 
# x  = dftest$hrmode
# y  = dftest$bmax
# ux = dftest$u_hrmode
# uy = dftest$u_bmax
# plot(x, y, pch=pch,
#      xlab = expression(mode),
#      ylab = expression(b[max]),
#      col = cols[icol]
# )
# grid()
# # segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
# # segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
# abline(a=0,b=1)
# box()
# 
# plot(ux, uy, pch=pch,
#      xlab = expression(u(mode)),
#      ylab = expression(u(b[max])),
#      col = cols[icol]
# )
# grid()
# abline(a=0,b=1)
# abline(a=0,b=1/2,lty=2)

## z-score ----
png(file='../article/fig_hrmode_vs_bmax.png', width=2400, height=1200)
cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1.1, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 4.5

par(mfrow = c(1,2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

x  = dftest$hrmode
y  = dftest$bmax
ux = dftest$u_hrmode
uy = dftest$u_bmax
z = (x-y)/sqrt(ux^2+uy^2)

hist(z, col = cols[6],
     xlab = expression(z[b]),
     main = '')
abline(v=-2,lty=2)
# box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

x  = dftest$gimc
y  = dftest$gmax
ux = dftest$u_gimc
uy = dftest$u_gmax
z = (x-y)/sqrt(ux^2+uy^2)
hist(z, col = cols[3],
     xlab = expression(z[G]),
     main = '')
abline(v=-2,lty=2)
# box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)
dev.off()

# Check ----
gfun = function (b,X)
  ErrViewLib::gini(abs(X -b ))
     
set = 'ZAS2019'
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

X = Errors[,'SLATM_L2']

x = seq(-1,1.5,by=0.05)
y = c()
for (i in 1:length(x)) 
  y[i] = ErrViewLib::gini(abs(X - x[i])) 
plot(x,y)
abline(v=hrmode(X))
abline(v = bmax(X),lty=2)
