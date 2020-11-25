resultsTab = file.path('..', 'results', 'tables', 'allStats.csv')

D = read.csv(resultsTab)

D[['z']]       = abs(D$mse / D$rmsd)
D[['logit_W']]  = boot::logit(D$W)

sel = D$Dataset != 'DAS2019' &
  substr(D$Dataset,1,7) != 'VER2019'
D = D[sel,]

# Correlations ####
pdf(file.path('..','results', 'figs', 'correls.pdf'),
    height = 15,
    width = 15)
par(mfrow = c(3, 3), pty = 's', cex=1.5)

for (stat1 in c('gini', 'lasym','pietra'))
  for (stat2 in c('z', 'skew', 'kurt', 'logit_W')) {
    xlim = range(D[[stat1]])
    ylim = range(D[[stat2]])
    sel = D$ctd < 0 & !D$out
    x = D[[stat1]][sel]
    y = D[[stat2]][sel]
    icol = factor(D$Dataset[sel])
    plot(
      x,
      y,
      pch = 16,
      xlab = stat1,
      xlim = xlim,
      ylab = stat2,
      ylim = ylim,
      col = icol,
      main = paste0('Raw errors / ', signif(cor(x, y), 2))
    )
    grid()
    abline(lm(y ~ x), col = 2)
    # legend('topright', title = paste(stat1, stat2), legend = '')
    box()

    sel = D$ctd == 0 & !D$out
    x = D[[stat1]][sel]
    y = D[[stat2]][sel]
    icol = factor(D$Dataset[sel])
    plot(
      x,
      y,
      pch = 16,
      xlab = stat1,
      xlim = xlim,
      ylab = stat2,
      ylim = ylim,
      col = icol,
      main = paste0('c-errors / ', signif(cor(x, y), 2))
    )
    grid()
    abline(lm(y ~ x), col = 2)
    box()

    sel = D$ctd == 1 & !D$out
    x = D[[stat1]][sel]
    y = D[[stat2]][sel]
    icol = factor(D$Dataset[sel])
    plot(
      x,
      y,
      pch = 16,
      xlab = stat1,
      xlim = xlim,
      ylab = stat2,
      ylim = ylim,
      col = icol,
      main = paste0('lc-errors / ', signif(cor(x, y), 2))
    )
    grid()
    abline(lm(y ~ x), col = 2)
    box()
  }
dev.off()

# Gini vs. CV ####
png(file='../article/corr_Gini_CV.png', width=1200, height=1200)

cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1.2, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 4

par(mfrow = c(1, 1),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

stat1 = 'gini'
stat2 = 'z'
xlim = c(0,1)
ylim = c(0,4)
sel = D$ctd < 0 & !D$out
x = D[[stat1]][sel]
y = D[[stat2]][sel]
s = D[['rmsd']][sel]
w = D[['W']][sel]
icol = factor(D$Dataset[sel])
plot(
  x,
  1/y/sqrt(pi),
  pch = 1,
  cex = cex + (x-0.7)*5,
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = expression(CV / sqrt(pi)),
  ylim = ylim, yaxs = 'i',
  col = cols[6],
  main = 'Raw errors'#paste0('Raw errors / ', signif(cor(x, y), 2))
)
grid()
abline(a=0,b=1,lty=2,col='gray70',lwd=2*lwd)
abline(v=0.4147,lty=2,col='gray70',lwd=2*lwd)
x = c()
y = c()
mu = seq(0,10,by=0.05)
for(i in 1:length(mu)) {
  X = abs(rnorm(50000,mu[i],1))
  x[i] = ineq::Gini(X)
  y[i] = 1/(mu[i]*sqrt(pi))
}
lines(x,y, lty = 1, col = cols[3], lwd=2*lwd)
box()
# mtext(
#   text = '(a)',
#   side = 3,
#   adj = 1,
#   cex = cex,
#   line = 0.3)

dev.off()

# Gini vs. Kurt ####
png(file='../article/corr_Gini_Kurt+W.png', width=2400, height=1200)

par(mfrow = c(1, 2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

stat1 = 'gini'
stat2 = 'kurt'
xlim = c(0,1)
ylim = c(0,20)
sel = D$ctd == 0 & !D$out
x = D[[stat1]][sel]
y = D[[stat2]][sel]
icol = factor(D$Dataset[sel])
plot(
  x,
  y,
  pch = 1,
  # cex = cex + (x-0.7)*5,
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = 'Kurtosis',
  ylim = ylim, yaxs = 'i',
  col = cols[6],
  main = 'Centered errors'#paste0('Raw errors / ', signif(cor(x, y), 2))
)
grid()
abline(lm(y~x),lty=2,col='gray70',lwd=2*lwd)
# abline(a=0,b=1,lty=2,col='gray70',lwd=2*lwd)
# abline(v=0.4147,lty=2,col='gray70',lwd=2*lwd)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

stat1 = 'gini'
stat2 = 'W'
xlim = c(0,1)
ylim = c(-1,6)
sel = D$ctd == 0 & !D$out
x = D[[stat1]][sel]
y = D[[stat2]][sel]
icol = factor(D$Dataset[sel])
plot(
  x,
  boot::logit(y),
  pch = 1,
  # cex = cex + (x-0.7)*5,
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = 'logit(W)',
  ylim = ylim, yaxs = 'i',
  col = cols[6],
  main = 'Centered errors'#paste0('Raw errors / ', signif(cor(x, y), 2))
)
grid()
abline(lm(boot::logit(y)~x),lty=2,col='gray70',lwd=2*lwd)
# abline(a=0,b=1,lty=2,col='gray70',lwd=2*lwd)
# abline(v=0.4147,lty=2,col='gray70',lwd=2*lwd)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()

# Gini vs. Pietra ####
png(file='../article/corr_Gini_Pietra.png', width=1200, height=1200)

par(mfrow = c(1, 1),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

# for(dat in unique(D$Dataset)) {
  stat1 = 'gini'
  stat2 = 'pietra'
  xlim = c(0,0.8)
  ylim = c(0,0.8)
  sel = D$ctd == -1 & !D$out #& D$Dataset == dat
  x = D[[stat1]][sel]
  ux = D[[paste0('u_',stat1)]][sel]
  y = D[[stat2]][sel]
  uy = D[[paste0('u_',stat2)]][sel]
  icol = as.numeric(factor(D$Dataset[sel]))
  cols2 = c(cols,cols)
  pch = c(rep(1,length(cols)),rep(2,length(cols)))
  plot(
    x,
    y,
    pch  = pch,
    xlab = 'Gini',
    xlim = xlim, xaxs = 'i',
    ylab = 'Pietra',
    ylim = ylim, yaxs = 'i',
    col  = cols2[icol],
    main = 'Raw errors'
  )
  # segments(x,y-2*uy,x,y+2*uy,col=cols2[icol])
  # segments(x-2*ux,y,x+2*ux,y,col=cols2[icol])
  grid()
  abline(lm(y~0+x),lty=2,col='gray70',lwd=2*lwd)
  points(
    x,
    y,
    pch  = pch,
    col  = cols2[icol]
  )
  legend(
    'topleft', bty='n',
    legend = unique(D$Dataset[sel]),
    pch = pch,
    col = unique(cols2[icol]),
    cex=0.75
  )
  box()
#   mtext(
#     text = '(a)',
#     side = 3,
#     adj  = 1,
#     cex  = cex,
#     line = 0.3)
#
# }

dev.off()

stop()

