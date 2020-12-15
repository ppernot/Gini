source('./0-Setup.R')

# Fig. 1 ----
cat('Fig. 1\n')
frow = 1; fcol = 2
png(file= file.path('..','results','figs','Fig_1.png'), 
    width= fcol * reso, height=frow * reso)
par(
  mfrow = c(frow, fcol),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

## CDF ----
c = curve(
  2 * pnorm(x) - 1,
  from = 0,
  to = 3,
  xaxs = 'i', xlab = '|Error|',
  yaxs = 'i', ylab = 'Probability',
  pty = 's',
  main = 'Cumulative Distribution'
)
polygon(
  c(c$x, 0),
  c(c$y, 1),
  density = 2,
  angle = 135,
  col = cols[3],
  border = NA
)

p = 0.9
qp = qnorm(p+(1-p)/2)
abline(h = p, lty = 2, col = cols[2])
abline(v = qp,lty = 2, col = cols[2])
mtext(
  "p'",
  side = 2,
  at = p,
  las = 1,
  col = cols[2],
  adj = 1.8,
  cex = cex
)
mtext(
  expression(Q[p~"'"]),
  side = 1,
  at = qp,
  las = 1,
  col = cols[2],
  padj = 0.3,
  cex = cex
)

sel = c$y <= p
polygon(
  c(0, c$x[sel], 0),
  c(p, c$y[sel], p),
  density = NULL,
  angle = 45,
  col = cols_tr[6],
  border = NA
)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

## Lorenz ----
lc = ineq::Lc(abs(rnorm(10000, 0, 1)))
plot(lc,
     lty =1, lwd = lwd,
     col = cols[1])
abline(v=p,col = cols[2], lty = 2)
mtext(
  "p'",
  side = 1,
  at = p,
  las = 0,
  col = cols[2],
  padj = 0.3,
  cex = cex
)
lpp = lc$L[which(lc$p>=p)[1]]
abline(h=lpp,col = cols[2], lty = 2)
mtext(
  expression(L[p~"'"]),
  side = 2,
  at = lpp,
  las = 2,
  col = cols[2],
  adj = 2,
  cex = cex
)
polygon(
  c(lc$p, 0),
  c(lc$L, 0),
  density = 2,
  angle = 90,
  col = cols[4],
  border = NA
)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()

# Fig. 2 ----
cat('Fig. 2\n')

## Generate plot data ----
N = 1e5
bl = 5
db = 0.1

gi = ku = bi = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b + rnorm(N)
  bi[i] = b
  gi[i] = ErrViewLib::gini(X)
  ku[i] = ErrViewLib::kurtcs(X)
}

gi0 = ku0 = bi0 = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b + runif(N,-1,1)
  bi0[i] = b
  gi0[i] = ErrViewLib::gini(X)
  ku0[i] = ErrViewLib::kurtcs(X)
}

nu = 2
gi1 = ku1 = bi1 = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b + rt(N, df = nu)
  bi1[i] = b
  gi1[i] = ErrViewLib::gini(X)
  ku1[i] = ErrViewLib::kurtcs(X)
}

mu = 1
si = 0.5
gi2 = ku2 = bi2 = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b + rlnorm(N,mu,si)
  bi2[i] = b
  gi2[i] = ErrViewLib::gini(X)
  ku2[i] = ErrViewLib::kurtcs(X)
}

m = exp(mu-si^2)
gi2mc = ku2mc = bi2mc = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b - m + rlnorm(N,mu,si)
  bi2mc[i] = b
  gi2mc[i] = ErrViewLib::gini(X)
  ku2mc[i] = ErrViewLib::kurtcs(X)
}

g = 1
h = 0
gi3 = ku3 = bi3 = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b + gh(N,g,h)
  bi3[i] = b
  gi3[i] = ErrViewLib::gini(X)
  ku3[i] = ErrViewLib::kurtcs(X)
}

m = hrmode(gh(1e6,g,0))
gi3mc = ku3mc = bi3mc = c()
i = 0
for (b in seq(-bl, bl, by = db)) {
  i = i + 1
  X = b - m + gh(N,g,h)
  bi3mc[i] = b
  gi3mc[i] = ErrViewLib::gini(X)
  ku3mc[i] = ErrViewLib::kurtcs(X)
}

## Plot ----

frow = 1; fcol = 2
png(file= file.path('..','results','figs','Fig_2.png'), 
    width= fcol * reso, height=frow * reso)

par(
  mfrow = c(frow, fcol),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

plot(
  bi, gi,
  col = cols[2],
  type = 'l',
  lty = 1,
  ylim = c(0.0, 0.7),
  ylab = expression(G[F]),
  xlim = c(-5, 5),
  xlab = 'Bias',
  xaxs = 'i',
  yaxs = 'i',
  main = ''
)
grid()
lines(bi0,
      gi0,
      lty = 5,
      col = cols[5])
lines(bi1,
      gi1,
      lty = 2,
      col = cols[6])
lines(bi2,
      gi2,
      lty = 3,
      col = cols[3])
lines(bi3,
      gi3,
      lty = 4,
      col = cols[4])
legend(
  'topright', cex=0.75,
  bty = 'o',
  box.col = 'white',
  legend = c('Normal', 'Student', 'logNormal','g-and-h','Uniform'),
  col = cols[c(2, 6, 3, 4, 5)],
  lty = 1:5
)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

plot(
  bi,
  gi,
  col = cols[2],
  type = 'l',
  lty = 1,
  ylim = c(0.0, 0.7),
  ylab = expression(G[F]),
  xlim = c(-5, 5),
  xlab = 'Bias',
  xaxs = 'i',
  yaxs = 'i',
  main = ''
)
grid()
lines(bi0,
      gi0,
      lty = 5,
      col = cols[5])
lines(bi1,
      gi1,
      lty = 2,
      col = cols[6])
lines(bi2mc,
      gi2mc,
      lty = 3,
      col = cols[3])
lines(bi3mc,
      gi3mc,
      lty = 4,
      col = cols[4])
legend(
  'topright', cex=0.75,
  bty = 'o',
  box.col = 'white',
  title = 'Mode-centered',
  legend = ''
)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()

# Fig. 3 ----
cat('Fig. 3\n')

## Read plot data ----

D = read.csv(resultsTab)
sets = D$ctd == -1 & # No linear correction
       ! D$out       # No outliers removed
d = D[sets,]

## Plot ----

frow = 1; fcol = 2
png(file= file.path('..','results','figs','Fig_3.png'), 
    width= fcol * reso, height=frow * reso)
par(
  mfrow = c(frow, fcol),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

x  = d$hrmode
y  = d$bmax
ux = d$u_hrmode
uy = d$u_bmax
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

x  = d$gimc
y  = d$gmax
ux = d$u_gimc
uy = d$u_gmax
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

# Fig. 4 ----
cat('Fig. 4\n')

## Read / Generate plot data ----

D = read.csv(resultsTab)
sets = D$ctd == -1 & # No linear correction
      !D$out         # No outliers removed
d = D[sets,]

### Calculate data for subplot (c)  ----
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

## Plot ----
frow = 2; fcol = 2
png(file = file.path('..','results','figs','Fig_4.png'), 
    width = fcol * reso, height = frow * reso)
par(
  mfrow = c(frow, fcol),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

# Color/pch scheme
icol = as.numeric(factor(d$Dataset)) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(d$Dataset))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17

### GMCF vs. G (a) ----

x = d$gini
y = d$gimc
ux = d$u_gini
uy = d$u_gimc

plot(x, y, pch=pch,
     xlim = c(0.1,0.7), xlab = expression(G[F]),
     ylim = c(0.1,0.7), ylab = expression(G[MCF]),
     col = cols[icol]
)
grid()
abline(a=0,b=1)

# pch scheme for legend
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

### u(GMCF) vs. u(G) (b) ----

plot(ux, uy, pch=pch,
     xlim = c(0.0,0.035), xlab = expression(u(G[F])),
     ylim = c(0.0,0.035), ylab = expression(u(G[MCF])),
     col = cols[icol]
)
grid()
abline(a=0,b=1)
abline(a=0,b=2,lty=2)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

### betaGM vs GMCF (c) ----

x  = d$gimc
y  = d$skewgm_mcf
ux = d$u_gimc
uy = d$u_skewgm_mcf

plot(x, y, pch=pch, type ='n',
     xlim = c(0.4,0.7), xlab = expression(G[MCF]),
     ylim = c(0.0,0.8), ylab = expression(beta[MCF]),
     col = cols[icol]
)
grid()

# Smooth guide line
reg = lm(sk ~ 1 + gi + I(gi^2))
pr = predict(reg)
io = order(pr)
lines(gi[io],pr[io],col='gray50',lty=2,lwd=5)

points(x, y, pch=pch, col = cols[icol])

box()
mtext(
  text = '(c)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

### u(betaGM) vs u(GMCF) (d) ----

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

# Fig. 5 ----
cat('Fig. 5\n')

## Read plot data ----
D = read.csv(resultsTab)
sets = D$ctd == -1 & # No linear correction
  ! D$out       # No outliers removed
d = D[sets,]

## Plot ----

frow = 1; fcol = 2
png(file= file.path('..','results','figs','Fig_5.png'), 
    width= fcol * reso, height=frow * reso)
par(
  mfrow = c(frow, fcol),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

x = d$gimc
y = abs(d$skewgm)
ux = d$u_gini
uy = d$u_skewgm

# Color/pch scheme
icol = as.numeric(factor(d$Dataset)) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(d$Dataset))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17

### (a) ----
plot(x, y, pch=pch,
     xlim = c(0.4,0.7), xlab = expression(G[MCF]),
     ylim = c(-0.02,0.6), ylab = expression(abs(beta[GM])),
     col = cols[icol]
)
grid()

# pch scheme for legend
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

### (b) ----

y =  d$kurtcs
uy = d$u_kurtcs
plot(x, y, pch=pch,
     xlim = c(0.4,0.7), xlab = expression(G[MCF]),
     ylim = c(-0.5,6.0), ylab = expression(kappa[CS]),
     col = cols[icol]
)
grid()
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()

# Fig. 6 ----
cat('Fig. 6\n')

## Read plot data ----
D = read.csv(resultsTab)

## Plot ----

frow = 1; fcol = 2
png(file= file.path('..','results','figs','Fig_6.png'), 
    width= fcol * reso, height=frow * reso)
par(
  mfrow = c(frow, fcol),
  pty = pty,
  mar = c(3, 5, 1.1, 1),
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

G_thresh = 0.5

### (a) ----

sets = D$ctd == -1 & # No linear correction
  ! D$out       # No outliers removed
d = D[sets,]

# MUE-Ranking
setsNames = unique(d$Dataset)
d[['rmue']] = d$mue
for(s in setsNames) {
  s1 = which(d$Dataset == s)
  d$rmue[s1] = rank(d$mue[s1])
}

# Colors
icol = as.numeric(factor(d$Dataset))
jcol= rep('gray80',length(sets))
s2 = pnorm(G_thresh,d$gimc,d$u_gimc) < 0.95 # GMCF >= 0.5
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
    expression(G[MCF] < 0.5), 
    expression(G[MCF] >= 0.5)
  ),
  pch = 19,
  col = c('gray80',cols[2]),
  y.intersp = 0.9
)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

### (b) ----
# Remove global outliers 
sets = D$ctd == -1 & 
  D$out 
d = D[sets,]

# MUE-Ranking
setsNames = unique(d$Dataset)
d[['rmue']] = d$mue
for(s in setsNames) {
  s1 = which(d$Dataset == s)
  d$rmue[s1] = rank(d$mue[s1])
}

# Colors
icol = as.numeric(factor(d$Dataset))
jcol = rep('gray80',length(sets))
s2 = pnorm(eps2,d$gimc,d$u_gimc) < 0.95 #d$gini >= 0.5
jcol[s2] = cols[2]

plot(d$rmue, icol, pch=19, col=jcol,
     xaxt = 'n', xlim = c(0.8, 10.2), xlab = 'rank(MUE)',
     yaxt = 'n', ylab ='', ylim = c(0.5,8.5))
axis(side = 1, at =1:10, gap.axis = 1/4)
mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
legend(
  3.5, 8.5,
  bty = 'n', cex = 0.9,
  title = 'Global outliers removed',
  legend = '',
  pch = -1
)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

dev.off()

# Fig. 7 ----
cat('Fig. 7\n')

frow = 1; fcol = 1
png(
  file = file.path('..', 'results', 'figs', 'Fig_7.png'),
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


mue = 1

mu0 = 0
s0  = sqrt(pi / 2)

mu1 = mue
s1  = 0.1

# Proba to exceed MUE for absolute errors errors
(p1 = 1-2*(pnorm(mue,mu0,s0)-0.5)) # 0.42 ; accounts for folding
(p2 = 1-pnorm(mue,mu1,s1)) # 0.5 ; no folding necessary

curve(dnorm(x, mu0, s0),
      from = -2.75,
      to = 2.75,
      lwd = lwd,
      col = cols[2],
      xlab = "Error",
      xlim = c(-2.75,2.75),
      xaxs = 'i',
      ylab = "Probability Density Function",
      ylim = c(0,4.5),
      yaxs = 'i')
grid()
abline(v = 0, lwd = 1.5)
abline(v = mue, col = cols[4], lty = 2)
mtext('MUE',
      side = 3,
      at = mue,
      cex = cex,
      col = cols[4])
curve(
  dnorm(x, mu1, s1),
  from = 0,
  to = 2,
  lwd = lwd,
  col = cols[6],
  add = TRUE
)
legend(
  'topleft', bty = 'o', box.col = 'white',
  legend = c('N(0,1.25)','N(1,0.1)'),
  lwd = lwd, lty = 1,
  col = cols[c(2,6)]
)
box()

dev.off()
