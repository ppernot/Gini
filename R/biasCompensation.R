
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


png(file = '../article/fig_biasCompensation.png', width=2400, height=1200)
cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1.5, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 5
cex  = 4.5

par(
  mfrow = c(1,2),
  mar = mar,
  pty = pty,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

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