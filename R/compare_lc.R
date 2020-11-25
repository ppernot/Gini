lcp = function(X, index = 1:length(X), p=0.95, ...) {
  # Lorenz curve
  X = sort(abs(X[index]))
  pr = (1:length(X)) / length(X)
  lc = cumsum(X)/sum(X)
  # linear interpolation
  iup = which(pr >= p)[1]
  ilo = iup - 1
  Lcp = lc[ilo]+ (p-pr[ilo])*
    (lc[iup]-lc[ilo])/(pr[iup]-pr[ilo])
  return(1-Lcp)
}

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

cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 2.5

# Shape ####
png(file='../article/fig_compare_lc.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

curve(
  dunif(x, 1,2),
  from = 0,
  to = 8,
  col = cols[6], lwd = lwd,
  xlab = 'Errors (a.u.)',
  ylab = 'PDF',
  yaxs = 'i', ylim = c(0,1.1)
)
grid(col='gray70')
curve(
  dnorm(x, 5, 1),
  from = 0,
  to = 10,
  col = cols[2], lwd = lwd,
  add = TRUE
)
box()

X1 = abs(runif(10000,1,2))
print(sd(X1)/mean(X1))
X2 = abs(rnorm(10000,5))
print(sd(X2)/mean(X2))

plot(ineq::Lc(X1),col=cols[6], lwd = lwd)
grid(col='gray70')
lines(ineq::Lc(X2), lty=2, col=cols[2],lwd = lwd)
legend(
  'topleft', bty = 'n',
  legend = c('Unif(1,2)','Norm(5,1)'),
  lty = c(1,2), lwd = lwd,
  col = cols[c(6,2)]
)
box()
dev.off()

# Bias ####
png(file='../article/fig_bias_lc.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)


mu = c(0,1,2,3,4,5,6)


curve(
  dnorm(x, mu[1],1),
  from = -3,
  to = 9,
  col = cols[2], lwd = lwd,
  xlab = 'Errors (a.u.)',
  ylab = 'PDF',
  yaxs = 'i', ylim = c(0,0.42)
)
grid(col = 'gray70')
abline(v=0,lty=2)
for (i in 2:length(mu))
  curve(
    dnorm(x, mu[i], 1),
    from = -3,
    to = 9,
    col = cols[1+i], lwd = lwd,
    add = TRUE
  )
box()

plot(ineq::Lc(abs(rnorm(10000, mu[1], 1))),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
for (i in 2:length(mu))
  lines(ineq::Lc(abs(rnorm(10000, mu[i], 1))),
        lty = 1, lwd = lwd,
        col = cols[i+1])
box()
legend(
  'topleft', bty = 'n',
  title = expression(mu),
  legend = mu,
  lty = 1, lwd = lwd,
  col = cols[2:7]
)
dev.off()

# GS ####
png(file='../article/fig_bias_tails_GS.png', width=1600, height=800)
lwd = 4
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

mu = seq(0,20,by=0.1)
S = G = CV = l95 = c()
for (i in 1:length(mu)) {
  X = abs(rnorm(10000, mu[i], 1))
  G[i] = ineq::Gini(X)
  S[i] = ineq::Lasym(X)
  CV[i] = 1/mu[i]
  l95[i] = lcp(X)
}

plot(
  mu, S,
  type = 'l',
  lty = 1,
  lwd = lwd,
  col = cols[6],
  xlab = expression(mu),
  xaxs = 'i',
  ylim = c(0, 1.1),
  yaxs = 'i',
  ylab = 'G, S'
)
grid(col='gray70')

lines(
  mu,CV/sqrt(pi),
  lty = 2,
  col = cols[3],
  lwd = lwd
)
lines(
  mu,G,
  lty = 1,
  lwd = lwd,
  col = cols[2],
)
lines(
  mu,l95,
  lty = 1,
  lwd = lwd,
  col = cols[4],
)
legend(
  'right',  bty = 'n',
  legend = c('G', 'S',
             expression(CV / sqrt(mu)),
             expression(l[95])),
  lty = c(1, 1, 2, 1),
  col = cols[c(2, 6, 3, 4)],
  lwd = lwd
)
box()

df = seq(1,20,by=1)
S = G = l95 = c()
for (i in 1:length(df)) {
  X = abs(rt(10000, df[i]))
  G[i] = ineq::Gini(X)
  S[i] = ineq::Lasym(X)
  l95[i] = lcp(X)
}
plot(
  df, S,
  type = 'l',
  lty = 1,
  lwd = lwd,
  col = cols[6],
  xaxs = 'i',
  xlab = expression(nu~', '~g^{-1}),
  ylim = c(0, 1.1),
  yaxs = 'i',
  ylab = 'G, S'
)
grid(col='gray70')

lines(
  df,G,
  lty = 1,
  lwd = lwd,
  col = cols[2],
)
lines(
  df,l95,
  lty = 1,
  lwd = lwd,
  col = cols[4],
)

g = seq(0,1,by=0.02)
S = G = l95 = c()
for (i in 1:length(g)) {
  X = abs(gh(100000,g[i],0))
  G[i] = ineq::Gini(X)
  S[i] = ineq::Lasym(X)
  l95[i] = lcp(X)
}
lines(
  1/g, S,
  lty = 2,
  lwd = lwd,
  col = cols[6],
)
lines(
  1/g,G,
  lty = 2,
  lwd = lwd,
  col = cols[2],
)
lines(
  1/g,l95,
  lty = 2,
  lwd = lwd,
  col = cols[4],
)

box()
dev.off()

# Tails ####
png(file='../article/fig_tails_lc.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)


df = c(1,2,4,8,100)

curve(
  dt(x, df[1]),
  from = -6,
  to = 6,
  col = cols[2], lwd = lwd,
  xlab = 'Errors (a.u.)',
  xlim = c(-5,5),
  xaxs = 'i',
  ylab = 'PDF',
  yaxs = 'i',
  ylim = c(0,0.42)
)
grid(col = 'gray70')
abline(v=0,lty=2)
for (i in 2:length(df))
  curve(
    dt(x, df[i]),
    from = -5,
    to = 5,
    col = cols[1+i], lwd = lwd,
    add = TRUE
  )
box()

plot(ineq::Lc(abs(rt(10000, df[1]))),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
for (i in 2:length(df))
  lines(ineq::Lc(abs(rt(10000, df[i]))),
        lty = 1, lwd = lwd,
        col = cols[i+1])
box()
legend(
  'topleft', bty = 'n',
  title = expression(nu),
  legend = df,
  lty = 1, lwd = lwd,
  col = cols[2:7]
)
dev.off()
for(i in 1:length(df)) {
  X = abs(rt(100000, df[i]))
  cat(' df = ',df[i],
      ' gini = ',gini(X),
      ' lasym = ',lasym(X),'\n')

}

# Tails-gh ####

sam = list()
i=0
dmax = 0
for (g in c(0,0.25,0.5,0.75,1)[1]) { # Asymmetry
  for (h in c(0,0.1,0.2,0.25)) { # Tails
    i=i+1
    X = gh(1000000,g,h)
    d = density(X)
    sam[[i]] = list(g=g,h=h,X=X,d=d)
    dmax = max(dmax,d$y)
  }
}

png(file='../article/fig_tailsgh_lc.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)


leg = c()
for(i in 1:length(sam)) {
  d = sam[[i]]$d
  if(i == 1) {
    plot(
      d,
      col = cols[2],
      lwd = lwd,
      xlab = 'Errors (a.u.)',
      xlim = c(-5,5),
      xaxs = 'i',
      ylab = 'PDF',
      yaxs = 'i',
      ylim = c(0,1.05*dmax),
      main=''
    )
    grid()
  } else {
    lines(
      d,
      col = cols[i+1],
      lwd = lwd
    )
  }
  leg[i] = paste0('g = ',sam[[i]]$g,', h = ',sam[[i]]$h)
}
box()


plot(ineq::Lc(abs(sam[[1]]$X)),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
for (i in 2:length(sam))
  lines(ineq::Lc(abs(sam[[i]]$X)),
        lty = 1, lwd = lwd,
        col = cols[i+1])
X = rnorm(10000)
lines(ineq::Lc(abs(X)),
      lty = 2, lwd = lwd,
      col = 'gray70')
box()

legend(
  'topleft', bty = 'n',
  cex=0.6,
  # title = 'ncp',
  legend = leg,
  lty = 1, lwd = lwd,
  col = cols[2:7]
)
dev.off()

# Asymmetry ####
sam = list()
i=0
dmax = 0
for (g in c(0,0.25,0.5,0.75,1)) { # Asymmetry
  for (h in c(0,0.2,0.4,0.6)[1]) { # Tails
    i=i+1
    X = 0 + gh(1000000,g,h)
    d = density(X,n=2048)
    gi = gini(X)
    p = pietra(X)
    la = lasym(X)
    sam[[i]] = list(g=g,h=h,X=X,d=d,gi=gi,p=p,la=la)
    dmax = max(dmax,d$y)
  }
}
for(i in 1:length(sam))
  cat(' g = ',sam[[i]]$g,
      ' gini = ',sam[[i]]$gi,
      ' lasym = ',sam[[i]]$la,'\n')

png(file='../article/fig_asym_lc.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)


leg = c()
for(i in 1:length(sam)) {
   d = sam[[i]]$d
   if(i == 1) {
     plot(
       d,
       col = cols[2],
       lwd = lwd,
       xlab = 'Errors (a.u.)',
       xlim = c(-5,5),
       xaxs = 'i',
       ylab = 'PDF',
       yaxs = 'i',
       ylim = c(0,1.05*dmax),
       main=''
     )
     grid()
   } else {
     lines(
       d,
       col = cols[i+1],
       lwd = lwd
     )
   }
   leg[i] = paste0('g = ',sam[[i]]$g,', h = ',sam[[i]]$h)
}
box()


plot(ineq::Lc(abs(sam[[1]]$X)),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
for (i in 2:length(sam))
  lines(ineq::Lc(abs(sam[[i]]$X)),
        lty = 1, lwd = lwd,
        col = cols[i+1])
X = rnorm(10000)
lines(ineq::Lc(abs(X)),
      lty = 2, lwd = lwd,
      col = 'gray70')
box()

legend(
  'topleft', bty = 'n',
  cex=1,
  # title = 'ncp',
  legend = leg,
  lty = 1, lwd = lwd,
  col = cols[2:7]
)
dev.off()

# Bias + Asymmetry ####
sam = list()
i=0
dmax = 0
for (g in c(0,0.25,0.5,0.75,1)) { # Asymmetry
  for (h in c(0,0.2,0.4,0.6)[1]) { # Tails
    i=i+1
    X = 1.2 + gh(1000000,g,h)
    d = density(X,n=2048)
    gi = gini(X)
    p = pietra(X)
    la = lasym(X)
    sam[[i]] = list(g=g,h=h,X=X,d=d,gi=gi,p=p,la=la)
    dmax = max(dmax,d$y)
  }
}
for(i in 1:length(sam))
  cat(' g = ',sam[[i]]$g,
      ' gini = ',sam[[i]]$gi,
      ' lasym = ',sam[[i]]$la,'\n')

png(file='../article/fig_bias+asym_lc.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)


leg = c()
for(i in 1:length(sam)) {
  d = sam[[i]]$d
  if(i == 1) {
    plot(
      d,
      col = cols[2],
      lwd = lwd,
      xlab = 'Errors (a.u.)',
      xlim = c(-5,5),
      xaxs = 'i',
      ylab = 'PDF',
      yaxs = 'i',
      ylim = c(0,1.05*dmax),
      main=''
    )
    grid()
  } else {
    lines(
      d,
      col = cols[i+1],
      lwd = lwd
    )
  }
  leg[i] = paste0('g = ',sam[[i]]$g,', h = ',sam[[i]]$h)
}
box()


plot(ineq::Lc(abs(sam[[1]]$X)),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
for (i in 2:length(sam))
  lines(ineq::Lc(abs(sam[[i]]$X)),
        lty = 1, lwd = lwd,
        col = cols[i+1])
X = rnorm(10000)
lines(ineq::Lc(abs(X)),
      lty = 2, lwd = lwd,
      col = 'gray70')
box()

legend(
  'topleft', bty = 'n',
  cex=1,
  # title = 'ncp',
  legend = leg,
  lty = 1, lwd = lwd,
  col = cols[2:7]
)
dev.off()
