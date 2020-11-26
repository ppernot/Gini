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

N = 100000

# Bias ####
png(file='../article/fig_NormBias.png', width=1600, height=800)
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

plot(ineq::Lc(abs(rnorm(N, mu[1], 1))),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
for (i in 2:length(mu))
  lines(ineq::Lc(abs(rnorm(N, mu[i], 1))),
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

# Student ####
png(file='../article/fig_Student.png', width=1600, height=800)
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

plot(ineq::Lc(abs(rt(N, df[1]))),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
X = rnorm(N)
lines(ineq::Lc(abs(X)),
      lty = 2, lwd = lwd,
      col = 'gray70')
for (i in 2:length(df))
  lines(ineq::Lc(abs(rt(N, df[i]))),
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

# StudentBias1 ####
png(file='../article/fig_StudentBias1.png', width=1600, height=800)
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
  dt(x-1, df[1]),
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
    dt(x-1, df[i]),
    from = -6,
    to = 6,
    col = cols[1+i], lwd = lwd,
    add = TRUE
  )
box()

plot(ineq::Lc(abs(1+rt(N, df[1]))),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
X = rnorm(N)
lines(ineq::Lc(abs(X)),
      lty = 2, lwd = lwd,
      col = 'gray70')
for (i in 2:length(df))
  lines(ineq::Lc(abs(1+rt(N, df[i]))),
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

# StudentBias2 ####
png(file='../article/fig_StudentBias2.png', width=1600, height=800)
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
  dt(x-2, df[1]),
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
    dt(x-2, df[i]),
    from = -6,
    to = 6,
    col = cols[1+i], lwd = lwd,
    add = TRUE
  )
box()

plot(ineq::Lc(abs(2+rt(N, df[1]))),
     lty =1, lwd = lwd,
     col = cols[2])
grid(col = 'gray70')
X = rnorm(N)
lines(ineq::Lc(abs(X)),
      lty = 2, lwd = lwd,
      col = 'gray70')
for (i in 2:length(df))
  lines(ineq::Lc(abs(2+rt(N, df[i]))),
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



# GAndH ####
sam = list()
i=0
dmax = 0
for (g in c(0,0.25,0.5,0.75,1)) { # Asymmetry
  for (h in c(0,0.2,0.4,0.6)[1]) { # Tails
    i=i+1
    X = 0 + gh(N,g,h)
    d = density(X,n=2048)
    gi = ErrViewLib::gini(X)
    p = ErrViewLib::pietra(X)
    la = ErrViewLib::lasym(X)
    sam[[i]] = list(g=g,h=h,X=X,d=d,gi=gi,p=p,la=la)
    dmax = max(dmax,d$y)
  }
}
# for(i in 1:length(sam))
#   cat(' g = ',sam[[i]]$g,
#       ' gini = ',sam[[i]]$gi,
#       ' lasym = ',sam[[i]]$la,'\n')

png(file='../article/fig_GAndH.png', width=1600, height=800)
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
     abline(v=0,lty=2)
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
X = rnorm(N)
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

# GAndHBias1 ####
sam = list()
i=0
dmax = 0
for (g in c(0,0.25,0.5,0.75,1)) { # Asymmetry
  for (h in c(0,0.2,0.4,0.6)[1]) { # Tails
    i=i+1
    X = 1 + gh(N,g,h)
    d = density(X,n=2048)
    gi = ErrViewLib::gini(X)
    p = ErrViewLib::pietra(X)
    la = ErrViewLib::lasym(X)
    sam[[i]] = list(g=g,h=h,X=X,d=d,gi=gi,p=p,la=la)
    dmax = max(dmax,d$y)
  }
}
# for(i in 1:length(sam))
#   cat(' g = ',sam[[i]]$g,
#       ' gini = ',sam[[i]]$gi,
#       ' lasym = ',sam[[i]]$la,'\n')

png(file='../article/fig_GANDHBias1.png', width=1600, height=800)
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
    abline(v=0,lty=2)
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

# GAndHBias1 ####
sam = list()
i=0
dmax = 0
for (g in c(0,0.25,0.5,0.75,1)) { # Asymmetry
  for (h in c(0,0.2,0.4,0.6)[1]) { # Tails
    i  = i+1
    X  = 2 + gh(N,g,h)
    d  = density(X,n=2048)
    gi = ErrViewLib::gini(X)
    p  = ErrViewLib::pietra(X)
    la = ErrViewLib::lasym(X)
    sam[[i]] = list(g=g,h=h,X=X,d=d,gi=gi,p=p,la=la)
    dmax = max(dmax,d$y)
  }
}
# for(i in 1:length(sam))
#   cat(' g = ',sam[[i]]$g,
#       ' gini = ',sam[[i]]$gi,
#       ' lasym = ',sam[[i]]$la,'\n')

png(file='../article/fig_GANDHBias2.png', width=1600, height=800)
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
    abline(v=0,lty=2)
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

# Literature datasets ####

dataSets = c(
  'BOR2019',
  'HAI2018',
  'NAR2019',
  'PER2018',
  'SCH2018',
  'THA2015', # Need relative errors
  'WU2015',  # Need relative errors
  'ZAS2019',
  'ZHA2018'
)
relSets = c('THA2015','WU2015') # Use relative errors
units = c('eV','Debye','kcal/mol','kcal/mol',
          'eV','%','%','kcal/mol','eV/atom')
names(units) = dataSets

tcols = rep(cols,6)

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
  
  png(file = file.path('..', 'article', paste0('fig_SI_',set,'.png')), 
      width=1600, height=800)
  par(
    mfrow = c(1, 2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex
  )

  plot(ecdf(abs(Errors[,1])),
       do.points=FALSE,
       lty =1, lwd = lwd,
       col = tcols[1],
       xlim = c(0,hd(abs(unlist(Errors)),0.99)),
       xlab = paste0('|Errors| (',units[set],')'),
       ylim = c(0,1), yaxs ='i', ylab ='Probability',
       main = 'ECDF')
  grid(col = 'gray70')
  # X = rnorm(10000)
  # lines(ecdf(abs(X)),
  #       do.points=FALSE,
  #       lty = 2, lwd = lwd,
  #       col = 'gray70')
  for (i in 2:ncol(Errors))
    lines(ecdf(abs(Errors[,i])),
          do.points=FALSE,
          lty = 1, lwd = lwd,
          col = tcols[i])
  box()
  
  plot(ineq::Lc(abs(Errors[,1])),
       lty =1, lwd = lwd,
       col = tcols[1])
  grid(col = 'gray70')
  X = rnorm(10000)
  lines(ineq::Lc(abs(X)),
        lty = 2, lwd = lwd,
        col = 'gray70')
  for (i in 2:ncol(Errors))
    lines(ineq::Lc(abs(Errors[,i])),
          lty = 1, lwd = lwd,
          col = tcols[i])
  box()
  
  dev.off()
}