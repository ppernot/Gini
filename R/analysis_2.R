# G scale thresholds ####
eps1 = 0.35
eps2 = 0.50

# Graphics parameters ####
showUnc = TRUE
xlim = c(0.1,0.7) # set for all four graphs
cols = rev(inlmisc::GetColors(8))[1:7]
pty  = 's'
mar  = c(3, 3, 1.1, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 4.5

# Get data ####
resultsTab = file.path('..', 'results', 'tables', 'allStats.csv')
D = read.csv(resultsTab)
D[['z']] = abs(D$mse / D$rmsd) # Relative bias
D[['u_z']] = D[['z']] * sqrt(D$u_mse^2/D$mse^2 + D$u_rmsd^2/D$rmsd^2)
# D[['logit_W']]  = boot::logit(D$W)


# Gini vs. Kurt ####
png(file='../article/fig_Gini_vs_CV.png', width=2400, height=2400)
#
par(mfrow = c(2, 2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

# (a) ####
stat1 = 'gini'
stat2 = 'z'
ylim = c(-0.1,5.1)
sel = D$ctd == -1 & 
  !D$out &
  substr(D$Dataset,1,4) == 'Ref_' &
  D$Methods != 'df=1' # Ignore, too unstable (Cauchy)

x = D[[stat1]][sel]
y = D[[stat2]][sel]

# Colors & symbols scheme
dataSets = substring(D$Dataset[sel],first=5)
setsNames = unique(dataSets)
setsCols = c(1,2,2,2,5,5,5)
names(setsCols) = setsNames
setsPchs = c(15,15,16,17,15,16,17)
names(setsPchs) = setsNames
icol = setsCols[dataSets]
pch  = setsPchs[dataSets]

plot(
  x,
  y,
  pch = pch,
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = '|MSE| / RMSD',
  ylim = ylim, yaxs = 'i',
  col = cols[icol],
  main = ''
)
grid()
curve(
  1/(sqrt(pi)*x),
  from=0.01,
  to=1,
  lty = 2,
  col='grey50',
  add = TRUE
  )
if(showUnc) {
  ux = D[[paste0('u_',stat1)]][sel]
  uy = D[[paste0('u_',stat2)]][sel]
  segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
  segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
}

legend(
  'topright', bty = 'o', box.col = 'white',
  legend = setsNames,
  col = cols[setsCols],
  pch = setsPchs,
  cex = 0.8
)
box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

# (b) ####
stat1 = 'gini'
stat2 = 'kurtcs'
x = D[[stat1]][sel]
y = D[[stat2]][sel]

ylim = c(-0.5,3)
plot(
  x,
  y,
  pch = pch,
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = 'Excess kurtosis',
  ylim = ylim, yaxs = 'i',
  col = cols[icol],
  main = ''
)
grid()

if(showUnc) {
  ux = D[[paste0('u_',stat1)]][sel]
  uy = D[[paste0('u_',stat2)]][sel]
  segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
  segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
}

sel = D$ctd == -1 &
  !D$out &
   D$Dataset == 'Ref_GandH'
x = D[[stat1]][sel][-1]
y = D[[stat2]][sel][-1]
abline(lm(y~x),lty=2,lwd=lwd,col='grey50')

box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

# 2nd row ####
# (c) ####

stat1 = 'gini'
stat2 = 'z'
ylim = c(-0.1,3)
sel = D$ctd == -1 &
  !D$out &
  substr(D$Dataset,1,4) == 'Ref_'

x = D[[stat1]][sel]
y = D[[stat2]][sel]
icol = factor(D$Dataset[sel])
plot(
  x,
  y,
  pch = 17, type = 'n',
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = '|MSE| / RMSD',
  ylim = ylim, yaxs = 'i',
  col = 'gray50',
  main = ''
)
grid()
rect(eps1,-2,eps2,100,col = 'gray90',border = NA)
curve(
  1/(sqrt(pi)*x),
  from=0.01,
  to=1,
  lty = 2,
  col='grey50',
  add = TRUE
)

sel = D$ctd == -1 &
  !D$out &
  substr(D$Dataset,1,4) != 'Ref_' &
  D$Dataset != 'HAI2018'
x = D[[stat1]][sel]
y = D[[stat2]][sel]
icol = as.numeric(factor(D$Dataset[sel])) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(D$Dataset[sel]))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
points(
  x, y,
  pch = pch,
  col = cols[icol]
)

if(showUnc) {
  ux = D[[paste0('u_',stat1)]][sel]
  uy = D[[paste0('u_',stat2)]][sel]
  segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
  segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
}

pch = unique(as.numeric(factor(D$Dataset[sel])))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
legend(
  'topright', bty = 'o', box.col = 'white', ncol = 1,
  legend = unique(D$Dataset[sel]),
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

# (d) ####
stat1 = 'gini'
stat2 = 'kurtcs'
sel = D$ctd == -1 &
  !D$out &
  substr(D$Dataset,1,4) == 'Ref_'
x = D[[stat1]][sel]
y = D[[stat2]][sel]
ylim = c(-1,6)
plot(
  x,
  y, type = 'n',
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = 'Excess kurtosis',
  ylim = ylim, yaxs = 'i',
  main = ''
)
grid()
rect(eps1,-2,eps2,100,col = 'gray90',border = NA)

sel = D$ctd == -1 &
  !D$out &
  D$Dataset == 'Ref_GandH'
x = D[[stat1]][sel][-1]
y = D[[stat2]][sel][-1]
abline(lm(y~x),lty=2,lwd=lwd,col='grey50')



sel = D$ctd == -1 &
  !D$out &
  substr(D$Dataset,1,4) != 'Ref_'&
  D$Dataset != 'HAI2018'
x = D[[stat1]][sel]
y = D[[stat2]][sel]
icol = as.numeric(factor(D$Dataset[sel])) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(D$Dataset[sel]))

sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
points(
  x, y,
  pch = pch,
  col = cols[icol]
)

if(showUnc) {
  ux = D[[paste0('u_',stat1)]][sel]
  uy = D[[paste0('u_',stat2)]][sel]
  segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
  segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
}


box()
mtext(
  text = '(d)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)


dev.off()



# Unbiased Gini vs. Kurt ####
png(file='../article/fig_Gini_vs_CV_unbiased.png', width=2400, height=1200)
#
par(mfrow = c(1, 2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

# (a) ####

sel = D$ctd == -1 & !D$out &
  substr(D$Dataset,1,4) != 'Ref_'&
  D$Dataset != 'HAI2018'
x1 = D$gini[sel]
u_x1 = D$u_gini[sel]
jcol= rep('gray80',length(sel))
s1 = pnorm(eps1,x1,u_x1) > 0.05 
s2 = pnorm(eps2,x1,u_x1) < 0.95 
jcol[s1] = cols[4]
jcol[s2] = cols[2]

sel = D$ctd == 0 & !D$out &
  substr(D$Dataset,1,4) != 'Ref_'&
  D$Dataset != 'HAI2018'
x2 = D$gini[sel]
u_x2 = D$u_gini[sel]
jcol2= rep('gray80',length(sel))
s12 = pnorm(eps1,x2,u_x2) > 0.05 
s22 = pnorm(eps2,x2,u_x2) < 0.95 
jcol2[s12] = cols[4]
jcol2[s22] = cols[2]

# Stats about shift sign
rho = 0.9
dif = x2 - x1
u_dif = sqrt(u_x1^2 + u_x2^2 - 2*rho*u_x1*u_x2)
print( mean(dif > 0) )
print( mean(dif < 0) )
print( mean( pnorm(0,dif,u_dif) < 0.05) ) # P(dif > 0)
print( mean( pnorm(0,dif,u_dif) > 0.95) ) # P(dif < 0)


icol = factor(D$Dataset[sel])

y1 = 3 + rnorm(length(x1),0,0.10)
y1[s1] = y1[s1] - 2
y1[s2] = y1[s2] + 2

y2 = 4 + rnorm(length(x2),0,0.10)
y2[s1] = y2[s1] - 2
y2[s2] = y2[s2] + 2

plot(
  x1, y1,
  pch = 16, type = 'p', cex=0.7,
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  yaxt = 'n', ylab = '',
  ylim = c(0.5,6.5), yaxs = 'i',
  col = jcol,
  main = ''
)
grid()
abline(v=c(eps1,eps2),lty=2,col=cols[1])
segments(x1,y1,x2,y2,col='gray70',lty=2)
points(
  x1,y1,pch=16,col=jcol, cex=0.7
)
points(
  x2,y2,pch=16,col=jcol2, cex=0.7
)
if(showUnc) {
  segments(x1 - 2 * u_x1, y1, x1 + 2 * u_x1, y1, col = jcol)
  segments(x2 - 2 * u_x2, y2, x2 + 2 * u_x2, y2, col = jcol2)
}
mtext(
  c('Raw','Centered','Raw','Centered','Raw','Centered'),
  side = 2,
  at = 1:6,
  las = 1,
  cex = 3,
  line= 0.2
  )

legend(
  'topleft',
  bty = 'n', cex = 0.7,
  legend = c(
    paste0('G < ',eps1), 
    paste0(eps1,' < G < ',eps2), 
    paste0('G > ',eps2)
  ),
  pch = 19,
  col = c(cols[4],'gray80',cols[2]),
  y.intersp = 0.9
)
# box()
mtext(
  text = '(a)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)

# (b) ####
stat1 = 'gini'
stat2 = 'kurtcs'
sel = D$ctd == -1 &
  !D$out &
  substr(D$Dataset,1,4) == 'Ref_'
x = D[[stat1]][sel]
y = D[[stat2]][sel]
ylim = c(-1,6)
plot(
  x,
  y, type = 'n',
  xlab = 'Gini',
  xlim = xlim, xaxs = 'i',
  ylab = 'Excess kurtosis',
  ylim = ylim, yaxs = 'i',
  main = ''
)
grid()
rect(eps1,-2,eps2,100,col = 'gray90',border = NA)

sel = D$ctd == -1 &
  !D$out &
  D$Dataset == 'Ref_GandH'
x = D[[stat1]][sel][-1]
y = D[[stat2]][sel][-1]
abline(lm(y~x),lty=2,lwd=lwd,col='grey50')

sel = D$ctd == 0 &
  !D$out &
  substr(D$Dataset,1,4) != 'Ref_'&
  D$Dataset != 'HAI2018'
x = D[[stat1]][sel]
y = D[[stat2]][sel]
icol = as.numeric(factor(D$Dataset[sel])) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(D$Dataset[sel]))

sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
points(
  x, y,
  pch = pch,
  col = cols[icol]
)

if(showUnc) {
  ux = D[[paste0('u_',stat1)]][sel]
  uy = D[[paste0('u_',stat2)]][sel]
  segments(x, y - 2 * uy, x, y + 2 * uy, col = cols[icol])
  segments(x - 2 * ux, y, x + 2 * ux, y, col = cols[icol])
}

pch = unique(as.numeric(factor(D$Dataset[sel])))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17
legend(
  'topleft', bty = 'o', box.col = 'white', ncol = 1,
  legend = paste0('c-',unique(D$Dataset[sel])),
  col = unique(cols[icol]),
  pch = pch,
  cex=0.8
)
box()
mtext(
  text = '(b)',
  side = 3,
  adj = 1,
  cex = cex,
  line = 0.3)


dev.off()

# MUE vs Gini ####

png(file='../article/fig_Gini_vs_Rank.png', width=2400, height=1200)
mar  = c(3, 5, 1.1, 1)
par(mfrow = c(1, 2),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    lwd = lwd,
    cex = cex)

sel = D$ctd == -1 & !D$out &
  substr(D$Dataset,1,4) != 'Ref_'&
  D$Dataset != 'HAI2018'
setsNames = unique(D$Dataset[sel])
icol = as.numeric(factor(D$Dataset[sel]))
d = D[sel,]
d[['rmue']] = d$mue
for(s in setsNames) {
  s1 = which(d$Dataset == s)
  d$rmue[s1] = rank(d$mue[s1])
}

jcol= rep('gray80',length(sel))
s1 = pnorm(eps1,d$gini,d$u_gini) > 0.05 #d$gini <= 0.35
s2 = pnorm(eps2,d$gini,d$u_gini) < 0.95 #d$gini >= 0.5
jcol[s1] = cols[4]
jcol[s2] = cols[2]

plot(d$rmue, icol, pch=19, col=jcol,
     xlim = c(0.8, 10.2), xlab = 'rank(MUE)',
     yaxt = 'n', ylab ='', ylim = c(0.5,8.5))
mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
legend(
  6.5, 8.8,
  bty = 'n', cex = 0.7,
  legend = c(
    paste0('G < ',eps1), 
    paste0(eps1,' < G < ',eps2), 
    paste0('G > ',eps2)
  ),
  pch = 19,
  col = c(cols[4],'gray80',cols[2]),
  y.intersp = 0.9
)
box()

# Remove global outliers ####
sel = D$ctd == -1 & D$out &
  substr(D$Dataset,1,4) != 'Ref_'&
  D$Dataset != 'HAI2018'
icol = as.numeric(factor(D$Dataset[sel]))
setsNames = unique(D$Dataset[sel])
d = D[sel,]
d[['rmue']] = d$mue
for(s in setsNames) {
  s1 = which(d$Dataset == s)
  d$rmue[s1] = rank(d$mue[s1])
}

jcol= rep('gray80',length(sel))
s1 = pnorm(eps1,d$gini,d$u_gini) > 0.05 #d$gini <= 0.35
s2 = pnorm(eps2,d$gini,d$u_gini) < 0.95 #d$gini >= 0.5
jcol[s1] = cols[4]
jcol[s2] = cols[2]

plot(d$rmue, icol, pch=19, col=jcol,
     xlim = c(0.8, 10.2), xlab = 'rank(MUE)',
     yaxt = 'n', ylab ='', ylim = c(0.5,8.5))
mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
legend(
  6, 8.5,
  bty = 'n', cex = 0.7,
  title = 'Removed outliers',
  legend = '',
  pch = -1
)
box()

# # Q95 vs Gini ####
# sel = D$ctd == -1 & !D$out &
#   substr(D$Dataset,1,4) != 'Ref_'
# setsNames = unique(D$Dataset[sel])
# d = D[sel,]
# d[['rmue']] = d$mue
# for(s in setsNames) {
#   s1 = which(d$Dataset == s)
#   d$rmue[s1] = rank(d$q95hd[s1])
# }
#
# jcol= rep('gray80',length(sel))
# s1 = pnorm(0.35,d$gini,d$u_gini) > 0.05 #d$gini <= 0.35
# s2 = pnorm(0.50,d$gini,d$u_gini) < 0.95 #d$gini >= 0.5
# jcol[s1] = cols[4]
# jcol[s2] = cols[2]
#
# plot(d$rmue, icol, pch=19, col=jcol,
#      xlim = c(0.8, 10.2), xlab = 'rank(Q95)',
#      yaxt = 'n', ylab ='', ylim = c(0.5,9.5))
# mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
# legend(
#   6.5, 9.8,
#   bty = 'n', cex = 0.7,
#   legend = c('G < 0.35', '0.35 < G < 0.5', 'G > 0.5'),
#   pch = 19,
#   col = c(cols[4],'gray80',cols[2]),
#   y.intersp = 0.9
# )
# box()
#
# # Remove global outliers ####
# sel = D$ctd == -1 & D$out &
#   substr(D$Dataset,1,4) != 'Ref_'
# icol = as.numeric(factor(D$Dataset[sel]))
# setsNames = unique(D$Dataset[sel])
# d = D[sel,]
# d[['rmue']] = d$mue
# for(s in setsNames) {
#   s1 = which(d$Dataset == s)
#   d$rmue[s1] = rank(d$q95hd[s1])
# }
#
# jcol= rep('gray80',length(sel))
# s1 = pnorm(0.35,d$gini,d$u_gini) > 0.05 #d$gini <= 0.35
# s2 = pnorm(0.50,d$gini,d$u_gini) < 0.95 #d$gini >= 0.5
# jcol[s1] = cols[4]
# jcol[s2] = cols[2]
#
# plot(d$rmue, icol, pch=19, col=jcol,
#      xlim = c(0.8, 10.2), xlab = 'rank(Q95)',
#      yaxt = 'n', ylab ='', ylim = c(0.5,9.5))
# mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
# legend(
#   6, 9.5,
#   bty = 'n', cex = 0.7,
#   title = 'Removed outliers',
#   legend = '',
#   pch = -1
# )
# box()

dev.off()

# # Q95 vs Gini ####
# sel = D$ctd == -1 &
#   !D$out &
#   substr(D$Dataset,1,4) != 'Ref_'
# icol = as.numeric(factor(D$Dataset[sel]))
# setsNames = unique(D$Dataset[sel])
# d = D[sel,]
# d[['rmue']] = d$mue
# for(s in setsNames) {
#   s1 = which(d$Dataset == s)
#   d$rmue[s1] = rank(d$q95hd[s1]) #/min(d$mue[sel])
# }
#
# png(file='../article/fig_Gini_vs_RankQ95.png', width=800, height=800)
# cex = 3
# mar  = c(3, 5, 1.1, 1)
# par(mfrow = c(1, 1),
#     pty = pty,
#     mar = mar,
#     mgp = mgp,
#     tcl = tcl,
#     lwd = lwd,
#     cex = cex)
#
# jcol= rep('gray80',length(sel))
# s1 = d$gini <= 0.35
# s2 = d$gini >= 0.5
# jcol[s1] = cols[4]
# jcol[s2] = cols[2]
#
# plot(d$rmue, icol, pch=19, col=jcol,
#      xlim = c(0.5, 6.5), xlab = 'rank(Q95)',
#      yaxt = 'n', ylab ='', ylim = c(0.5,9.5))
# mtext(setsNames,side =2, at =1:max(icol), las=1, adj=1.1, cex=cex)
# legend(
#   4, 10,
#   bty = 'n', cex = 0.7,
#   legend = c('G < 0.35', '0.35 < G < 0.5', 'G > 0.5'),
#   pch = 19,
#   col = c(cols[4],'gray80',cols[2]),
#   y.intersp = 0.9
# )
# box()
# dev.off()
