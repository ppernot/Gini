cols     = rev(inlmisc::GetColors(8))[1:7]
cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.1))[1:7]
cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.4))[1:7]
pty      = 's'
mar      = c(3, 3, 1, 1)
mgp      = c(2, .75, 0)
tcl      = -0.5
lwd      = 4
cex      = 2.5

png(file='../article/principle.png', width=1600, height=800)
par(
  mfrow = c(1, 2),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)

# CDF ####
c = curve(
  2 * pnorm(x) - 1,
  from = 0,
  to = 3,
  xaxs = 'i', xlab = '|Error|',
  yaxs = 'i', ylab = 'Probability',
  pty = 's',
  main = 'Cumulative Distribution'
)
# grid()
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

# Lorenz ####
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
dev.off()

png(file='../article/principle_l95.png', width=800, height=800)
par(
  mfrow = c(1, 1),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex
)


# Lorenz ####
lc = ineq::Lc(abs(rnorm(1e5, 0, 1)))
clc = lc
clc$p = 1-lc$p
clc$L = 1-lc$L

plot(clc,
     lty = 2,
     lwd = lwd,
     cex = cex,
     col = 'gray70',
     xlim =c(0,0.2), xlab = 'p',
     ylim = c(0,0.8),ylab = 'l(p) = 1 - L(1-p)',
     main = ''
     )
grid()
p=0.05
abline(v=p,col = cols[6], lty = 2)
# mtext(
#   p,
#   side = 1,
#   at = p,
#   las = 0,
#   col = cols[2],
#   padj = 0.3,
#   cex = cex
# )
lpp = clc$L[which(clc$p<=p)[1]]
abline(h=lpp,col = 'gray50', lty = 2)
mtext(
  signif(lpp,2),
  side = 2,
  at = lpp,
  las = 2,
  col = 'gray50',
  adj = 1.5,
  cex = cex
)
df = c(1,2,3,5,10)
icol=1
for(nu in df) {
  icol = icol + 1
  lc = ineq::Lc(abs(rt(10000, df=nu )))
  clc = lc
  clc$p = 1-lc$p
  clc$L = 1-lc$L

  lines(clc,
        lty = 1, lwd = lwd,
        col = cols[icol]
  )
  lpp = clc$L[which(clc$p<=p)[1]]
  abline(h=lpp,col = cols[icol], lty = 2)
}
box()
legend(
  'topright', bty = 'o', bg = 'white', box.col = 'white',
  title = expression(nu),
  legend = c(df,expression(infinity)),
  lty = c(rep(1,5),2),
  lwd = lwd,
  col = c(cols[2:6],'gray50')
)
dev.off()
