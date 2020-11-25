cols     = rev(inlmisc::GetColors(8))[1:7]
cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.1))[1:7]
cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.4))[1:7]
pty      = 's'
mar      = c(3, 3, 1, 1)
mgp      = c(2, .75, 0)
tcl      = -0.5
lwd      = 4
cex      = 2.5

png(file='../article/fig_ambiguity.png', width=800, height=800)
par(
  mfrow = c(1, 1),
  pty = pty,
  mar = mar,
  mgp = mgp,
  tcl = tcl,
  lwd = lwd,
  cex = cex,
  xaxs = "i",
  yaxs = "i"
)

mue = 1

mu1 = 1
s1  = 0.1

(p1 = 1-2*(pnorm(1,0,sqrt(pi / 2))-0.5))

curve(dnorm(x, 0, sqrt(pi / 2)),
      from = -2.75,
      to = 2.75,
      lwd = lwd,
      col = cols[2],
      xlab = "Error",
      ylab = "Probability Density Function",
      ylim = c(0,4.5))
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
