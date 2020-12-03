library(rgl)

stat1 = 'gini'
stat2 = 'kurtcs'
stat3 = 'skewgm'
sel = D$ctd == -2 & !D$out &
  substr(D$Dataset, 1, 4) == 'Ref_' &
  D$Methods != 'df=1'
x = D[[stat1]][sel]
y = D[[stat2]][sel]
z = abs(D[[stat3]][sel])

dataSets = substring(D$Dataset[sel], first = 5)
setsNames = unique(dataSets)
setsCols = c(1, 2, 2, 2, 5, 5, 5)
names(setsCols) = setsNames
setsPchs = c(15, 15, 16, 17, 15, 16, 17)
names(setsPchs) = setsNames
icol = setsCols[dataSets]
pch  = setsPchs[dataSets]

p = plot3d(
  x, y, z,
  pch = pch,
  col = cols[icol],
  xlab = 'G',
  ylab = expression(kappa[CS]),
  zlab = expression(abs(beta[GM])),
  type = 's',
  size = 2
)
grid3d(side = 'x+')
grid3d(side = 'y+')
grid3d(side = 'z')

rgl.snapshot('../article/plot3D_Ref.png')

sel = D$ctd == -2 & !D$out &
  substr(D$Dataset, 1, 4) != 'Ref_' 
x = D[[stat1]][sel]
y = D[[stat2]][sel]
z = abs(D[[stat3]][sel])

icol = as.numeric(factor(D$Dataset[sel])) %% length(cols)
icol[icol==0] = length(cols)
pch = as.numeric(factor(D$Dataset[sel]))
sel1 = pch<=length(cols)
pch[sel1] = 16
pch[!sel1]= 17

p = plot3d(
  x, y, z,
  pch = pch,
  col = cols[icol],
  xlab = 'G',
  ylab = expression(kappa[CS]),
  zlab = expression(abs(beta[GM])),
  type = 's',
  size = 2
)
grid3d(side = 'x+')
grid3d(side = 'y+')
grid3d(side = 'z')

rgl.snapshot('../article/plot3D_Lit.png')
