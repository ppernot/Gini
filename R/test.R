df = c(1,2,3,5,10)
for(nu in df) {
  X = 1+ rt(1e5, df=nu )
  X = X -median(X)
  print(kurtcs(X))
  X = abs(X)
  print(q95hd(X)/median(X)-2.905)
}
X = 1+rnorm(1e5)
X = X -median(X)
print(kurtcs(X))
X = abs(X)
print(q95hd(X)/median(X)-2.905)

plot(ineq::Lc(runif(1e4,-1,1)),col=2)
