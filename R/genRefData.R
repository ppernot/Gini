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

# Generate reference data sets
N = 1000

# Shifted normal ####
mu = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
tab = matrix(NA,ncol=length(mu),nrow=N)
colnames(tab)=paste0('mu=',mu)
for(i in 1:length(mu)) {
  tab[,i] = -rnorm(N,mean=mu[i]) # Change sign to get correct errors sign
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_NormBias_Data.csv')
)

# Student's-t ####
df = c(1,2,3,4,5,6,7,8,9,10,20,50,100)
tab = matrix(NA,ncol=length(df),nrow=N)
colnames(tab)=paste0('df=',df)
for(i in 1:length(df)) {
  tab[,i] = -rt(N,df=df[i])
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_Student_Data.csv')
)

# Biased Student's-t ####
bias = 1
df = c(1,2,3,4,5,6,7,8,9,10,20,50,100)
tab = matrix(NA,ncol=length(df),nrow=N)
colnames(tab)=paste0('df=',df)
for(i in 1:length(df)) {
  tab[,i] = -(bias + rt(N,df=df[i]))
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_StudentBias1_Data.csv')
)

# Biased Student's-t ####
bias = 2
df = c(1,2,3,4,5,6,7,8,9,10,20,50,100)
tab = matrix(NA,ncol=length(df),nrow=N)
colnames(tab)=paste0('df=',df)
for(i in 1:length(df)) {
  tab[,i] = -(bias + rt(N,df=df[i]))
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_StudentBias2_Data.csv')
)

# g-and-h ####
bias = 0
g = seq(0.1,1,by=0.1)
tab = matrix(NA,ncol=length(g),nrow=N)
colnames(tab)=paste0('g=',g)
for(i in 1:length(g)) {
  tab[,i] = -(bias + gh(N,g[i],0))
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_GandH_Data.csv')
)

# Biased g-and-h ####
bias = 1
g = seq(0.1,1,by=0.1)
tab = matrix(NA,ncol=length(g),nrow=N)
colnames(tab)=paste0('g=',g)
for(i in 1:length(g)) {
  tab[,i] = -(bias + gh(N,g[i],0))
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_GandHBias1_Data.csv')
)

# Biased g-and-h ####
bias = 2
g = seq(0.1,1,by=0.1)
tab = matrix(NA,ncol=length(g),nrow=N)
colnames(tab)=paste0('g=',g)
for(i in 1:length(g)) {
  tab[,i] = -(bias + gh(N,g[i],0))
}
tab = cbind(Ref=rep(0,N),tab)
write.csv(
  tab,
  file = file.path('..', 'data', 'Ref_GandHBias2_Data.csv')
)
