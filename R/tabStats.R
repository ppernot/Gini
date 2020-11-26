numDig = 2

# Generate table of statistics for datasets

resultsTab = file.path('..', 'results', 'tables', 'allStats.csv')
D = read.csv(resultsTab)
D[['1/cv']] = abs(D$mse / D$rmsd)
D[['u_1/cv']] = abs(D$mse / D$rmsd) *
  sqrt(D$u_mse^2/D$mse^2 + D$u_rmsd^2/D$rmsd^2)

# Reference datasets ####

sets = D$ctd == -1 &
  D$out &
  substr(D$Dataset,1,4) == 'Ref_'&
  D$Methods != 'df=1' # Ignore, too unstable (Cauchy)

stats = c('mue','mse','rmsd','q95hd','kurtcs','gini','1/cv')
u_stats = paste0('u_',stats)
cols = c('Dataset','Methods',stats,u_stats)
d = D[sets,cols]
d$Dataset = substring(d$Dataset,first=5)

df = data.frame(Dataset = d$Dataset, Methods = d$Methods)
for(s in stats) {
  v = d[[s]]
  vu = matrix(apply(cbind(v, d[[paste0('u_', s)]]), 1,
                    function(x)
                      prettyUnc(x[1], x[2], numDig = numDig)),
              ncol = 1)
  colnames(vu) = s
  df = cbind(df, vu)
}

sink(file.path('..', 'results', 'tables', 'refStats.tex'))
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Statistics of the reference datasets',
    label = "tab:statsRef"
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()

# Literature datasets ####

sets = D$ctd == -1 &
  D$out &
  substr(D$Dataset,1,4) != 'Ref_'

stats = c('mue','mse','rmsd','q95hd','kurtcs','gini','1/cv')
u_stats = paste0('u_',stats)
cols = c('Dataset','Methods',stats,u_stats)
d = D[sets,cols]
# d$Dataset = substring(d$Dataset,first=5)

df = data.frame(Dataset = d$Dataset, Methods = d$Methods)
for(s in stats) {
  v = d[[s]]
  vu = matrix(apply(cbind(v, d[[paste0('u_', s)]]), 1,
                    function(x)
                      prettyUnc(x[1], x[2], numDig = numDig)),
              ncol = 1)
  colnames(vu) = s
  df = cbind(df, vu)
}

sink(file.path('..', 'results', 'tables', 'litStats.tex'))
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Statistics of the literature datasets',
    label = "tab:statsLit"
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()
