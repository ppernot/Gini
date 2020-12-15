source('./0-Setup.R')

numDig = 2 # Number of digits to report in uncertainty

# Read data ----
D = read.csv(resultsTab)

# Table 2 ----
sets = D$ctd == -1 & # No linear correction
  ! D$out       # No outliers removed
d = D[sets,]

io = order(d$gimc, decreasing = TRUE)[1:10]
df = data.frame(
  Dataset = d$Dataset[io], 
  Methods = d$Methods[io])
stats = c('gimc','skewgm','kurtcs')
for(s in stats) {
  v = d[[s]][io]
  vu = matrix(apply(cbind(v, d[[paste0('u_', s)]][io]), 1,
                    function(x)
                      ErrViewLib::prettyUnc(
                        x[1], x[2], numDig = numDig)),
              ncol = 1)
  colnames(vu) = s
  df = cbind(df, vu)
}

sink(file.path('..', 'results', 'tables', 'Table_2.tex'))
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'The ten methods with the largest GMCF values, 
               and the corresponding skewness and kurtosis.',
    label = "tab:statsLit"
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()

# Table 3 ----
sets = D$Dataset == 'ZHA2018' & 
  !D$out         # No outliers removed
d = D[sets,]

stats = c('mue','mse','q95hd','skewgm','kurtcs','gimc')

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

sink(file.path('..', 'results', 'tables', 'Table_3.tex'))
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Statistics for the methods of the ZHA2018 dataset.',
    label = "tab:statsZHA"
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()


