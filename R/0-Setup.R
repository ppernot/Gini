# Libraries ----

libs <- c(
  "devtools",
  "Hmisc",
  "rlist",
  "boot",
  "lmtest",
  "inlmisc",
  "distillery",
  "ineq",
  "BiocManager"
)
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(lib,
                     # repos = "https://cran.univ-paris1.fr",
                     dependencies = TRUE )
                     
  }
  library(lib, character.only = TRUE, quietly = TRUE)
}

lib = "ErrViewLib"
# if(!require(lib,character.only = TRUE))
#   devtools::install_github(paste0("ppernot/",lib))
library(lib, character.only = TRUE)
lib = "genefilter"
# if(!require(lib,character.only = TRUE))
#   BiocManager::install("genefilter", version="3.10")
library(lib, character.only = TRUE)

# Graphical parameters ----

cols = rev(inlmisc::GetColors(8))[1:7]
cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.1))[1:7]
cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.4))[1:7]
pty  = 's'
mar  = c(3, 3, 1.5, 1)
mgp  = c(2, .75, 0)
tcl  = -0.5
lwd  = 4
cex  = 4.5
reso = 1200

# Global variables ----

resultsTab = file.path('..', 'results', 'tables', 'allStats.csv')

dataSets = c('BOR2019',
             'NAR2019',
             'PER2018',
             'SCH2018',
             'THA2015', # Needs relative errors
             'WU2015',  # Needs relative errors
             'ZAS2019',
             'ZHA2018')

relSets = c('THA2015', 'WU2015') # Use relative errors

units = c('eV',
          'Debye',
          'kcal/mol',
          'kcal/mol',
          'eV',
          '%',
          '%',
          'kcal/mol',
          'eV/atom')
names(units) = dataSets

# Functions ----

skewgm_mcf = function(X, index = 1:length(X), ...) {
  X = X[index]
  X = abs(X - hrmode(X)) # Absolute mode-centered errors
  return(ErrViewLib::skewgm(X))
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

bmax = function(X, index = 1:length(X), ...) {
  X = X[index]
  m = hrmode(X)
  gfun = function (b, X)
    ErrViewLib::gini(abs(X - b))
  o = optim(
    par   = m,
    fn    = gfun,
    # method = "Brent",
    lower = min(X),
    upper = max(X),
    control = list(fnscale = -1),
    X     = X
  )
  return(o$par)
}

gmax = function(X, index = 1:length(X), ...) {
  X = X[index]
  m = hrmode(X)
  gfun = function (b, X)
    ErrViewLib::gini(abs(X - b))
  o = optim(
    par   = m,
    fn    = gfun,
    # method = "Brent",
    lower = min(X),
    upper = max(X),
    control = list(fnscale = -1),
    X     = X
  )
  return(ErrViewLib::gini(abs(X - o$par)))
}