#--- Set up and import libraries ----
#packrat::init()
#packrat::status()
#)

#library(pdp) #https://bgreenwell.github.io/pdp/articles/pdp.html for partial dependence — but very slow
# install.packages("data.table", repos="https://Rdatatable.github.io/data.table")

#devtools::install_github("ropensci/drake")
suppressPackageStartupMessages(library(drake))
suppressPackageStartupMessages(library(reprex))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(profvis))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(rprojroot))
suppressPackageStartupMessages(library(bigstatsr))
suppressPackageStartupMessages(library(paradox))
suppressPackageStartupMessages(library(ranger))
suppressPackageStartupMessages(library(edarf)) 
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(mlr3verse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(future.callr))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(quantreg))


#remotes::install_github("https://github.com/mlr-org/mlr3spatiotempcv")
suppressPackageStartupMessages(library(mlr3spatiotempcv))

#remotes::install_github("mlr3learners/mlr3learners.partykit", INSTALL_opts = c("--no-multiarch"), force=T)
suppressPackageStartupMessages(library(mlr3learners.partykit))

#devtools::install_github("zmjones/edarf", subdir = "pkg")
suppressPackageStartupMessages(library(edarf))

# suppressPackageStartupMessages(library(bigstatsr))
# suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(doParallel))