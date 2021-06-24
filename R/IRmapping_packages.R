#--- Set up and import libraries ----
#renv::status()
#renv::hydrate()
#renv::snapshot()

#To write Rtools path to .Renviron to build packages
# if (!(file.exists("~/.Renviron"))) {
#   writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# } else {
#   renvf <- setDT(read.table("~/.Renviron"))
#   if (!(any(renvf[, grepl('[{]RTOOLS40_HOME[}]', V1)]))) {
#     print(".Renviron file exists in ")
#   }
# }

#library(pdp) #https://bgreenwell.github.io/pdp/articles/pdp.html for partial dependence — but very slow
#install.packages("data.table", repos="https://Rdatatable.github.io/data.table")

#remotes::install_github("https://github.com/mlr-org/mlr3spatiotempcv")
#Conflicts with mlr3pipelines https://github.com/mlr-org/mlr3pipelines/pull/371
suppressPackageStartupMessages(library(mlr3spatiotempcv))

#remotes::install_github("mlr3learners/mlr3learners.partykit", INSTALL_opts = c("--no-multiarch"), force=T)
suppressPackageStartupMessages(library(mlr3learners.partykit))

#devtools::install_github("zmjones/edarf", subdir = "pkg")
suppressPackageStartupMessages(library(edarf))

#remotes::install_github("zeehio/facetscales")
suppressPackageStartupMessages(library(facetscales))

suppressPackageStartupMessages(library(bigstatsr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dataRetrieval))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(drake))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(future.callr))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(gdalUtils))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(gstat))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(maps))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(mlr3))
suppressPackageStartupMessages(library(mlr3learners))
suppressPackageStartupMessages(library(mlr3pipelines))
suppressPackageStartupMessages(library(mlr3tuning))
suppressPackageStartupMessages(library(mlr3viz))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(olsrr))
suppressPackageStartupMessages(library(paradox))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(poweRlaw))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(ranger))
suppressPackageStartupMessages(library(rbin))
suppressPackageStartupMessages(library(reprex))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(rgeos))
suppressPackageStartupMessages(library("rnaturalearth"))
suppressPackageStartupMessages(library("rnaturalearthdata"))
suppressPackageStartupMessages(library(rprojroot))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(spdep))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(profvis))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(raster))
#suppressPackageStartupMessages(library(tidyhydat)) #Only use for manual inspection of Canadian gauges
suppressPackageStartupMessages(library(visNetwork))



#suppressPackageStartupMessages(library(viridis))
# suppressPackageStartupMessages(library(bigstatsr))
# suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(doParallel))
