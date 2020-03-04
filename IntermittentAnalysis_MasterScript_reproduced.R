#Intermittent river analysis master script 
#Author: Mathis L. Messager
#Contact info: mathis.messager@mail.mcgill.ca
#Affiliation: 
#Global HydroLAB, Department of Geography, McGill University
#EcoFlows Lab, RiverLy Research Unit, INRAE Lyon



#---- Set up and import libraries ----
#packrat::init()
packrat::status()

library(reprex)
library(tictoc)
library(profvis)
library(stringr)
# install.packages("data.table", repos="https://Rdatatable.github.io/data.table")
# packrat::snapshot()
library(plyr)
library(ggplot2)
library(gridExtra)
library(sf)
library(rprojroot)
require(bigstatsr)
require(parallel)
require(doParallel)
library(paradox)
library(ranger)
#library(pdp) #https://bgreenwell.github.io/pdp/articles/pdp.html for partial dependence — but very slow
library(edarf) 
library(xlsx)
library(ggpubr)
require(mlr3verse)
library(data.table)
#library(vroom)

#### -------------------------- Define directory structure --------------------------------
#Get mDur and mFreq values for all stations 
rootdir <- find_root(has_dir("src"))
datdir <- file.path(rootdir, 'data')
resdir <- file.path(rootdir, 'results')

outgdb <- file.path(resdir, 'spatialoutputs.gdb')

#Import GRDC stations and only keep those < 200 m from a HydroSHEDS reach
GRDCpjoin <- st_read(dsn=outgdb, layer='GRDCstations_riverjoin') %>%
  .[.$station_river_distance<200,]

#Get data paths of daily records for GRDC stations
fileNames <- file.path(datdir, 'GRDCdat_day', paste(GRDCpjoin$GRDC_NO,  ".txt", sep=""))
which(do.call(rbind, lapply(fileNames, file.exists))) #Check whether any file does not exist


#### -------------------------- Define functions ------------------------------------------
diny <- function(year) {
  #Takes a year and returns the number of days within it (366 for leap years)
  365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
}

zero.lomf <- function(x, first=TRUE) {
  #Finds the index, for each row, of the previous row with a non-zero value 
  #Takes in a univariate time series with zeros'''
  #If first value is zero, will asign that day as the previous day of flow 
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if(first==T & x[1]==0)
      non.zero.idx=c(1,non.zero.idx)
    rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1))) #Repeat index of previous row with non-NA as many times gap until next non-NA values
  }
}

durfreq_indiv <- function(path, maxgap, monthsel=NULL) {
  #For a given GRDC gauge discharge file (given by path), format and compute 
  # firstYear      : num first year on full record
  # lastYear       : num last year on full record
  # totalYears     : int total number of years on full record
  # firstYear_kept : num first year on record with < maxgap missing days
  # lastYear_kept  : num first year on record with < maxgap missing days
  # totalYears_kept: int total number of years with < maxgap missing days
  # totaldays      : num total number of days with discharge data
  # sumDur         : int total number of days with discharge = 0
  # mDur           : num mean number of days/year with discharge = 0
  # mFreq          : num mean number of periods with discharge = 0 (with at least one day of flow between periods)
  
  gaugetab <- cbind(fread(path, header=T, skip = 40, sep=";", colClasses=c('character', 'character', 'numeric', 'numeric', 'integer')),
                    GRDC_NO = strsplit(basename(path), '[.]')[[1]][1])%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)
  
  gaugetab[, year:=as.numeric(substr(dates, 1, 4))] %>% #Create year column and total number of years on record
    .[, prevflowdate := gaugetab[zero.lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA]
  
  gaugetab[!(Original == -999 | is.na(Original)), `:=`(missingdays = diny(year)-.N,
                                                       datadays = .N), by= 'year'] #Compute number of missing days per year
  
  
  gaugetab[,month:=as.numeric(substr(dates, 6, 7))] #create a month column
  
  if (gaugetab[(missingdays < maxgap), .N>0]) { #Make sure that there are years with sufficient data
    monthlyfreq <- gaugetab[(missingdays < maxgap), 
                            .(GRDC_NO = unique(GRDC_NO),
                              monthrelfreq = length(unique(na.omit(prevflowdate)))/length(unique(year))), by='month'] %>% #Count proportion of years with zero flow occurrence per month
      dcast(GRDC_NO~month, value.var='monthrelfreq') %>%
      setnames(as.character(seq(1,12)), c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  } else {
    monthlyfreq <- data.table(GRDC_NO = gaugetab[, unique(GRDC_NO)], month=1:12, monthrelfreq=rep(NA,12)) %>%
      dcast(GRDC_NO~month, value.var='monthrelfreq') %>%
      setnames(as.character(seq(1,12)), c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  }
  
  #If analysis is only performed on a subset of months
  if (!is.null(monthsel)) {
    gaugetab <- gaugetab[month %in% monthsel, ]
  }
  
  gaugetab_yearly <- merge(gaugetab[, .(missingdays=max(missingdays, na.rm=T), 
                                        datadays=max(datadays, na.rm=T)), by='year'],
                           gaugetab[Original == 0, .(dur=.N,
                                                     freq=length(unique(prevflowdate))), by='year'], #Compute number of days of zero flow for years with number of gap days under threshold
                           by = 'year', all.x = T) %>%
    .[!is.finite(missingdays), missingdays := diny(year)] %>%
    .[!is.finite(datadays), datadays := diny(year)] %>%
    .[dur==diny(year), freq:=1] %>%
    .[is.na(dur), `:=`(dur=0, freq=0)]
  
  
  gaugetab_all <- cbind(monthlyfreq,
                        gaugetab_yearly[, .(firstYear=min(year), 
                                            lastYear=max(year), 
                                            totalYears=length(unique(year)))],
                        gaugetab_yearly[missingdays <= maxgap, .(firstYear_kept=min(year), 
                                                                 lastYear_kept=max(year), 
                                                                 totalYears_kept=length(unique(year)),
                                                                 totaldays = sum(datadays), 
                                                                 sumDur = sum(dur),
                                                                 mDur = mean(dur),
                                                                 mFreq = mean(freq),
                                                                 intermittent = factor(ifelse(sum(dur)>0, 1, 0), levels=c('0','1')))]
  )
  return(gaugetab_all)
}

durfreq_parallel <- function(pathlist, maxgap, monthsel_list=NULL, reverse=FALSE) {
  #Run durfreq_indiv in parallel
  #Can compute mean number of zero flow days, frequency of zero flow periods, and intermittency on selected months
  #monthsel_list: named list with names being GRDC_NO and values the selected months
  #reverse: whether to use monthsel_list as excluding months rather than subsetting
  
  cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)
  
  return(as.data.table(rbindlist(foreach(j=pathlist, 
                                         .packages = c("data.table", "magrittr"), 
                                         .export=c('durfreq_indiv', 'zero.lomf', 'diny')) 
                                 %dopar% {
                                   if (is.null(monthsel_list)) {
                                     return(durfreq_indiv(j, maxgap = 20))
                                   } else {
                                     monthsublist = monthsel_list[[strsplit(basename(j), '[.]')[[1]][1]]]
                                     
                                     #Run duration and frequency computation for subset of months
                                     return(durfreq_indiv(j, maxgap = maxgap, 
                                                          monthsel= if (reverse) setdiff(1:12, monthsublist) else monthsublist)
                                     )
                                   }
                                 }
  )))
}

pdtiles_grid <- function(mod, dt, colnums) {
  "Create a plot matrix of partial dependence interactions"
  
  allvar <- mod$variable.importance[order(-mod$variable.importance)]
  selcols <- names(allvar)[1:colnums]
  pdout <- edarf::partial_dependence(mod, vars= selcols, 
                                     n=c(10,10), interaction=TRUE, data=as.data.frame(dt)) %>% #Warning: doe snot work with data_table
    setDT 
  
  #Iterate over every pair of PCs
  tileplots_l <- sapply(seq_len(length(selcols)-1), function(i) {
    sapply(seq(2, length(selcols)), function(j) {
      print(paste0('Generating plot', i , j))
      xvar <- selcols[i]
      yvar <- selcols[j]
      
      if (i != j) {
        ggplotGrob(
          ggplot(pdout,  aes_string(x=xvar, y=yvar)) +
            geom_tile(aes(fill = `1`)) + 
            scale_fill_distiller(palette='Spectral') + 
            geom_jitter(data=dt, aes(color=intermittent)) +
            theme_bw()
        )
      } else {
        #Write proportion of variance explained by PC
        text_grob(round(allvar[1:colnums]), 3)
      }
    })
  })
  #Plot, only keeping unique combinations by grabbing lower triangle of plot matrix
  do.call("grid.arrange", list(grobs=tileplots_l[lower.tri(tileplots_l, diag=T)], ncol=colnums-1)) 
}

fread_cols <- function(file_name, colsToKeep) {
  header <- fread(file_name, nrows = 1, header = FALSE)
  keptcols <- colsToKeep[colsToKeep %chin% unlist(header)]
  missingcols <- colsToKeep[!(colsToKeep %chin% unlist(header))]
  paste('Importing', file_name, 'with ', length(keptcols), 'columns out of ', length(colsToKeep), 'supplied column names')
  fread(input=file_name, header=TRUE, select=keptcols, verbose=TRUE)
}

#### -------------------------- Format and inspect data -------------------------------------------------
#---- Get mean drying duration and frequency per year ----
tic()
gaugestats <- durfreq_parallel(pathlist=fileNames, maxgap=20) 
toc()

#Join intermittency statistics to predictor variables and subset to only include those gauges with at least
gaugestats_join <- gaugestats[as.data.table(GRDCpjoin), on='GRDC_NO'] %>%
  .[!is.na(totalYears_kept) & totalYears_kept>10,] # Only keep stations with at least 10 years of data

#---- Compute derived predictor variables ----
gaugestats_join[, `:=`(pre_mm_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), 
                       pre_mm_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
                .SDcols= paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Compute minimum and maximum catchment precipitation
  .[, `:=`(pre_mm_cvar=pre_mm_cmn/pre_mm_cmx, #min/max monthly catchment precip
           dis_mm_pvar=ifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx), #min/max monthly watershed discharge
           dis_mm_pvaryr=ifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),#min monthly/average yearly watershed discharge
           ele_pc_rel = ifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav))] #catchment average elv - watershec average elev

gaugestats_join[, cmi_ix_cmn := do.call(pmin, c(.SD)),
                .SDcols= paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))] %>%
  .[, swc_pc_cmn := do.call(pmin, c(.SD)),
    .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))]

#Add new variables to meta_format (variable labels)
addedvars <- data.table(varname=c('Precipitation catchment Annual minimum', 'Precipitation catchment Annual maximum', 
                                  'Precipitation catchment Annual min/max', 
                                  'Climate Moisture Index catchment Annual minimum',
                                  'Discharge watershed Annual min/max', 'Discharge watershed Annual min/average',
                                  'Elevation catchment average - watershed average',
                                  'Soil Water Content catchment Annual minimum'), 
                        varcode=c('pre_mm_cmn', 'pre_mm_cmx', 'pre_mm_cvar', 
                                  'cmi_ix_cmn',
                                  'dis_mm_pvar', 'dis_mm_pvaryr',
                                  'ele_pc_rel',
                                  'swc_pc_cmn'))

#---- Convert -9999 values to NA ----
for (j in which(sapply(gaugestats_join,is.numeric))) { #Iterate through numeric column indices
  set(gaugestats_join,which(gaugestats_join[[j]]==-9999),j, NA)} #Set those to 0 if -9999


#---- Select predictor variables ----
rfpredcols<- c('dis_m3_pyr',
               'dis_m3_pmn',
               'dis_m3_pmx',
               'dis_mm_pvar',
               'dis_mm_pvaryr',
               'run_mm_cyr',
               'inu_pc_umn',
               'inu_pc_umx',
               'inu_pc_cmn',
               'lka_pc_cse',
               'lka_pc_use',
               'dor_pc_pva',
               'gwt_cm_cav',
               'ele_pc_rel',
               # 'sgr_dk_rav', #Don't use stream gradient as data are missing for all of Greenland due to shitty DEM
               'clz_cl_cmj',
               'tmp_dc_cyr',
               'tmp_dc_cmn',
               'tmp_dc_cmx',
               'tmp_dc_uyr',
               'pre_mm_uyr',
               'pre_mm_cvar',
               'pre_mm_cmn',
               'pet_mm_uyr',
               'ari_ix_uav',
               'ari_ix_cav',
               'cmi_ix_uyr',
               'cmi_ix_cmn',
               'snw_pc_uyr',
               'snw_pc_cyr',
               'snw_pc_cmx',
               'glc_cl_cmj',
               'pnv_cl_cmj',
               'wet_pc_cg1',
               'wet_pc_cg2',
               'wet_pc_ug1',
               'wet_pc_ug2',
               'for_pc_use',
               'for_pc_cse',
               'ire_pc_use',
               'ire_pc_cse',
               'gla_pc_use',
               'gla_pc_cse',
               'prm_pc_use',
               'prm_pc_cse',
               'cly_pc_uav',
               'cly_pc_cav',
               'slt_pc_uav',
               'slt_pc_cav',
               'snd_pc_uav',
               'snd_pc_cav',
               'soc_th_uav',
               'soc_th_cav',
               'swc_pc_uyr',
               'swc_pc_cyr',
               'swc_pc_cmn',
               'lit_cl_cmj',
               'kar_pc_use',
               'kar_pc_cse')

#Check that all columns are in dt
rfpredcols[!(rfpredcols %in% names(gaugestats_join))]

#---- Associate column names with variables names ----

#Get predictor variable names
riveratlas_metapath <- file.path(datdir, 'HydroATLAS', 'HydroATLAS_metadata_MLM.xlsx')
metaall <- read.xlsx(riveratlas_metapath, sheetName='Overall') %>% setDT 
metascale <- read.xlsx(riveratlas_metapath, sheetName='scale') %>% setDT %>% setnames('Key', 'Keyscale')
metastat <- read.xlsx(riveratlas_metapath, sheetName='stat') %>% setDT %>% setnames('Key', 'Keystat')

meta_format <- as.data.table(expand.grid(Column.s.=metaall$Column.s., Keyscale=metascale$Keyscale, Keystat=metastat$Keystat))[
  metaall, on='Column.s.'][
    metascale, on = 'Keyscale'][
      metastat, on = 'Keystat', allow.cartesian=TRUE]

meta_format[, `:=`(varcode = paste0(gsub('[-]{3}', '', Column.s.), Keyscale, Keystat),
                   varname = paste(Attribute, Spatial.representation, Temporal.or.statistical.aggregation.or.other.association))]

rfpredcols_dt <- merge(data.table(varcode=rfpredcols), rbind(meta_format, addedvars, fill=T), by='varcode', all.y=F)

#### -------------------------- Run RF model -------------------------------------------------
#---- Tune model ----
# rf_formula <- as.formula(paste0('intermittent~', 
#                                 paste(rfpredcols, collapse="+"), 
#                                 collapse=""))

#Create taxk
task_inter <- mlr3::TaskClassif$new(id ='inter_basic',
                                    backend = gaugestats_join[!is.na(cly_pc_cav), c('intermittent', rfpredcols),with=F],
                                    target = "intermittent")

#Create learner
lrn_ranger <- mlr3::lrn('classif.ranger', num.trees=500, replace=F, predict_type="prob")
lrn_ranger$param_set

#Define parameter ranges to tune on
tune_ranger = ParamSet$new(list(
  ParamInt$new("mtry", lower = 1, upper = 11),
  ParamInt$new("min.node.size", lower = 1, upper = 10),
  ParamDbl$new("sample.fraction", lower = 0.1, upper = 0.9)
))
tune_ranger

#Define resampling strategy
rcv_ranger = rsmp("repeated_cv", repeats=10, folds=5) #5-fold aspatial CV repeated 10 times
#Define performance measure
measure_ranger = msr("classif.bacc") #use balanced accuracy to account for imbalanced dataset
#Define termination rule 
evals20 = term("evals", n_evals = 20) #termine tuning after 20 rounds

#Create tuning insance
instance_ranger = TuningInstance$new(
  task = task_inter,
  learner = lrn_ranger,
  resampling = rcv_ranger,
  measures = measure_ranger,
  param_set = tune_ranger,
  terminator = evals20
)
print(instance_ranger)

#Define tuning algorithm
randomtuner = tnr("random_search", batch_size = 4L)

#Launch tuning
result = randomtuner$tune(instance_ranger)
print(result)

#Use output tuning parameters for final model training
lrn_ranger$param_set$values =  mlr3misc::insert_named(
  instance_ranger$result$params, 
  list(importance='permutation')
)

lrn_ranger$train(task_inter)

#Get output trained ranger model
ranger_basictrained <- lrn_ranger$model




#### -------------------------- Diagnose initial model -------------------------------------------------
# ---- Variable importance ----

varimp_basic <- data.table(varcode=names(ranger_basictrained$variable.importance), 
                           importance=ranger_basictrained$variable.importance)[
                             rfpredcols_dt, on='varcode' ]%>%
  .[, varname := factor(varname, levels=.[,varname[order(-importance)]])] %>%
  setorder(-importance)

ggplot(varimp_basic[1:30,],aes(x=varname, y=importance)) + 
  geom_bar(stat = 'identity') +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=8))

# ---- Check misclassification rate with different threshold probabilities ----
gaugestats_join[!is.na(cly_pc_cav), intermittent_predprob := ranger_basictrained$predictions[,'1']]

threshold_misclass <- ldply(seq(0,1,0.01), function(i) {
  gaugestats_join[!is.na(cly_pc_cav), intermittent_predcat := ifelse(intermittent_predprob>=i, 1, 0)]
  confumat <- gaugestats_join[!is.na(cly_pc_cav), .N, by=c('intermittent', 'intermittent_predcat')]
  outvec <- data.table(i, 
                       `Misclassification rate`=gaugestats_join[intermittent!=intermittent_predcat,.N]/gaugestats_join[,.N],
                       `True positive rate (sensitivity)` = confumat[intermittent=='1' & intermittent_predcat==1, N]/confumat[intermittent=='1', sum(N)],
                       `True negative rate (specificity)` = confumat[intermittent=='0' & intermittent_predcat==0, N]/confumat[intermittent=='0', sum(N)])
  return(outvec)
}) %>% setDT

ggplot(melt(threshold_misclass, id.vars='i'), aes(x=i, y=value, color=variable)) + 
  geom_line(size=1.2) + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme_bw()

threshold_misclass[i %in% c(0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.35, 0.5, 0.45, 0.5),]
gaugestats_join[, intermittent_predcat := ifelse(intermittent_predprob>=0.27, 1, 0)]

# ---- Partial dependence plots ----
pdtiles_grid(mod = ranger_basictrained, dt = gaugestats_join, colnums = 5)

# ---- Output GRDC predictions as points ----
st_write(obj=merge(GRDCpjoin, gaugestats_join[, !(colnames(GRDCpjoin)[colnames(GRDCpjoin) != 'GRDC_NO']) , with=F], by='GRDC_NO'),
         dsn=file.path(resdir, 'GRDCstations_predbasic800.gpkg'), driver = 'gpkg', delete_dsn=T)


# ---- Generate predictions to map on river network ----
#Read in river network attribute table
riveratlas <- fread_cols(file.path(resdir, 'RiverATLAS_v10tab.csv'), 
                         colsToKeep = c("HYRIV_ID", rfpredcols, 'ele_mt_cav', 'ele_mt_uav',
                                        paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
                                        paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
                                        paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))))

#Inspect -9999 values for cmi_ix_uyr
check <- riveratlas[cmi_ix_uyr == -9999,] #All have precipitation = 0
riveratlas[cmi_ix_uyr == -9999, cmi_ix_uyr := 0]

check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
riveratlas[snw_pc_cyr == -9999, snw_pc_cyr:=0]
riveratlas[snw_pc_cmx == -9999, snw_pc_cmx:=0]

check <- riveratlas[is.na(sgr_dk_rav),]

#Convert -9999 values to NA
colNAs<- riveratlas[, lapply(.SD, function(x) sum(is.na(x) | x==-9999))]

for (j in which(sapply(riveratlas,is.numeric))) { #Iterate through numeric column indices
  set(riveratlas,which(riveratlas[[j]]==-9999),j, NA)} #Set those to 0 if -9999

#Compute derived variables
riveratlas[, `:=`(pre_mm_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), 
                  pre_mm_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
           .SDcols= paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Compute minimum and maximum catchment precipitation
  .[, `:=`(pre_mm_cvar=ifelse(pre_mm_cmx==0, 1, pre_mm_cmn/pre_mm_cmx), #min/max monthly catchment precip
           dis_mm_pvar=ifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx), #min/max monthly watershed discharge
           dis_mm_pvaryr=ifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),#min monthly/average yearly watershed discharge
           ele_pc_rel = ifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav))] #catchment average elv - watershec average elev

riveratlas[, cmi_ix_cmn := do.call(pmin, c(.SD)),
           .SDcols= paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Compute minimum catchment monthly moisture index
  .[is.na(cmi_ix_cmn), cmi_ix_cmn := 0] %>% #if cmin_ix_cmn is na (because pre_cmn == 0, set it to 0)
  .[, swc_pc_cmn := do.call(pmin, c(.SD)), #Compute catchment annual minimum soil water content
    .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))]

#Predict model (should take ~10-15 minutes for 9M rows x 40 cols — chunk it up to avoid memory errors)
for (clz in unique(riveratlas$clz_cl_cmj)) {
  print(clz)
  tic()
  riveratlas[!is.na(cly_pc_cav) & !is.na(cly_pc_uav) & clz_cl_cmj == clz, 
             predbasic800 := predict(ranger_basictrained, type='response', data=.SD)$predictions[,'1'], 
             .SDcols = c("HYRIV_ID", rfpredcols)]
  toc()
}

riveratlas[, predbasic800cat := ifelse(predbasic800>=0.3, 1, 0)]
fwrite(riveratlas[, c('HYRIV_ID', 'predbasic800cat'), with=F], file.path(resdir, 'RiverATLAS_predbasic800.csv'))










############# TO DO ##############
## Map uncertainty in predictions for each gauge — relate to length of record and environmental characteristics
## Compare gauge environmental characteristics compared to full river network, looking at confusion matrix results
## Understand variable associations with gauges
#Identify gaps in reference streamgauge data (multidimensional space)

######## Model improvements
#Implement nested resampling
#Check spatial resampling
#Implement conditional inference forest
#Use CAST package for variable selection









##################################################################################################################################################
############################ MORE ALREADY SPED UP/TRANSCRIBED ####################################################################################
#Set default ranger learner with explicit parameter set
# rangerlrn <- mlr3::lrn('classif.ranger', id = 'ranger', 
#                        num.trees = 500,
#                        #mtry = , #Default is square root of number of variables
#                        min.node.size = 1, #Default is 1 for classification, 5 for regression and 10 for probability
#                        replace = FALSE, #in ranger package, default replace is True but Boulesteix et al. 2012 show that less biased
#                        sample.fraction = 0.632, #Default for sampling without replacement
#                        split.select.weights = ,
#                        always.split.variables= ,
#                        respect.unordered.factors = 'ignore',
#                        importance='permutation',
#                        write.forest=TRUE, 
#                        scale.permutation.importance = FALSE,
#                        num.threads = bigstatsr::nb_cores(),
#                        save.memory = FALSE,
#                        verbose = TRUE,
#                        splitrule = "gini",
#                        num.random.splits = 1L,
#                        keep.inbag = FALSE,
#                        predict_type = 'prob',
#                        max.depth = NULL #Unlimited depth default
# )



#---- Get mDur and mFreq in winter and non-winter periods to assess whether intermittency is due to freezing ----

#Compute 3-month period of minimum temperature
window <- 3
mintemplist<- dcast(
  melt(
    reaches_with_points[, c('GRDC_NO', grep('^tmp.*_c[0-9]{2}', colnames(reaches_with_points), value = T)), with=F], 
    id.var = 'GRDC_NO'),
  variable~GRDC_NO, value.var='value')[,-1] %>% #Pre-format temperature
  frollapply(n=window, FUN=mean, align="center") %>% #Compute 3-month rolling mean
  lapply(function(x) {
    mincenter <- which(x==min(x,na.rm=T)) #Get center month index of 3-month period with minimum temp
    return(seq(mincenter-floor(window/3), mincenter+ceiling(window/3)))}) #Recreate window of months
names(mintemplist) <- reaches_with_points$GRDC_NO

#Compute mdur and mfreq for winter months
gaugestats_winter <- durfreq_parallel(pathlist=fileNames, maxgap=20, monthsel_list=mintemplist, reverse=FALSE)

#Compute mdur and mfreq for non-winter months
gaugestats_nowinter <- durfreq_parallel(pathlist=fileNames, maxgap=20, monthsel_list=mintemplist, reverse=TRUE)

#---- Check month of intermittency ----
monthinter_melt <- melt(gaugestats[, c('GRDC_NO', "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), with=F],
                        id.var='GRDC_NO')
ggplot(monthinter_melt, aes(x=variable, y=value, fill=GRDC_NO)) +
  geom_area(stat='identity')

#Base ranger training
#To adapt
OOBfulldat <- OOBCurve(mod=mod_basic, measures = list(mmce, auc, brier), task = intermittent.task, data = gaugestats_join[!is.na(cly_pc_cav),]) %>%
  setDT %>%
  .[, numtrees := .I]

#Plot it
ggplot(melt(OOBfulldat, id.vars = 'numtrees'), aes(x=numtrees, y=value, color=variable)) + 
  geom_point(size=1, alpha=0.5) + 
  scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.1)) + 
  scale_x_continuous(expand=c(0,0), breaks = seq(0,1000,100)) + 
  theme_bw()