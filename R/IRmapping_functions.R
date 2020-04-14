#Functions for intermittent river analysis

##### -------------------- Utility functions ---------------------------------
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

comp_derivedvar <- function(dt) {
  #---- Inspect and correct -9999 values ----
  print('Inspect and correct -9999 values')
  #check <- riveratlas[cmi_ix_uyr == -9999,] #All have precipitation = 0
  dt[cmi_ix_uyr == -9999, cmi_ix_uyr := 0]
  
  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  dt[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  dt[snw_pc_cmx == -9999, snw_pc_cmx:=0]
  
  #check <- dt[is.na(sgr_dk_rav),]
  
  print('Number of NA values per column')
  colNAs<- dt[, lapply(.SD, function(x) sum(is.na(x) | x==-9999))]
  print(colNAs)
  
  #Convert -9999 values to NA
  for (j in which(sapply(dt,is.numeric))) { #Iterate through numeric column indices
    set(dt,which(dt[[j]]==-9999),j, NA)} #Set those to 0 if -9999
  
  
  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  pre_mcols <- paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)) #Monthly precipitation columns
  dt[, `:=`(pre_mm_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
            pre_mm_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
     .SDcols= pre_mcols] %>% 
    #       min/max monthly catchment precip (while dealing with times when )
    .[, `:=`(pre_mm_cvar= fifelse(pre_mm_cmx==0, 0, pre_mm_cmn/pre_mm_cmx), 
             #min/max monthly watershed discharge
             dis_mm_pvar=fifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx), 
             #min monthly/average yearly watershed discharge
             dis_mm_pvaryr=fifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr), 
             #catchment average elv - watershec average elev
             ele_pc_rel = fifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav))] 
  
  dt[, cmi_ix_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
     .SDcols= paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly cmi
    .[, swc_pc_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
      .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] #Get minimum monthly swc
  
  return(dt)
}

threshold_misclass <- function(i, in_gaugestats) {
  in_gaugestats[!is.na(cly_pc_cav), intermittent_predcat := fifelse(intermittent_predprob>=i, 1, 0)]
  confumat <-in_gaugestats[!is.na(cly_pc_cav), .N, by=c('intermittent', 'intermittent_predcat')]
  outvec <- data.table(
    i, 
    misclas=in_gaugestats[intermittent!=intermittent_predcat,.N]/
      in_gaugestats[,.N],
    sens = confumat[intermittent=='1' & intermittent_predcat==1, N]/
      confumat[intermittent=='1', sum(N)],
    spec  = confumat[intermittent=='0' & intermittent_predcat==0, N]/
      confumat[intermittent=='0', sum(N)])
  return(outvec)
}

weighted_sd <- function(x, w) {
  #Compute weighted standard deviation
  return(sqrt(sum((w) * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)))
}

extract_impperf_nestedrf <- function(nested_rflearner, imp=T, perf=T) {
  #Get variable importance and performance measure for one instance of a resampled rf learner 
  return(c(if (imp) {nested_rflearner$model$learner$importance()}, 
           if (perf) {nested_rflearner$tuning_result$perf}))
}

weighted_vimportance_nestedrf <- function(rfresamp) {
  #Compute weighted mean and standard deviation of variable importance based on
  #accuracy measure (not error measure — could be amended to be based on 1-error)
  
  varnames <- names(rfresamp$learners[[1]]$model$learner$importance())
  
  out_vimportance <- lapply(rfresamp$data$learner, #Extract vimp and perf for each resampling instance
                            extract_impperf_nestedrf) %>% 
    do.call(rbind, .) %>% 
    as.data.table %>%
    melt(id.var='classif.bacc') %>%
    .[, list(imp_wmean = weighted.mean(value, classif.bacc), #Compute weighted mean for each variable
             imp_wsd =  weighted_sd(value, classif.bacc)),  #Compute weighted sd for each variable
      by=variable] 
  
  return(out_vimportance)
}

extract_pd_nestedrf <- function(learner_id, in_rftuned, datdf, selcols, ngrid) {
  "Create a plot matrix of partial dependence interactions"
  
  in_mod <- in_rftuned$learners[[learner_id]] 
  in_fit <- in_mod$learner$model 
  foldperf <- extract_impperf_nestedrf(in_mod, imp=F, perf=T)
  
  # selcols <- in_vimp_plot$data %>% #Can use that if extracting from tunredrf is expensive
  #   setorder(-imp_wmean) %>%
  #   .[colnums, variable]
  
  pdout <- edarf::partial_dependence(in_fit, vars = selcols, n = ngrid,
                                     interaction = TRUE, data = datdf) %>% #Warning: does not work with data_table
    setDT %>%
    .[,(names(foldperf)) := foldperf]
  
  return(pdout)
}

fread_cols <- function(file_name, cols_tokeep) {
  header <- fread(file_name, nrows = 1, header = FALSE)
  keptcols <- cols_tokeep[cols_tokeep %chin% unlist(header)]
  missingcols <- cols_tokeep[!(cols_tokeep %chin% unlist(header))]
  paste('Importing', file_name, 'with ', length(keptcols), 
        'columns out of ', length(cols_tokeep), 'supplied column names')
  fread(input=file_name, header=TRUE, select=keptcols, verbose=TRUE)
}

get_oversamp_ratio <- function(in_task) {
  #When given an mlr3 task with binary classification target, gets which class is minority and the ratio
  return(
    in_task$data()[, .N,
                   by=get(rftuned$task_inter$target_names)] %>%
      setorder(N) %>%
      .[, list(minoclass=get[1], ratio=N[2]/N[1])]
  )
}

po_over <- po("classbalancing", id = "oversample", adjust = "minor", 
              reference = "minor", shuffle = FALSE, 
              ratio = get_oversamp_ratio(rftuned$task_inter)$ratio)
table(po_over$train(list(rftuned$task_inter))$output$truth()) 

#Not used
durfreq_parallel <- function(pathlist, maxgap, monthsel_list=NULL, 
                             reverse=FALSE) {
  #Run durfreq_indiv in parallel
  #Can compute mean number of zero flow days, frequency of zero flow periods, and intermittency on selected months
  #monthsel_list: named list with names being GRDC_NO and values the selected months
  #reverse: whether to use monthsel_list as excluding months rather than subsetting
  
  cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)
  
  return(as.data.table(rbindlist(foreach(j=pathlist, 
                                         .packages = c("data.table", "magrittr"), 
                                         .export=c('comp_durfreq', 'zero.lomf', 'diny')) 
                                 %dopar% {
                                   if (is.null(monthsel_list)) {
                                     return(comp_durfreq(j, maxgap = 20))
                                   } else {
                                     monthsublist = monthsel_list[[strsplit(basename(j), '[.]')[[1]][1]]]
                                     
                                     #Run duration and frequency computation for subset of months
                                     return(comp_durfreq(j, maxgap = maxgap, 
                                                         monthsel= if (reverse) setdiff(1:12, monthsublist) else monthsublist)
                                     )
                                   }
                                 }
  )))
}

##### -------------------- Workflow functions ---------------------------------
def_filestructure <- function() {
  # Get main directory for project
  rootdir <- find_root(has_dir("src")) 
  #Directory where raw data are located
  datdir <-  file.path(rootdir, 'data')
  # Directory where results and figures are written
  resdir <- file.path(rootdir, 'results')
  # File geodatabase to write outputs
  outgdb <- file.path(resdir, 'spatialoutputs.gdb')
  # hydrometric stations that have been joined to RiverATLAS
  in_gaugep <- file.path(outgdb, 'GRDCstations_riverjoin') 
  # Directory containing hydrometric data
  in_gaugedir <-  file.path(datdir, 'GRDCdat_day')
  # River atlas formatted variables
  in_riveratlas_meta <- file.path(datdir, 'HydroATLAS', 'HydroATLAS_metadata_MLM.xlsx')
  # Output geopackage of hydrometric stations with appended predicted intermittency class
  out_gauge <- file.path(resdir, 'GRDCstations_predbasic800.gpkg') 
  # River atlas attribute data1
  in_riveratlas <- file.path(resdir, 'RiverATLAS_v10tab.csv') 
  
  return(c(rootdir=rootdir, datdir=datdir, resdir=resdir, outgdb=outgdb, 
           in_gaugep=in_gaugep, in_gaugedir=in_gaugedir, in_riveratlas_meta=in_riveratlas_meta, 
           out_gauge=out_gauge, in_riveratlas=in_riveratlas))
}

read_gaugep <- function(in_filestructure, dist) {
  #Import gauge stations and only keep those < dist m from a HydroSHEDS reach
  return(st_read(dsn=dirname(in_filestructure['in_gaugep']),
                 layer=basename(in_filestructure['in_gaugep'])) %>%
           .[.$station_river_distance<dist,]
  )
}

read_gauged_paths <- function(in_filestructure, in_gaugep) {
  #Get data paths of daily records for gauge stations
  fileNames <- file.path(in_filestructure['in_gaugedir'], paste(in_gaugep$GRDC_NO,  ".txt", sep=""))
  #Check whether any GRDC record does not exist
  print(paste(length(which(do.call(rbind, lapply(fileNames, file.exists)))),
              'GRDC records do not exist...'))
  return(fileNames)
}

comp_durfreq <- function(path, maxgap, monthsel=NULL) {
  #For a given gauge discharge file (given by path), format and compute 
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
  
  gaugetab <- cbind(fread(path, header=T, skip = 40, sep=";", 
                          colClasses=c('character', 'character', 'numeric', 'numeric', 'integer')),
                    GRDC_NO = strsplit(basename(path), '[.]')[[1]][1])%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)
  
  gaugetab[, year:=as.numeric(substr(dates, 1, 4))] %>% #Create year column and total number of years on record
    .[, prevflowdate := gaugetab[zero.lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA]
  
  #Compute number of missing days per year
  gaugetab[!(Original == -999 | is.na(Original)), `:=`(missingdays = diny(year)-.N,
                                                       datadays = .N), by= 'year'] 
  
  
  gaugetab[,month:=as.numeric(substr(dates, 6, 7))] #create a month column
  
  if (gaugetab[(missingdays < maxgap), .N>0]) { #Make sure that there are years with sufficient data
    monthlyfreq <- gaugetab[
      (missingdays < maxgap), 
      .(GRDC_NO = unique(GRDC_NO),
        monthrelfreq = length(unique(na.omit(prevflowdate)))/length(unique(year))), by='month'] %>% #Count proportion of years with zero flow occurrence per month
      dcast(GRDC_NO~month, value.var='monthrelfreq') %>%
      setnames(as.character(seq(1,12)),
               c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  } else {
    monthlyfreq <- data.table(GRDC_NO = gaugetab[, unique(GRDC_NO)], month=1:12, monthrelfreq=rep(NA,12)) %>%
      dcast(GRDC_NO~month, value.var='monthrelfreq') %>%
      setnames(as.character(seq(1,12)),
               c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  }
  
  #If analysis is only performed on a subset of months
  if (!is.null(monthsel)) {
    gaugetab <- gaugetab[month %in% monthsel, ]
  }
  #Compute number of days of zero flow for years with number of gap days under threshold
  gaugetab_yearly <- merge(gaugetab[, .(missingdays=max(missingdays, na.rm=T), 
                                        datadays=max(datadays, na.rm=T)), by='year'],
                           gaugetab[Original == 0, .(dur=.N,
                                                     freq=length(unique(prevflowdate))), by='year'], 
                           by = 'year', all.x = T) %>%
    .[!is.finite(missingdays), missingdays := diny(year)] %>%
    .[!is.finite(datadays), datadays := diny(year)] %>%
    .[dur==diny(year), freq:=1] %>%
    .[is.na(dur), `:=`(dur=0, freq=0)]
  
  
  gaugetab_all <- cbind(monthlyfreq,
                        gaugetab_yearly[, .(firstYear=min(year), 
                                            lastYear=max(year), 
                                            totalYears=length(unique(year)))],
                        gaugetab_yearly[missingdays <= maxgap,
                                        .(firstYear_kept=min(year), 
                                          lastYear_kept=max(year), 
                                          totalYears_kept=length(unique(year)),
                                          totaldays = sum(datadays), 
                                          sumDur = sum(dur),
                                          mDur = mean(dur),
                                          mFreq = mean(freq),
                                          intermittent = factor(fifelse(sum(dur)>0, 1, 0), levels=c('0','1')))]
  )
  return(gaugetab_all)
}

format_gaugestats <- function(in_gaugestats, in_gaugep) {
  #Format gaugestats into data.table
  
  #Join intermittency statistics to predictor variables and subset to only include those gauges with at least
  gaugestats_join <- do.call(rbind, in_gaugestats) %>%
    setDT %>%
    .[as.data.table(in_gaugep), on='GRDC_NO'] %>%
    .[!is.na(totalYears_kept) & totalYears_kept>10,] %>% # Only keep stations with at least 10 years of data
    comp_derivedvar #Compute derived variables and remove -9999
  return(gaugestats_join)
}

selectformat_predvars <- function(in_filestructure, in_gaugestats) {
  #---- List predictor variables ----
  predcols<- c('dis_m3_pyr',
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
  message(paste(length(predcols[!(predcols %in% names(in_gaugestats))]), 
                'variables are missing from formatted gauge dataset')) 
  
  #---- Associate HydroATLAS column names with variables names ----
  
  #Get predictor variable names
  metaall <- read.xlsx(in_filestructure['in_riveratlas_meta'], 
                       sheetName='Overall') %>%
    setDT 
  
  metascale <- read.xlsx(in_filestructure['in_riveratlas_meta'], 
                         sheetName='scale') %>%
    setDT %>% 
    setnames('Key', 'Keyscale')
  
  metastat <- read.xlsx(in_filestructure['in_riveratlas_meta'], 
                        sheetName='stat') %>%
    setDT %>% 
    setnames('Key', 'Keystat')
  
  meta_format <- as.data.table(expand.grid(Column.s.=metaall$Column.s., 
                                           Keyscale=metascale$Keyscale, 
                                           Keystat=metastat$Keystat)) %>%
                .[metaall, on='Column.s.'] %>%
                .[metascale, on = 'Keyscale'] %>%
                .[metastat, on = 'Keystat',
                  allow.cartesian=TRUE]
  
  meta_format[, `:=`(
    varcode = paste0(gsub('[-]{3}', '', Column.s.),
                     Keyscale, 
                     Keystat),
    varname = paste(Attribute, 
                    Spatial.representation, 
                    Temporal.or.statistical.aggregation.or.other.association))]
  
  #Add newly generated variables to meta_format (variable labels)
  addedvars <- data.table(varname=c('Precipitation catchment Annual min/max', 
                                    'Discharge watershed Annual min/max', 
                                    'Discharge watershed Annual min/average',
                                    'Elevation catchment average - watershed average'), 
                          varcode=c('pre_mm_cvar', 
                                    'dis_mm_pvar', 'dis_mm_pvaryr',
                                    'ele_pc_rel'))
  
  predcols_dt <- merge(data.table(varcode=predcols), 
                       rbind(meta_format, addedvars, fill=T), 
                       by='varcode', all.y=F)
  
  return(predcols_dt)
}

tune_rf <- function(in_gaugestats, in_predvars, 
                    insamp_nfolds, insamp_neval, insamp_nbatch,
                    outsamp_nrep, outsamp_nfolds) {
  #Create task
  datsel <- in_gaugestats[!is.na(cly_pc_cav), 
                          c('intermittent',in_predvars$varcode),
                          with=F]
  task_inter <- mlr3::TaskClassif$new(id ='inter_basic',
                                      backend = datsel,
                                      target = "intermittent")
  
  #Create learner
  lrn_ranger <- mlr3::lrn('classif.ranger', 
                          num.trees=500, 
                          replace=F, 
                          predict_type="prob", 
                          importance = "permutation")
  #print(lrn_ranger$param_set)
  
  #Create mlr3 pipe operator to oversample minority class based on major/minor ratio
  #https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
  #https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
  #Sampling happens only during training phase.
  po_over <- po("classbalancing", id = "oversample", adjust = "minor", 
                reference = "minor", shuffle = TRUE, 
                ratio = get_oversamp_ratio(rftuned$task_inter)$ratio)
  #table(po_over$train(list(rftuned$task_inter))$output$truth()) #Make sure that oversampling worked
  
  #Create a graph learner so that oversampling happens systematically upstream of all training
  lrn_ranger_over <- GraphLearner$new(po_over %>>% lrn_ranger)
  
  #Define parameter ranges to tune on
  tune_ranger <- ParamSet$new(list(
    ParamInt$new("mtry", lower = 1, upper = 11),
    ParamInt$new("min.node.size", lower = 1, upper = 10),
    ParamDbl$new("sample.fraction", lower = 0.1, upper = 0.95)
  ))
  
  #Define inner resampling strategy
  rcv_ranger = rsmp("cv", folds=insamp_nfolds) #5-fold aspatial CV repeated 10 times
  #Define performance measure
  measure_ranger = msr("classif.bacc") #use balanced accuracy to account for imbalanced dataset
  #Define termination rule 
  evals20 = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds
  
  #Define hyperparameter tuner wrapper for inner sampling
  learns = list(
    at_ranger <- AutoTuner$new(learner= lrn_ranger,
                               resampling = rcv_ranger, 
                               measures = measure_ranger,
                               tune_ps = tune_ranger, 
                               terminator = evals20,
                               tuner =  tnr("random_search", 
                                            batch_size = insamp_nbatch)), #batch_size determines level of parallelism
    
    at_ranger_over <- AutoTuner$new(learner= lrn_ranger_over,
                                    resampling = rcv_ranger, 
                                    measures = measure_ranger,
                                    tune_ps = tune_ranger, 
                                    terminator = evals20,
                                    tuner =  tnr("random_search", 
                                                 batch_size = insamp_nbatch)) 
  )

  
  #Perform outer resampling, keeping models for diagnostics later
  outer_resampling = rsmp("repeated_cv", 
                          repeats = outsamp_nrep, 
                          folds = outsamp_nfolds)
  
  nestedresamp_bmrdesign <- benchmark_grid(tasks = task_inter, 
                                                 learners = learns, 
                                                 resamplings = outer_resampling)
  
  nestedresamp_bmrout = benchmark(nestedresamp_bmrdesign,
                                  store_models = TRUE)
  
  #Train full model
  at_ranger$train(task_inter)
  
  print(nestedresamp_ranger$aggregate())
  print(nestedresamp_ranger$prediction()$confusion)
  return(list(rf_outer = nestedresamp_ranger, #Resampling results
              rf_inner = at_ranger, #Core learner (with hyperparameter tuning)
              task_inter = task_inter)) #Task
}

ggvimp <- function(in_rftuned, in_predvars) {
  varimp_basic <- weighted_vimportance_nestedrf(in_rftuned$rf_outer) %>% 
    merge(., in_predvars, by.x='variable', by.y='varcode') %>%
    .[, varname := factor(varname, levels=varname[order(-imp_wmean)])]
  
  ggplot(varimp_basic[1:30,],aes(x=varname)) + 
    geom_bar(aes(y=imp_wmean), stat = 'identity') +
    geom_errorbar(aes(ymin=imp_wmean-2*imp_wsd, ymax=imp_wmean+2*imp_wsd)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    theme_classic() +
    theme(axis.text.x = element_text(size=8))
}

ggmisclass <-  function(in_gaugestats, in_rftuned) {
  #Get predicted probabilities of intermittency for each gauge
  in_gaugestats[!is.na(cly_pc_cav), intermittent_predprob := 
                  as.data.table(in_rftuned$rf_outer$prediction())[order(row_id), mean(prob.1), by=row_id]$V1]
  #Get misclassification error, sensitivity, and specificity for different classification thresholds 
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class
  threshold_confu_dt <- ldply(seq(0,1,0.01), threshold_misclass, in_gaugestats) %>%
    setDT
  
  #Get classification threshold at which sensitivity and specificity are the most similar
  balanced_thresh <- threshold_confu_dt[which.min(abs(spec-sens)),] 
  print(paste('Sensitivity =', round(balanced_thresh$sens,2),
              'and Specificity =', round(balanced_thresh$spec,2),
              'at a classification threshold of', balanced_thresh$i))
  
  #Plot it
  return(
    ggplot(melt(threshold_confu_dt, id.vars='i'), aes(x=i, y=value, color=variable)) + 
      geom_line(size=1.2) + 
      geom_vline(xintercept=balanced_thresh$i) +
      scale_x_continuous(expand=c(0,0)) + 
      scale_y_continuous(expand=c(0,0)) +
      theme_bw()
  )
}

ggpd <- function (in_rftuned, in_predvars, colnums, ngrid, parallel=T) {
  
  #Get partial dependence across all folds
  nlearners <- in_rftuned$rf_outer$resampling$param_set$values$folds
  datdf <- as.data.frame(in_rftuned$rf_outer$data$task[[1]]$data()) #This may be shortened
  varimp <- weighted_vimportance_nestedrf(in_rftuned$rf_outer)
  selcols <- as.character(varimp$variable[colnums])
  
  if (parallel) {
    print(paste("Computing partial dependence with future.apply across", nlearners,
                "CV folds"))
    pd <- future.apply::future_lapply(seq_len(nlearners),
                                      extract_pd_nestedrf,
                                      in_rftuned = in_rftuned$rf_outer,
                                      datdf = datdf,
                                      selcols = selcols,
                                      ngrid = ngrid,
                                      future.scheduling = structure(TRUE,ordering = "random"), 
                                      future.packages = c("data.table","edarf","ranger"))
    
  } else {
    print(paste("Computing partial dependence iteratively across", nlearners,
                "CV folds"))
    pd <- lapply(seq_len(nlearners),
                 extract_pd_nestedrf,
                 in_rftuned = in_rftuned$rf_outer,
                 datdf = datdf,
                 selcols = selcols,
                 ngrid = ngridb)
  }
  
  #Get weighted mean
  pdformat <- do.call(rbind, pd) %>%
    setDT %>%
    .[, list(mean1 = weighted.mean(`1`, classif.bacc)),
      by= selcols]
  
  vargrid <- t(combn(1:length(selcols), 2))
  #leglims <- pdformat[, c(min(mean1), max(mean1))]
  
  #Iterate over every pair of variables
  tileplots_l <- mapply(function(i, j) {
    print(paste0('Generating plot', i , j))
    xvar <- selcols[i]
    yvar <- selcols[j]
    
    pdbivar <- pdformat[, list(bimean = mean(mean1)), by=c(xvar, yvar)]
    
    ggplotGrob(
      ggplot(pdbivar, aes_string(x=xvar, y=yvar)) +
        geom_tile(aes(fill = bimean)) + 
        scale_fill_distiller(palette='Spectral') + # , limits=leglims
        geom_jitter(data=datdf, aes(color=intermittent)) +
        labs(x=stringr::str_wrap(in_predvars[varcode==xvar, varname], 
                                 width = 30),
             y=stringr::str_wrap(in_predvars[varcode==yvar, varname],
                                 width = 30)) +
        theme_bw()
    )
  }, vargrid[,1], vargrid[,2] 
  )
  #Plot, only keeping unique combinations by grabbing lower triangle of plot matrix
  return(do.call("grid.arrange", list(grobs=tileplots_l))) 
}

write_preds <- function(in_filestructure, in_gaugep, in_gaugestats, in_rftuned, predvars) {
  # ---- Output GRDC predictions as points ----
  in_gaugestats[!is.na(cly_pc_cav), 
                IRpredprob := in_rftuned$rf_inner$predict(in_rftuned$task_inter)$prob[,2]]
  
  cols_toditch<- colnames(in_gaugep)[colnames(in_gaugep) != 'GRDC_NO']
  st_write(obj=merge(in_gaugep, 
                     in_gaugestats[, !cols_toditch, with=F], by='GRDC_NO'),
           dsn=in_filestructure['out_gauge'], 
           driver = 'gpkg', 
           delete_dsn=T)
  
  # ---- Generate predictions to map on river network ----
  cols_tokeep <-  c("HYRIV_ID", predvars[!is.na(ID),varcode],
                    'ele_mt_cav','ele_mt_uav',
                    paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
                    paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
                    paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0)))
  
  riveratlas <- fread_cols(in_filestructure['in_riveratlas'], 
                           cols_tokeep = cols_tokeep) %>%
    comp_derivedvar
  
  #Predict model (should take ~10-15 minutes for 9M rows x 40 cols — chunk it up by climate zone to avoid memory errors)
  for (clz in unique(riveratlas$clz_cl_cmj)) {
    print(clz)
    tic()
    riveratlas[!is.na(cly_pc_cav) & !is.na(cly_pc_uav) & clz_cl_cmj == clz,
               predbasic800 := in_rftuned$rf_inner$predict_newdata(as.data.frame(.SD))$prob[,2],]
    toc()
  }
  
  riveratlas[, predbasic800cat := ifelse(predbasic800>=0.5, 1, 0)]
  fwrite(riveratlas[, c('HYRIV_ID', 'predbasic800', 'predbasic800cat'), with=F],
         file.path(in_filestructure['resdir'], 'RiverATLAS_predbasic800.csv'))
}

