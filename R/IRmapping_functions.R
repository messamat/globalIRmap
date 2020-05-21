#Functions for intermittent river analysis

##### -------------------- Utility functions ---------------------------------

#------ diny -----------------
#' Number of days in the year
#'
#' Computes the number of days in a given year, taking in account leap years
#'
#' @param year Numeric or integer vector of the year (format: 4-digit year, %Y).
#'   
#' @return Number of days in that year.
#'
#' @examples
#' ```{r}
#' diny(1999)
#' diny(2000)
#' diny(2004)
#' diny(2100)
#' diny(1600)
#' ```
#' @export
diny <- function(year) {
  365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
}

#------ zero_lomf -----------------
#' Last \[non-zero\] Observation Moved Forward (lomf)
#'
#' Finds the index, for each row, of the previous row with a non-zero value
#'
#' @param x Numeric vector.
#' @param first (logical) Whether to consider first value as a non-zero value
#'   whose index is moved forward even if it is zero. This prevents having NAs
#'   in the results and somewhat assumes that, for a time series, the day prior
#'   to the first value is non-zero.
#'
#' @return Numeric vector of the indices of the previous non-zero for each
#'   element of the input vector.
#'
#' @examples
#' ```{r}
#' test1 <- c(1,1,1,0,0,0,0,1,1)
#' zero_lomf(test)
#' test2 <- c(0,0,0,0,0,1,1,0,1)
#' zero_lomf(test2, first=FALSE)
#' zero_lomf(test2, first=TRUE)
#' ```
#' @export
zero_lomf <- function(x, first=TRUE) {
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if(first==T & x[1]==0)
      non.zero.idx=c(1,non.zero.idx)
    #Repeat index of previous row with non-zero as many times gap until next non-zero values
    rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1))) 
  }
}

#------ comp_derivedvar -----------------
#' Compute derived variables
#'
#' Format and compute derived a set of environmental variables from input
#' dataset that contains
#' \href{https://www.hydrosheds.org/page/hydroatlas}{HydroATLAS} variables for
#' the purpose of predicting intermittency.
#'
#' @param dt data.table which contains rows as records and HydroATLAS
#'   environmental variables as columns.
#' @param copy (logical) whether to make a copy of \code{dt} (if TRUE) or to
#'   modify it in place (if FALSE)
#'
#' @details This function only formats and compute derived variables deemed
#'   relevant for the intermittency analysis. \n
#'   Steps include: \n
#'   1. Convert some -9999 that should be 0 after checking them
#'   2. Count the number of -9999 values per column
#'   3. Convert the rest of the -9999 values to NA
#'   4. Compute a set of derived variables based on existing variables in 
#'   HydroATLAS 
#'   (e.g. \code{pre_mm_cvar= fifelse(pre_mm_cmx==0, 0, pre_mm_cmn/pre_mm_cmx)})
#'   5. Correct scaling from RiverATLAS 
#'   (e.g. Degree of regulation is out of 10,000 in HydroATLAS)
#'
#' ```
#' @source Linke, S., Lehner, B., Dallaire, C. O., Ariwi, J., Grill, G., Anand,
#'   M., ... & Tan, F. (2019). Global hydro-environmental sub-basin and river
#'   reach characteristics at high spatial resolution. Scientific Data, 6(1),
#'   1-15.
#'   
#' @export

comp_derivedvar <- function(in_dt, copy=FALSE) {
  if (copy) {
    in_dt2 <- copy(in_dt)
  } else {
    in_dt2 <- in_dt
  }
  
  #---- Inspect and correct -9999 values ----
  print('Inspect and correct -9999 values')
  #check <- riveratlas[cmi_ix_uyr == -9999,] #All have precipitation = 0
  in_dt2[cmi_ix_uyr == -9999, cmi_ix_uyr := 0]
  
  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  in_dt2[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  in_dt2[snw_pc_cmx == -9999, snw_pc_cmx:=0]
  
  #check <- in_dt2[is.na(sgr_dk_rav),]
  
  print('Number of NA values per column')
  colNAs<- in_dt2[, lapply(.SD, function(x) sum(is.na(x) | x==-9999))]
  print(colNAs)
  
  #Convert -9999 values to NA
  for (j in which(sapply(in_dt2,is.numeric))) { #Iterate through numeric column indices
    set(in_dt2,which(in_dt2[[j]]==-9999),j, NA)} #Set those to 0 if -9999
  
  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  pre_mcols <- paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)) #Monthly precipitation columns
  in_dt2[, `:=`(pre_mm_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
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
  
  in_dt2[, cmi_ix_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
     .SDcols= paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly cmi
    .[, swc_pc_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
      .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] #Get minimum monthly swc
  
  
  #Scale variables based on HydroATLAS documentation
  in_dt2[, `:=`(
    ari_ix_cav = ari_ix_cav/100,
    ari_ix_uav = ari_ix_uav/100,
    cmi_ix_cmn = cmi_ix_cmn/100,
    cmi_ix_uyr = cmi_ix_uyr/100,
    dor_pc_pva = dor_pc_pva/100,
    lka_pc_cse = lka_pc_cse/10,
    lka_pc_use = lka_pc_use/10,
    tmp_dc_cmn =  tmp_dc_cmn/10,
    tmp_dc_cmx =  tmp_dc_cmx/10,
    tmp_dc_cyr =  tmp_dc_cyr/10,
    tmp_dc_uyr =  tmp_dc_uyr/10,
    gwt_m_cav = gwt_cm_cav/100
  )]
  
  return(in_dt2)
}

threshold_misclass <- function(i, in_preds) {
  if (inherits(in_preds, 'PredictionClassif')) {
    confu <- as.data.table(in_preds$set_threshold(1-i)$confusion)
  } 
  
  if (is.data.table(in_preds)) {
    if (in_preds[task_type == 'classif',.N] > 0) {
      confu <- in_preds[, as.data.table(pred[[1]]$set_threshold(1-i)$confusion), #Threshold is based on prob.0, not prob.1, so need to use 1-i
                        by=outf] %>% 
        .[, .(N=sum(N)), by=.(response, truth)]
    } 
    
    if (in_preds[task_type == 'regr', .N] > 0) {
      confu <- in_preds[, response := fifelse(prob.1>=i, '1', '0')] %>% 
        .[, truth := as.character(truth)] %>%
        .[, .N, by=.(response, truth)]
    }
  }
  
  outvec <- data.table(
    i, 
    misclas = confu[truth != response, sum(N)] / confu[, sum(N)],
    sens = confu[truth == '1' & response == '1', N] / confu[truth=='1', sum(N)],
    spec  = confu[truth=='0' & response==0, N]/confu[truth=='0', sum(N)])
  return(outvec)
}


get_outerrsmp <- function(in_rftuned, spatial_rsp=FALSE) {
  #Adapt whether direct resample result or output from selecttrain_rf
  if (inherits(in_rftuned, 'list')) {
    #If there is more than one type of outer resampling
    if (length(in_rftuned$rf_outer) > 0) {
      #Check which resampling is spatial — only works if one is spatial
      sp_i <- which(unlist(lapply(in_rftuned$rf_outer, function(x) {
        grepl('.*Resampling.*Sp.*', x$resampling$format())
      })))
      
      #If user request that spatial resampling be used
      if (spatial_rsp==TRUE) {
        
        if (length(sp_i)>0) {
          rsmp_res <- in_rftuned$rf_outer[[sp_i]] 
        } else { #But if there is no spatial resampling provided
          stop("spatial_rsp==TRUE but the in_rftuned does not include 
               any Spatial Resampling")
        }
        #If user didn't request spatial resampling to be used, grab the first 
        #resampling that is not spatial
      } else {
        rsmp_res <- in_rftuned$rf_outer[[
          min((1:length(in_rftuned$rf_outer))[-sp_i])]]  
      }
      
      #If there is only one type of outer resampling
    } else {
      rsmp_res <- in_rftuned$rf_outer
    }
  } else if (inherits(in_rftuned, "ResampleResult")) {
    rsmp_res <- in_rftuned 
  }
  return(rsmp_res)
}


weighted_sd <- function(x, w) {
  #Compute weighted standard deviation
  return(sqrt(sum((w) * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)))
}

extract_impperf_nestedrf <- function(nested_rflearner, 
                                     imp = TRUE, perf = TRUE, 
                                     pvalue = TRUE, pvalue_permutn = 100) {
  #Get variable importance and performance measure for one instance of a resampled rf learner 
  if (inherits(nested_rflearner, "AutoTuner")) {
    sublrn <- nested_rflearner$model$learner
  } else {
    sublrn <- nested_rflearner
  }
  
  print(paste0("Computing variable importance p_values", 
           "for resampling instance hash #", sublrn$hash))
  
  return(cbind(if (imp) {
    ####################### IF GraphLearner ####################################
    if (inherits(sublrn, "GraphLearner")) { 
      
      if ('classif.ranger' %in% names(sublrn$model)) {
        
        if (pvalue == TRUE) {
          in_task <- nested_rflearner$model$tuning_instance$task
          in_formula <- as.formula(paste0(in_task$target_names, '~.'))
          
          importance_pvalues(
            sublrn$model$classif.ranger$model,
            method = "altmann",
            num.permutations = pvalue_permutn,
            data = in_task$data(),
            formula= in_formula
          )
          
        } else {
          data.table(importance=sublrn$model$classif.ranger$model$variable.importance)
        }
      } 
      
      else if ('classif.cforest' %in% names(sublrn$model)) {
        if (pvalue == TRUE) {
          warning("p_value calculation is only available for ranger classification rf, ignoring p_value")
        }
        data.table(importance= 
                     partykit::varimp(sublrn$model$classif.cforest$model,
                                      nperm = 1L, 
                                      OOB = TRUE, 
                                      risk = c("loglik", "misclassification"),
                                      conditional = FALSE, 
                                      threshold = .2))
      }
    } else { ####################### IF direct model ####################################
      if (pvalue == TRUE) { #If want pvalue associated with predictor variables
        in_task <- nested_rflearner$model$task
        in_formula <- as.formula(paste0(in_task$target_names, '~.'))
        
        importance_pvalues(
          nested_rflearner$model,
          method = "altmann",
          num.permutations = pvalue_permutn,
          data = in_task$data(),
          formula= in_formula
        )
      } else { #If pvalue == FALSE
        data.table(importance= nested_rflearner$model$learner$importance())
      }
      
    }
  },
  if (perf) {
    outperf <- nested_rflearner$tuning_result$perf
    data.table(outperf) %>% setnames(names(outperf))
  }
  ))
}


weighted_vimportance_nestedrf <- function(rfresamp,
                                          pvalue = TRUE, pvalue_permutn = 100) {
  #Compute weighted mean and standard deviation of variable importance based on
  #accuracy measure (not error measure — could be amended to be based on 1-error)
  varnames <- rfresamp$task$feature_names
  
  out_vimportance <- lapply(rfresamp$learners, #Extract vimp and perf for each resampling instance
                            extract_impperf_nestedrf, 
                            imp=T, perf=T, pvalue=pvalue, pvalue_permutn) %>% 
    do.call(rbind, .) %>% 
    cbind(., varnames) %>%
    .[,
      list(imp_wmean = weighted.mean(importance, classif.bacc), #Compute weighted mean
           imp_wsd =  weighted_sd(importance, classif.bacc),
           imp_pvalue = weighted.mean(pvalue, classif.bacc)), #Compute weighted sd
      by=varnames] 
  
  return(out_vimportance)
}

extract_pd_nestedrf <- function(learner_id, in_rftuned, datdf, selcols, ngrid) {
  "Create a plot matrix of partial dependence interactions"
  
  in_mod <- in_rftuned$learners[[learner_id]] 
  #in_mod <- nestedresamp_ranger$learners[[1]] 
  
  if (inherits(in_mod$learner, "GraphLearner")) {
    in_fit <- in_mod$learner$model$classif.ranger$model
  } else {
    in_fit <- in_mod$learner$model 
  }
  
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
    in_task$data()[, .N, by=get(in_task$target_names)] %>%
      setorder(N) %>%
      .[, list(minoclass=get[1], ratio=N[2]/N[1])]
  )
}

convert_clastoregrtask <- function(in_task, in_id, oversample=FALSE) {
  if (oversample) {
    oversamp_ratio <- get_oversamp_ratio(in_task)
    mino_subdat <- in_task$data()[
      eval(in_task$target_names) == oversamp_ratio$minoclass,] #Get part of the underlying dataset of the minority class
    nmino <- mino_subdat[,.N]
    oversamp_index <- sample(nmino, nmino*(oversamp_ratio$ratio-1), replace=T)
    
    newdat <- rbind(in_task$data(),
                    mino_subdat[oversamp_index,]) %>%
      .[, eval(in_task$target_names):= as.numeric(as.character(get(in_task$target_names)))]
  } else {
    newdat <- in_task$data()[, eval(in_task$target_names) := 
                               as.numeric(as.character(get(in_task$target_names)))]
  }
  
  return(mlr3::TaskRegr$new(id =in_id,
                            backend = newdat,
                            target = in_task$target_names))
}

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
                                         .export=c('comp_durfreq', 'zero_lomf', 'diny')) 
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
    .[, prevflowdate := gaugetab[zero_lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
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
                                          intermittent = factor(fifelse(mean(dur)>=1, 1, 0), levels=c('0','1')))]
  )
  return(gaugetab_all)
}

format_gaugestats <- function(in_gaugestats, in_gaugep) {
  #Format gaugestats into data.table
  
  #Join intermittency statistics to predictor variables and subset to only include those gauges with at least
  gaugestats_join <- do.call(rbind, in_gaugestats) %>%
    setDT %>%
    .[as.data.table(in_gaugep), on='GRDC_NO'] %>%
    .[!is.na(totalYears_kept) & totalYears_kept>=10,] %>% # Only keep stations with at least 10 years of data
    .[, c('X', 'Y') := as.data.table(sf::st_coordinates(Shape))] %>%
    comp_derivedvar #Compute derived variables, rescale some variables, remove -9999
  
  return(gaugestats_join)
}

selectformat_predvars <- function(in_filestructure, in_gaugestats) {
  #---- List predictor variables ----
  predcols<- c(
    'ORD_STRA',
    'dis_m3_pyr',
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
    'gwt_m_cav',
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
    'kar_pc_cse',
    
    'ppd_pk_cav',
    'ppd_pk_uav',
    'urb_pc_cse',
    'urb_pc_use',
    'hft_ix_c93',
    'hft_ix_u93',
    'hft_ix_c09',
    'hft_ix_u09',
    'gdp_ud_cav',
    'hdi_ix_cav')
  
  
  
  
  
  #Check that all columns are in dt
  message(paste(length(predcols[!(predcols %in% names(in_gaugestats))]), 
                'variables are missing from formatted gauge dataset')) 
  
  #---- Associate HydroATLAS column names with variables names ----
  
  #Get predictor variable names
  metaall <- readxl::read_xlsx(in_filestructure['in_riveratlas_meta'], 
                       sheet='Overall') %>%
    setDT 
  
  metascale <- readxl::read_xlsx(in_filestructure['in_riveratlas_meta'], 
                         sheet='scale') %>%
    setDT %>% 
    setnames(c('Key','Spatial representation'),
             c('Keyscale', 'Spatial.representation'))
  
  metastat <- readxl::read_xlsx(in_filestructure['in_riveratlas_meta'], 
                        sheet='stat') %>%
    setDT %>% 
    setnames(c('Key','Temporal or statistical aggregation or other association'),
             c('Keystat', 'Temporal.or.statistical.aggregation.or.other.association'))
  
  meta_format <- as.data.table(expand.grid(`Column(s)`=metaall$`Column(s)`, 
                                           Keyscale=metascale$Keyscale, 
                                           Keystat=metastat$Keystat)) %>%
                .[metaall, on='Column(s)'] %>%
                .[metascale, on = 'Keyscale'] %>%
                .[metastat, on = 'Keystat',
                  allow.cartesian=TRUE]
  
  
  
  meta_format[, `:=`(
    unit = substr(`Column(s)`, 5, 6),
    varcode = paste0(gsub('[-]{3}', '', `Column(s)`),
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

create_tasks <- function(in_gaugestats, in_predvars) {
  #Create subset of gauge data for analysis (in this case, remove records with missing soil data)
  datsel <- in_gaugestats[!is.na(cly_pc_cav), 
                          c('intermittent',in_predvars$varcode, 'X', 'Y'),
                          with=F]
  
  #Basic task for classification
  task_classif <- mlr3spatiotempcv::TaskClassifST$new(
    id="inter_basicsp", 
    backend = datsel, 
    target = "intermittent",
    coordinate_names = c("X", "Y"))
  
  #Basic task for regression without oversampling
  task_regr <- convert_clastoregrtask(in_task = task_classif,
                                      in_id = 'inter_regr',
                                      oversample=FALSE) 
  
  #Basic task for regression with oversampling to have the same number of minority and majority class
  task_regrover <- convert_clastoregrtask(in_task = task_classif,
                                          in_id = 'inter_regrover',
                                          oversample=TRUE)
  return(list(classif=task_classif, regr=task_regr, regover=task_regrover))
}

benchmark_classif <- function(in_tasks, 
                              insamp_nfolds, insamp_neval, insamp_nbatch,
                              outsamp_nrep, outsamp_nfolds) {
  
  #Create task vars
  task_classif <- in_tasks$classif
  
  #---------- Create learners --------------------------------------------------
  #Compute ratio of intermittent to perennial observations
  imbalance_ratio <- get_oversamp_ratio(task_classif)$ratio
  
  #Create basic learner
  lrn_ranger <- mlr3::lrn('classif.ranger',
                          num.trees = 500, 
                          replace = FALSE, 
                          splitrule = 'gini',
                          predict_type = 'prob', 
                          importance = 'impurity_corrected',
                          respect.unordered.factors = 'order')
  
  #print(lrn_ranger$param_set)
  
  #Create a conditional inference forest learner with default parameters
  # mtry = sqrt(nvar), fraction = 0.632
  lrn_cforest <- mlr3::lrn('classif.cforest',
                           ntree = 500,
                           alpha = 0.05,
                           replace = F,
                           predict_type = "prob")
  
  #Create mlr3 pipe operator to oversample minority class based on major/minor ratio
  #https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
  #https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
  #Sampling happens only during training phase.
  po_over <- mlr3pipelines::po("classbalancing", id = "oversample", adjust = "minor", 
                              reference = "minor", shuffle = TRUE, 
                              ratio = imbalance_ratio)
  #table(po_over$train(list(task_classif))$output$truth()) #Make sure that oversampling worked
  
  #Create mlr3 pipe operator to put a higher class weight on minority class
  po_classweights <- mlr3pipelines::po("classweights", minor_weight = imbalance_ratio)
  
  #Create graph learners so that oversampling happens systematically upstream of all training
  lrn_ranger_overp <- mlr3pipelines::GraphLearner$new(po_over %>>% lrn_ranger)
  lrn_cforest_overp <- mlr3pipelines::GraphLearner$new(po_over %>>% lrn_cforest)
  
  #Create graph learners so that class weighin happens systematically upstream of all training
  lrn_ranger_weight <- mlr3pipelines::GraphLearner$new(po_classweights %>>% lrn_ranger)
  lrn_cforest_weight <- mlr3pipelines::GraphLearner$new(po_classweights  %>>% lrn_cforest)
  
  #---------- Set up inner resampling ------------------------------------------
  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)
    
    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T), 
                   lower = 5, 
                   upper = floor(length(task_classif$feature_names)/2)), #Half number of features
      ParamDbl$new(grep(".*fraction", prmset, value=T), 
                   lower = 0.2, 
                   upper = 0.8)
    ))
    
    in_split =in_lrn$param_set$get_values()[
      grep(".*split(rule|stat)", prmset, value=T)]
    
    if (in_split == 'maxstat') {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)], 
                     lower = 0.01, upper = 0.1)
      )
      
    } else if (any(grepl(".*min.node.size", prmset))) {
      tune_rf$add(
        ParamInt$new(prmset[grep(".*min.node.size", prmset)], 
                     lower = 1, upper = 10)
      )
    } else if (any(grepl(".*splitstat", prmset))) {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)], 
                     lower = 0.01, upper = 0.1)
      )
    }
  }
  
  #Define inner resampling strategy
  rcv_rf = rsmp("cv", folds=insamp_nfolds) #aspatial CV repeated 10 times
  
  #Define performance measure
  measure_rf_class = msr("classif.bacc") #use balanced accuracy as objective function
  
  #Define termination rule 
  evalsn = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds
  
  #Define hyperparameter tuner wrapper for inner sampling
  learns_classif = list(
    #Standard ranger rf without oversampling
    AutoTuner$new(learner= lrn_ranger,
                  resampling = rcv_rf, 
                  measures = measure_rf_class,
                  tune_ps = regex_tuneset(lrn_ranger), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch)), #batch_size determines level of parallelism
    
    
    #Standard ranger rf with oversampling
    AutoTuner$new(learner= lrn_ranger_overp,
                  resampling = rcv_rf, 
                  measures = measure_rf_class,
                  tune_ps = regex_tuneset(lrn_ranger_overp), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch)),
    
    #Standard ranger rf with class weights
    AutoTuner$new(learner= lrn_ranger_weight,
                  resampling = rcv_rf, 
                  measures = measure_rf_class,
                  tune_ps = regex_tuneset(lrn_ranger_weight), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch)), #batch_size determines level of parallelism
    
    #CIF rf without oversampling
    lrn_cforest,
    
    #CIF rf with oversampling
    lrn_cforest_overp,
    
    #CIF rf with class 
    lrn_cforest_weight
  )
  names(learns_classif) <-mlr3misc::map(learns_classif, "id")
  
  #---------- Set up outer resampling benchmarking -----------------------------
  
  #Perform outer resampling, keeping models for diagnostics later
  outer_resampling = rsmp("repeated_cv", 
                          repeats = outsamp_nrep, 
                          folds = outsamp_nfolds)
  
  #Run outer resampling and benchmarking on classification learners
  nestedresamp_bmrdesign_classif <- benchmark_grid(
    tasks = task_classif, 
    learners = learns_classif, 
    resamplings = outer_resampling)
  
  nestedresamp_bmrout_classif <- benchmark(
    nestedresamp_bmrdesign_classif, store_models = TRUE)
  
  return(
    list(
      bm_classif = nestedresamp_bmrout_classif,
      bm_tasks = list(task_classif=task_classif),
      measure_classif = measure_rf_class
    )
  )
}

benchmark_regr <- function(in_tasks, 
                           insamp_nfolds, insamp_neval, insamp_nbatch,
                           outsamp_nrep, outsamp_nfolds) {
  
  task_regr <- in_tasks$regr
  task_regrover <- in_tasks$regover
  
  #---------- Create learners --------------------------------------------------
  #Create regression learner with maxstat. Represents an approximation of 
  #classification learner with probabilities as the prediction type and a simili-
  #conditional-inference forest tweak to correct for the over-representation
  #of variables with many values over those with few categories
  lrn_ranger_maxstat <- mlr3::lrn('regr.ranger', 
                                  num.trees=500, 
                                  replace=FALSE, 
                                  splitrule = 'maxstat',
                                  importance = 'impurity_corrected',
                                  respect.unordered.factors = 'order')
  
  #---------- Set up inner resampling ------------------------------------------
  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)
    
    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T), 
                   lower = 1, upper = 11),
      ParamDbl$new(grep(".*fraction", prmset, value=T), 
                   lower = 0.2, upper = 0.8)
    ))
    
    in_split =in_lrn$param_set$get_values()[
      grep(".*split(rule|stat)", prmset, value=T)]
    
    if (in_split == 'maxstat') {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)], 
                     lower = 0.01, upper = 0.1)
      )
      
    } else if (any(grepl(".*min.node.size", prmset))) {
      tune_rf$add(
        ParamInt$new(prmset[grep(".*min.node.size", prmset)], 
                     lower = 1, upper = 10)
      )
    } else if (any(grepl(".*splitstat", prmset))) {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)], 
                     lower = 0.01, upper = 0.1)
      )
    }
  }
  
  #Define inner resampling strategy
  rcv_rf = rsmp("cv", folds=insamp_nfolds) #aspatial CV repeated 10 times
  
  #Define performance measure
  measure_rf_reg = msr("regr.mae") 
  
  #Define termination rule 
  evalsn = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds
  
  #Define hyperparameter tuner wrapper for inner sampling
  learns_regr = list(
    #Regression rf with MAXSTAT
    AutoTuner$new(learner= lrn_ranger_maxstat,
                  resampling = rcv_rf, 
                  measures = measure_rf_reg,
                  tune_ps = regex_tuneset(lrn_ranger_maxstat), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch))
  )
  
  #---------- Set up outer resampling benchmarking -----------------------------
  
  #Perform outer resampling, keeping models for diagnostics later
  outer_resampling = rsmp("repeated_cv", 
                          repeats = outsamp_nrep, 
                          folds = outsamp_nfolds)
  #Run outer resampling and benchmarking on regression learners
  nestedresamp_bmrdesign_regr <- benchmark_grid(
    tasks = list(task_regr, task_regrover), 
    learners = learns_regr, 
    resamplings = outer_resampling)
  
  nestedresamp_bmrout_regr <- benchmark(
    nestedresamp_bmrdesign_regr, store_models = TRUE)
  
  return(
    list(
      bm_regr = nestedresamp_bmrout_regr,
      bm_tasks = list(task_regr=task_regr,
                      task_regrover=task_regrover),
      measure_regr = measure_rf_reg
    )
  )
}


analyze_benchmark <- function(in_bm, in_measure) {
  
  print(in_bm$aggregate(in_measure))
  boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure)
  
  print(paste('It took', 
              in_bm$aggregate(msr('time_both'))$time_both,
              'seconds to train and predict with the',
              in_bm$aggregate(msr('time_both'))$learner_id,
              'model...'))
  
  bmdt <- as.data.table(in_bm)
  
  if (in_bm$task_type == 'regr') {
    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- bmdt$prediction[[rsmp_i]]$test %>%
        as.data.table %>%
        .[, `:=`(outf = bmdt$iteration[[rsmp_i]],
                 task = bmdt$task[[rsmp_i]]$id,
                 task_type = in_bm$task_type,
                 learner = bmdt$learner[[rsmp_i]]$id)]
      return(preds)
    }) %>%
      do.call(rbind, .) 
    
    if (!('prob.1' %in% names(preds)) & 'response' %in% names(preds)) {
      preds[, prob.1 := response]
    }
  }
  
  
  if (in_bm$task_type == 'classif') {
    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- data.table(outf = bmdt$iteration[[rsmp_i]],
                          task = bmdt$task[[rsmp_i]]$id,
                          learner = bmdt$learner[[rsmp_i]]$id,
                          task_type = in_bm$task_type,
                          pred = list(bmdt$prediction[[rsmp_i]]$test))
      return(preds)
    }) %>%
      do.call(rbind, .) 
  }
  
  tasklearner_unique <- preds[, expand.grid(unique(task), unique(learner))] %>%
    setnames(c('task', 'learner'))
  
  glist <- lapply(1:nrow(tasklearner_unique), function(tsklrn) {
    print(tasklearner_unique[tsklrn,])
    subpred <- preds[task ==tasklearner_unique$task[tsklrn] &
                       learner == tasklearner_unique$learner[tsklrn],]
    return(ggplotGrob(
      ggmisclass(in_predictions = subpred) + 
        ggtitle(paste(tasklearner_unique$task[tsklrn], 
                      tasklearner_unique$learner[tsklrn]))
    )
    ) 
  })
  
  return(list(bm_misclasscomp=do.call("grid.arrange", list(grobs=glist)),
              bm_boxcomp = boxcomp))
} 

benchmark_featsel <- function(in_rf, in_task, in_measure,
                              pcutoff = 0.1,
                              insamp_nfolds =  NULL, insamp_nevals = NULL,
                              outsamp_nrep = NULL, outsamp_nfolds =  NULL,
                              outsamp_nfolds_sp =  NULL) {
  #pcutoff is the p_value above which features are removes
  
  #Apply feature/variable selection
  vimp <- weighted_vimportance_nestedrf(
    rfresamp = in_rf$resample_result(
      uhash=unique(as.data.table(in_rf)$uhash))) %>%
    .[,imp_wmeanper := imp_wmean/sum(imp_wmean)]
  
  task_featsel <- in_task$clone()$select(
    vimp[imp_pvalue <= pcutoff, as.character(varnames)])
  task_featsel$id <- paste0(in_task$id, '_featsel')
  
  if (is.null(outsamp_nfolds_sp)) {
    outsamp_nfolds_sp = outsamp_nfolds
  }
  
  #Set up outer resampling including 
  outer_resampling = rsmp("repeated_cv", 
                          repeats = outsamp_nrep, 
                          folds = outsamp_nfolds)
  
  outer_resamplingsp = rsmp("repeated-spcv-coords",
                            repeats = outsamp_nrep,
                            folds = outsamp_nfolds_sp) #Create 20 folds (visualize)
  
  #Run outer resampling and benchmarking on classification learners
  bmrdesign_featsel <- benchmark_grid(
    tasks = list(in_task, task_featsel), 
    learners = in_rf$learners$learner[[1]], 
    resamplings = list(outer_resampling, outer_resamplingsp))
  
  bmrout_featsel <- benchmark(
    bmrdesign_featsel, store_models = TRUE)
  
  bm_analysis <- analyze_benchmark(bmrout_featsel, in_measure = in_measure)
  
  return(list(
    bm_classif = bmrout_featsel,
    bm_tasks = list(task_classif=in_task, 
                    task_classif_featsel = task_featsel),
    measure_classif = in_measure,
    bm_analysis = bm_analysis)
  )
}


selecttrain_rf <- function(in_rf, in_task, 
                           insamp_nfolds =  NULL, insamp_nevals = NULL) {
  #Prepare autotuner for full training
  lrn_autotuner <- in_rf$learners$learner[[1]]
  subbm <- in_rf$clone()$filter(task_ids = in_task$id)
  
  if (!is.null(insamp_nfolds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- insamp_nfolds
  }
  
  if (!is.null(insamp_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- insamp_nevals
  }
  
  #Train learners
  lrn_autotuner$train(in_task)
  
  #Return outer sampling object for selected model (or list of outer sampling objects)
  uhashes <- unique(as.data.table(subbm)$uhash)
  if (length(uhashes) == 1) {
    outer_resampling_output <- subbm$resample_result(uhash=uhashes)
  } else {
    outer_resampling_output <- lapply(uhashes, function(x) {
      subbm$resample_result(uhash=x)
    }) 
  }

  return(list(rf_outer = outer_resampling_output, #Resampling results
              rf_inner = lrn_autotuner, #Core learner (with hyperparameter tuning)
              task = in_task)) #Task
}

ggvimp <- function(in_rftuned, in_predvars, spatial_rsp=FALSE) {
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)

  #Get variable importance and format them
  varimp_basic <- weighted_vimportance_nestedrf(rfresamp = rsmp_res,
                                                pvalue = FALSE) %>% 
    merge(., in_predvars, by.x='variable', by.y='varcode') %>%
    .[, varname := factor(varname, levels=varname[order(-imp_wmean)])]
  
  #Plot 'em
  ggplot(varimp_basic[1:30,],aes(x=varname)) + 
    geom_bar(aes(y=imp_wmean), stat = 'identity') +
    geom_errorbar(aes(ymin=imp_wmean-2*imp_wsd, ymax=imp_wmean+2*imp_wsd)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    theme_classic() +
    theme(axis.text.x = element_text(size=8))
}

ggmisclass <-  function(in_predictions) {
  #Get predicted probabilities of intermittency for each gauge
  # in_gaugestats[!is.na(cly_pc_cav), intermittent_predprob := 
  #                 as.data.table(in_predictions)[order(row_id), mean(prob.1), by=row_id]$V1]
  #Get misclassification error, sensitivity, and specificity for different classification thresholds 
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class

  threshold_confu_dt <- ldply(seq(0,1,0.01), threshold_misclass, in_predictions) %>%
    setDT
  
  #Get classification threshold at which sensitivity and specificity are the most similar
  balanced_thresh <- threshold_confu_dt[which.min(abs(spec-sens)),] 
  print(paste('Sensitivity =', round(balanced_thresh$sens,2),
              'and Specificity =', round(balanced_thresh$spec,2),
              'at a classification threshold of', balanced_thresh$i))
  
  gout <- ggplot(melt(threshold_confu_dt, id.vars='i'), aes(x=i, y=value, color=variable)) + 
    geom_line(size=1.2) + 
    geom_vline(xintercept=balanced_thresh$i) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) +
    theme_bw()

  #Plot it
  return(gout)
}

ggpd <- function (in_rftuned, in_predvars, colnums, ngrid, 
                  parallel=T, spatial_rsp=FALSE) {
  
  #Get outer resampling of interest
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
  
  #Get partial dependence across all folds
  nlearners <- rsmp_res$resampling$param_set$values$folds
  datdf <- as.data.frame(rsmp_res$task$data()) #This may be shortened
  varimp <- weighted_vimportance_nestedrf(rsmp_res)
  selcols <- as.character(varimp$variable[colnums])
  
  if (parallel) {
    print(paste("Computing partial dependence with future.apply across", nlearners,
                "CV folds"))
    pd <- future.apply::future_lapply(seq_len(nlearners),
                                      extract_pd_nestedrf,
                                      in_rftuned = rsmp_res,
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
                 in_rftuned = rsmp_res,
                 datdf = datdf,
                 selcols = selcols,
                 ngrid = ngrid)
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
                                 width = 20),
             y=stringr::str_wrap(in_predvars[varcode==yvar, varname],
                                 width = 20)) +
        theme_bw() +
        theme(text = element_text(size=12))
    )
  }, vargrid[,1], vargrid[,2] 
  )
  #Plot, only keeping unique combinations by grabbing lower triangle of plot matrix
  return(do.call("grid.arrange", list(grobs=tileplots_l))) 
}

gguncertainty <- function(in_rftuned, in_gaugestats, in_predvars) {
  #Get average predictions for oversampled rows
  gaugepred <-  in_rftuned$rf_outer$prediction() %>%
    as.data.table %>%
    .[, list(truth=first(truth), prob.1=mean(prob.1)), by=row_id] %>%
    setorder(row_id)
  
  predattri <- cbind(in_gaugestats[!is.na(cly_pc_cav),], gaugepred) %>%
    .[, `:=`(preduncert = prob.1-as.numeric(as.character(intermittent)),
             yearskeptratio = totalYears_kept/totalYears)]
  
  #Plot numeric variables
  predmelt_num <- predattri[, which(as.vector(unlist(lapply(predattri, is.numeric)))), with=F] %>%
    cbind(predattri[, c('GRDC_NO', 'intermittent'), with=F]) %>%
    melt(id.vars=c('GRDC_NO', 'intermittent', 'prob.1', 'preduncert'))
  
  
  #Set variable labels
  varlabels <- setkey(in_predvars[, .(varcode, varname)], varcode)[
    levels(predmelt_num$variable)] %>%
    .[is.na(varname), varname := varcode]
  predmelt_num[, variable := factor(variable, labels = varlabels$varname)]
  
  varstoplot <- merge(in_predvars[, .(varcode, varname)], 
                      data.table(varcode = 
                                   c('totalYears_kept', 'yearskeptratio', 
                                     'mDur', 'mFreq', 'station_river_distance', 
                                     'UPLAND_SKM', 'ORD_STRA', 
                                     'dis_m3_pyr', 'dor_pc_pva',
                                     'cmi_ix_uyr','ari_ix_uav')), 
                      all.x =F, all.y=T)
  plotdt <- predmelt_num[variable %in% varstoplot$varcode | 
                           variable %in% varstoplot$varname ,]
  
  uncertainty_numplot <- 
    ggplot(plotdt, aes(x=value, y=preduncert, color=intermittent)) + 
    geom_rect(xmin=-Inf, xmax=Inf, ymin=-0.5, ymax=0.5, 
              fill='#d9d9d9', color='#d9d9d9') +
    geom_point(alpha = 1/4) + 
    geom_hline(yintercept=0, alpha=1/2) +
    geom_smooth(method='gam', formula = y ~ s(x, k=3)) +
    annotate("text", x = Inf-5, y = 0.5, angle = 90, 
             label = "Pred:Int, Obs:Per",
             color = '#1f78b4') +
    annotate("text", x = Inf-5, y = -0.5, angle = 90, 
             label = "Pred:Per, Obs:Int",
             color = '#ff7f00') + 
    scale_color_manual(values=c('#1f78b4', '#ff7f00'),
                       name='Observed regime', 
                       labels = c('Perennial', 'Intermittent')) + 
    #scale_x_sqrt(expand=c(0,0)) +
    coord_cartesian(clip='off') + 
    facet_wrap(~variable, scales='free', labeller=label_value) + 
    theme_classic() +
    theme(legend.position = c(0.8, 0.1))
  
  
  #Plot numeric variables
  predmelt_cat <- predattri[, c('GRDC_NO', 'intermittent', 'preduncert',
                                'ENDORHEIC', 'clz_cl_cmj'), with=F] %>%
    melt(id.vars=c('GRDC_NO', 'intermittent', 'preduncert'))
  
  uncertainty_catplot <- 
    ggplot(predmelt_cat, aes(x=as.factor(value), y=preduncert, 
                             fill=intermittent, color=intermittent)) + 
    geom_rect(xmin=-Inf, xmax=Inf, ymin=-0.5, ymax=0.5, 
              fill='#d9d9d9', color='#d9d9d9', alpha=1/2) +
    #geom_boxplot(alpha = 0.75) +
    geom_violin(alpha=0.75, color=NA) +
    geom_hline(yintercept=0, alpha=1/2) +
    coord_cartesian(clip='off') + 
    scale_fill_manual(values=c('#1f78b4', '#ff7f00'),
                      name='Observed regime', 
                      labels = c('Perennial', 'Intermittent')) + 
    scale_color_manual(values=c('#175885', '#9e3f00'),
                       name='Observed regime', 
                       labels = c('Perennial', 'Intermittent')) + 
    facet_wrap(~variable, scales='free', labeller=label_value) + 
    theme_bw() +
    theme(legend.position = c(0.8, 0.1))
  
  return(list(uncertainty_numplot, uncertainty_catplot))
}

write_preds <- function(in_filestructure, in_gaugep, in_gaugestats, in_rftuned, predvars) {
  # ---- Output GRDC predictions as points ----
  in_gaugestats[!is.na(cly_pc_cav),  #SHould change that to not have to adjust subselection in multiple spots
                IRpredprob := in_rftuned$rf_inner$predict(in_rftuned$task)$prob[,2]]
  
  cols_toditch<- colnames(in_gaugep)[colnames(in_gaugep) != 'GRDC_NO']
  

  out_gaugep <- merge(in_gaugep, 
                      in_gaugestats[, !cols_toditch, with=F], by='GRDC_NO')
  
  st_write(obj=out_gaugep,
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
  
  # --------- Return data for plotting ------------------------
  return(out_gaugep)
}

