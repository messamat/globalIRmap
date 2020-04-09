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
          ggplot(edarf_outdt,  aes_string(x=xvar, y=yvar)) +
            geom_tile(aes(fill = `1`)) + 
            scale_fill_distiller(palette='Spectral') + 
            geom_jitter(data=pointdt, aes(color=intermittent)) +
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
  paste('Importing', file_name, 'with ', length(keptcols), 
        'columns out of ', length(colsToKeep), 'supplied column names')
  fread(input=file_name, header=TRUE, select=keptcols, verbose=TRUE)
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

format_gaugestats <- function(in_gaugestats, in_gaugep) {
  #Format gaugestats into data.table
  
  
  #Join intermittency statistics to predictor variables and subset to only include those gauges with at least
  gaugestats_join <- do.call(rbind, in_gaugestats) %>%
    setDT %>%
    .[as.data.table(in_gaugep), on='GRDC_NO'] %>%
    .[!is.na(totalYears_kept) & totalYears_kept>10,] # Only keep stations with at least 10 years of data
  
  #---- Compute derived predictor variables ----
  pre_mcols <- paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)) #Monthly precipitation columns
  gaugestats_join[, `:=`(pre_mm_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
                         pre_mm_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
                  .SDcols= pre_mcols] %>% 
    .[, `:=`(pre_mm_cvar=pre_mm_cmn/pre_mm_cmx, #min/max monthly catchment precip
             dis_mm_pvar=ifelse(dis_m3_pmx==0, 
                                1, 
                                dis_m3_pmn/dis_m3_pmx), #min/max monthly watershed discharge
             dis_mm_pvaryr=ifelse(dis_m3_pyr==0, 
                                  1, 
                                  dis_m3_pmn/dis_m3_pyr),#min monthly/average yearly watershed discharge
             ele_pc_rel = ifelse(ele_mt_uav==0, 
                                 0, 
                                 (ele_mt_cav-ele_mt_uav)/ele_mt_uav))] #catchment average elv - watershec average elev
  
  gaugestats_join[, cmi_ix_cmn := do.call(pmin, c(.SD)),
                  .SDcols= paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly cmi
    .[, swc_pc_cmn := do.call(pmin, c(.SD)),
      .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] #Get minimum monthly swc
  

  #---- Convert -9999 values to NA ----
  for (j in which(sapply(gaugestats_join,is.numeric))) { #Iterate through numeric column indices
    set(gaugestats_join,which(gaugestats_join[[j]]==-9999),j, NA)} #Set those to 0 if -9999
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

tune_rf <- function(in_gaugestats, in_predvars) {
  #---- Tune model ----
  predcols <- in_predvars$varcode
  # rf_formula <- as.formula(paste0('intermittent~', 
  #                                 paste(predcols, collapse="+"), 
  #                                 collapse=""))
  
  task_inter <- mlr3::TaskClassif$new(
    id ='inter_basic',
    backend = in_gaugestats[!is.na(cly_pc_cav), 
                            c('intermittent', predcols),with=F],
    target = "intermittent")
  
  #Set default ranger learner with explicit parameter set
  print(mlr3::lrn('classif.ranger')$param_set)
  
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
  
  # mlr_measures
  # measure = lapply(c("classif.sensitivity", 'classif.specificity'), msr)
  # basic800pred$score(measure)
  # 
  # 
  # tune_ps = ParamSet$new(list(
  #   ParamDbl$new("cp", lower = 0.001, upper = 0.1),
  #   ParamInt$new("minsplit", lower = 1, upper = 10)
  # ))
  # tune_ps
}







