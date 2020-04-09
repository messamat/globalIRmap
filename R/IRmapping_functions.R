#Functions for intermittent river analysis

##### -------------------- Internal functions ---------------------------------
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
  paste('Importing', file_name, 'with ', length(keptcols), 'columns out of ', length(colsToKeep), 'supplied column names')
  fread(input=file_name, header=TRUE, select=keptcols, verbose=TRUE)
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

read_gaugep <- function(filestructure, dist) {
  #Import gauge stations and only keep those < dist m from a HydroSHEDS reach
  return(st_read(dsn=dirname(filestructure['in_gaugep']),
                 layer=basename(filestructure['in_gaugep'])) %>%
           .[.$station_river_distance<dist,]
  )
}
  
read_gauged_paths <- function(filestructure, in_gaugep) {
  #Get data paths of daily records for gauge stations
  fileNames <- file.path(filestructure['in_gaugedir'], paste(in_gaugep$GRDC_NO,  ".txt", sep=""))
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



