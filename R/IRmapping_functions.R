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
#' diny(1999)
#' diny(2000)
#' diny(2004)
#' diny(2100)
#' diny(1600)
#'
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
#' test1 <- c(1,1,1,0,0,0,0,1,1)
#' zero_lomf(test)
#' test2 <- c(0,0,0,0,0,1,1,0,1)
#' zero_lomf(test2, first=FALSE)
#' zero_lomf(test2, first=TRUE)
#'
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

#------ readformatGRDC -----------------
readformatGRDC<- function(path) {
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]

  gaugetab <- cbind(fread(path, header=T, skip = 40, sep=";",
                          colClasses=c('character', 'character', 'numeric',
                                       'numeric', 'integer')),
                    GRDC_NO = gaugeno)%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)

  gaugetab[, year:=as.numeric(substr(dates, 1, 4))] %>% #Create year column and total number of years on record
    .[, prevflowdate := gaugetab[zero_lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA] #If non-zero flow, set prevlowfate to NA

  #Compute number of missing days per year
  gaugetab[!(Original %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Original)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']

  gaugetab[,`:=`(month = as.numeric(substr(dates, 6, 7)), #create a month column
                 integervalue = fifelse(Original == round(Original), 1, 0) #Flag integer discharge values
  )]

  return(gaugetab)
}

#------ readformatGSIMmon -----------------
readformatGSIMmon <- function(path) {
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  gaugetab <- cbind(fread(path, header=T, skip = 21, sep=",",
                          colClasses=c('Date', rep('numeric', 8),
                                       'integer', 'integer')),
                    gsim_no = gaugeno) %>%
    setnames(new=gsub('[\\\t]|["]', '', names(.))) %>%
    setorder(date) %>%
    .[, ':='(year = as.numeric(substr(date, 1, 4)),
             month = as.numeric(substr(date, 6, 7)))] %>%
    .[data.table(month=c(12, 1:11),
                 season=rep(c('DJF', 'MAM', 'JJA', 'SON'), each=3)),
      on='month']

  #Compute minimum number of zero-flow days for each month
  gaugetab[, mDur_minmo := fcase(MAX==0L, n.available,
                                 MIN7==0L, 7L,
                                 MIN==0L, 1L,
                                 default = 0L)]

  #Compute number of missing days per year
  gaugetab[, `:=`(missingdays = diny(year) - sum(n.available),
                  datadays = sum(n.available)),
           by= 'year']

  return(gaugetab)
}

#------ readformatGSIMsea -----------------
readformatGSIMsea <- function(path) {
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  gaugetab <- cbind(fread(path, header=T, skip = 21, sep=",",
                          colClasses=c('Date', rep('numeric', 17),
                                       'integer', 'integer')),
                    gsim_no = gaugeno) %>%
    setnames(new=gsub('[\\\t]|["]', '', names(.))) %>%
    setorder(date) %>%
    .[, ':='(year = as.numeric(substr(date, 1, 4)),
             month = as.numeric(substr(date, 6, 7)))] %>%
    .[data.table(month=c(2, 5, 8, 11),
                 season=c('DJF', 'MAM', 'JJA', 'SON')),
      on='month']

  #Compute minimum number of zero-flow days for each season
  gaugetab[, mDur_minsea := fcase(
    P90==0, as.integer(floor(0.9*n.available)),
    P80==0, as.integer(floor(0.8*n.available)),
    P70==0, as.integer(floor(0.7*n.available)),
    P60==0, as.integer(floor(0.6*n.available)),
    P50==0, as.integer(floor(0.5*n.available)),
    P40==0, as.integer(floor(0.4*n.available)),
    P30==0, as.integer(floor(0.3*n.available)),
    P20==0, as.integer(floor(0.2*n.available)),
    P10==0, as.integer(floor(0.1*n.available)),
    MIN7==0L, 7L,
    MIN==0L, 1L,
    default = 0L
  )]
  return(gaugetab)
}
#------ flagGRDCoutliers ------
flagGRDCoutliers <- function(in_gaugetab) {
  in_gaugetab %>%
    .[Original %in% c(-999, -99, -9999, 999, 9999), Original := NA] %>%
    .[Calculated %in% c(-999, -99, -9999, 999, 9999), Calculated := NA] %>%
    .[, `:=`(jday = format(as.Date(dates), '%j'),#Julian day
             q_rleid = rleid(Original),#Identify each group of consecutive values
             flag_mathis = 0)] #Create flag field)

  #Flag negative values
  in_gaugetab[Original < 0, flag_mathis := flag_mathis + 1]

  #Flag when more than 10 identical values in a row or when a single zero-flow
  in_gaugetab[, flag_mathis := flag_mathis +
                ((Original > 0) & (.N > 10)) +
                ((Original == 0) & (.N == 1)),
              by=q_rleid]

  #Flag |log(Q + 0.01) - mean| > 6SD for julian day mean and SD of 5d mean of log(Q + 0.01)
  in_gaugetab[, logmean5d := frollapply(log(Original + 0.01), n = 5, align='center',
                                        FUN=mean, na.rm = T)] %>% #Compute 5-day mean of log(Q+0.01)
    .[, `:=`(jdaymean = mean(logmean5d, na.rm = T),
             jdaysd = sd(logmean5d, na.rm = T)),
      by = jday] %>% #Compute mean and SD of 5-day mean of log(Q + 0.01) by Julian day
    .[abs(log(Original + 0.01) - jdaymean) > (6 * jdaysd),
      flag_mathis := flag_mathis + 1]

  #Flag values where Original != Calculated discharge
  in_gaugetab[Original != Calculated, flag_mathis := flag_mathis + 1]

  return(in_gaugetab)
}

#------ plotGRDCtimeseries ----------------------
plotGRDCtimeseries <- function(GRDCgaugestats_record,
                               outpath=NULL, maxgap = 366,  showmissing = FALSE) {
  #Read and format discharge records
  if (GRDCgaugestats_record[,.N>1] &
      ('Original' %in% names(GRDCgaugestats_record))) {
    gaugetab <- GRDCgaugestats_record
  } else {
    gaugetab <- readformatGRDC(GRDCgaugestats_record$path) %>%
      flagGRDCoutliers %>%
      .[, dates := as.Date(dates)] %>%
      .[!is.na(Original), missingdays := diny(year)-.N, by= 'year']
  }

  #Plot time series
  qtiles <- union(gaugetab[, min(Original, na.rm=T)],
                  gaugetab[, quantile(Original, probs=seq(0, 1, 0.1), na.rm=T)])

  subgaugetab <- gaugetab[missingdays < maxgap,]

  rawplot <- ggplot(subgaugetab[missingdays < maxgap,],
                    aes(x=dates, y=Original)) +
    geom_line(color='#045a8d', size=1, alpha=1/5) +
    geom_point(data = subgaugetab[flag_mathis == 0 & Original > 0,],
               color='#045a8d', size=1, alpha=1/3) +
    geom_point(data = subgaugetab[flag_mathis > 0 & Original > 0,],
               color='green') +
    geom_point(data = subgaugetab[flag_mathis == 0 & Original == 0,],
               color='red') +
    geom_point(data = subgaugetab[flag_mathis > 0 & Original == 0,],
               color='black') +
    scale_y_sqrt(breaks=qtiles, labels=qtiles) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(y='Discharge (m3/s)',
         title=paste0('GRDC: ', GRDCgaugestats_record$GRDC_NO)) +
    coord_cartesian(expand=0, clip='off')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text())

  if (showmissing) {
    rawplot <- rawplot +
      geom_point(data=gaugetab[missingdays >= maxgap,], color='black', alpha=1/10)
  }

  if (!is.null(outpath)) {
    if (!(file.exists(outpath))) {
      ggsave(filename = outpath, plot = rawplot, device = 'png',
             width = 10, height = 10, units='in', dpi = 300)
    }
  } else {
    return(rawplot)
  }
}

#------ plotGSIMtimeseries ----------------------
plotGSIMtimeseries <- function(GSIMgaugestats_record, outpath=NULL, maxgap=366,
                               showmissing = FALSE) {
  #Read and format discharge records
  gaugetab <- readformatGSIMmon(GSIMgaugestats_record$path) %>%
    .[!is.na(MEAN),] %>%
    setorder(date)

  #Format for plotting: compute MEAN - 2*SD if > MIN, otherwise MIN
  gaugetab[, `:=`(ribbonlow = max(c(MEAN-2*SD, MIN)),
                  ribbonhigh = min(c(MEAN+2*SD, MAX))
  ), by=date]

  #Plot time series
  qtiles <- gaugetab[, quantile(union(MIN, MAX), probs=seq(0, 1, 0.1))]

  rawplot <- ggplot(gaugetab[missingdays < maxgap,], aes(x=date, y=MEAN)) +
    geom_line(color='#045a8d', size=1, alpha=1/3) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MEAN>0),],
               color='#045a8d', size=1, alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MEAN==0),],
               color='red', size=1, alpha=1/2) +
    geom_ribbon(aes(ymin = ribbonlow, ymax = ribbonhigh), color='lightblue', alpha=1/4) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MIN>0),],
               aes(y = MIN), color='black', alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MIN==0),],
               aes(y = MIN), color='darkred', alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MAX>0),],
               aes(y = MAX), color='black', alpha=1/2) +
    geom_point(data=gaugetab[(missingdays < maxgap) & (MAX==0),],
               aes(y = MAX), color='purple', alpha=1/2) +
    scale_y_sqrt(breaks=qtiles, labels=qtiles) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(y='Discharge (m3/s)',
         title=paste0('GSIM: ', GSIMgaugestats_record$gsim_no)) +
    coord_cartesian(expand=0, clip='off')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text())

  if (showmissing) {
    gaugetab_removed <- gaugetab[missingdays > maxgap,]
    rawplot <- rawplot +
      geom_point(data=gaugetab[MEAN>0,],
                 color='#045a8d', size=1, alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MEAN==0),],
                 color='red', size=1, alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MIN>0),],
                 aes(y = MIN), color='black', alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MIN==0),],
                 aes(y = MIN), color='darkred', alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MAX>0),],
                 aes(y = MAX), color='black', alpha=1/5) +
      geom_point(data=gaugetab[(missingdays < maxgap) & (MAX==0),],
                 aes(y = MAX), color='purple', alpha=1/5)
  }

  if (!is.null(outpath)) {
    if (!(file.exists(outpath))) {
      ggsave(filename = outpath, plot = rawplot, device = 'png',
             width = 10, height = 10, units='in', dpi = 300)
    }
  } else {
    return(rawplot)
  }
}
#------ checkGRDCzeroes --------
#Function to plot and subset two weeks on each side of every zero-flow period in record
checkGRDCzeroes <- function(GRDCstatsdt, in_GRDC_NO, period=15, yearthresh,
                            maxgap, in_scales='free_x', labelvals) {
  check <- readformatGRDC(GRDCstatsdt[GRDC_NO ==  in_GRDC_NO, path]) %>%
    flagGRDCoutliers %>%
    .[, dates := as.Date(dates)] %>%
    .[!is.na(Original), missingdays := diny(year)-.N, by= 'year'] %>%
    .[missingdays < maxgap,] %>%
    merge(check[, list(dates=seq(min(dates), max(dates), by='day'))],
          by='dates', all.x=T, all.y=T)

  checkdates <- unique(do.call("c",
                               lapply(
                                 as.Date(check[Original==0, dates]),
                                 function(i) seq.Date(i-period, i+period, by='day'))))

  zeroes <- check[dates %in% checkdates,]
  zeroes[, zerogrp := rleid(round(difftime(dates, lag(dates))))] %>%
    .[, grpN := .N, by=zerogrp]

  p <- plotGRDCtimeseries(zeroes[grpN > 1 & year > yearthresh,],
                          outpath=NULL, maxgap=maxgap, showmissing=T) +
    scale_y_sqrt() +
    scale_x_date(breaks='1 month') +
    facet_wrap(~zerogrp, scales=in_scales)

  if (labelvals == T) {
    p <- p + geom_text(data=zeroes[grpN > 1 & Original >0,],
                       aes(label=Original), vjust=1, alpha=1/2)
  }

  print(p)
  return(zeroes)
}

#------ ggalluvium_gaugecount -------
ggalluvium_gaugecount <- function(dtformat, alluvvar) {
  p <- ggplot(dtformat,
              aes(x = variable, stratum = value, alluvium = get(alluvvar),
                  y = 1, fill = value, label =count)) +
    scale_x_discrete(labels=c('Full record', 'post-1961', 'post-1971')) +
    scale_y_continuous(name='Number of gauges') +
    scale_fill_discrete(labels=c('Perennial', 'Intermittent', '< 10 years of data')) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    coord_cartesian(expand=c(0,0), clip='off') +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank())
  return(p)
}

#------ plot_winterir -------
plot_winterir <- function(dt, dbname, inp_resdir, yearthresh, plotseries = TRUE) {
  #Get data subset
  wintergauges <- dt[get(paste0('winteronlyir_o', yearthresh)) == 1 &
                       get(paste0('totalYears_kept_o', yearthresh)) >= 10,]

  #Create output directory
  resdir_winterirplots <- file.path(
    inp_resdir,
    paste0(dbname, 'winterir_rawplots_o', yearthresh, '_',
           format(Sys.Date(), '%Y%m%d')
    )
  )

  if (!(dir.exists(resdir_winterirplots))) {
    print(paste0('Creating ', resdir_winterirplots ))
    dir.create(resdir_winterirplots )
  }

  #Generate plots depending on database
  if (plotseries) {
    if (str_to_upper(dbname) == 'GRDC') {
      lapply(wintergauges$GRDC_NO, function(gauge) {
        print(gauge)
        plotGRDCtimeseries(GRDCgaugestats_record =  wintergauges[GRDC_NO == gauge,],
                           outpath = file.path(resdir_winterirplots,
                                               paste0('GRDC', gauge, '.png')))
      })
    } else if (str_to_upper(dbname) == 'GSIM') {
      lapply(wintergauges$gsim_no, function(gauge) {
        print(gauge)
        plotGSIMtimeseries(GSIMgaugestats_record = wintergauges[gsim_no == gauge,],
                           outpath = file.path(resdir_winterirplots,
                                               paste0('GSIM', gauge, '.png')))
      })
    }
  }

  return(wintergauges)
}

#------ plot_coastalir -------
plot_coastalir <- function(in_gaugep = in_gaugep, dt = GRDCstatsdt, yearthresh,
                           dbname = 'grdc', inp_resdir = inp_resdir,
                           plotseries = TRUE) {

  idno <- fifelse(str_to_upper(dbname) == 'GRDC', 'GRDC_NO', 'gsim_no')

  coastalgauges <- in_gaugep[!is.na(in_gaugep$class99_19_rsp9_buf3k1) &
                               !is.na(as.data.frame(in_gaugep)[,idno]),] %>%setDT

  coastalir <- dt[get(paste0('winteronlyir_o', yearthresh)) == 1 &
                    get(paste0('totalYears_kept_o', yearthresh)) >= 10 &
                    get(idno) %in% coastalgauges[,get(idno)],]

  resdir_coastalirplots <- file.path(inp_resdir,
                                     paste0(str_to_upper(dbname),
                                            'coastalir_rawplots_o',
                                            yearthresh, '_',
                                            format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_coastalirplots))) {
    print(paste0('Creating ', resdir_coastalirplots ))
    dir.create(resdir_coastalirplots )
  }

  #Generate plots depending on database
  if (plotseries) {
    if (str_to_upper(dbname) == 'GRDC') {
      lapply(coastalir$GRDC_NO, function(gauge) {
        print(gauge)
        plotGRDCtimeseries(GRDCgaugestats_record = coastalir[GRDC_NO == gauge,],
                           outpath = file.path(resdir_coastalirplots, paste0('GRDC', gauge, '.png')))
      })
    } else if (str_to_upper(dbname) == 'GSIM') {
      lapply(coastalir$gsim_no, function(gauge) {
        print(gauge)
        plotGSIMtimeseries(GSIMgaugestats_record = coastalir[gsim_no == gauge,],
                           outpath = file.path(resdir_coastalirplots, paste0('GSIM', gauge, '.png')))
      })
    }
  }


  return(coastalir)
}

#------ comp_ymean -------------
#' Compute annual mean
#'
#' Compute annual mean of a variable based on columns of monthly means.
#'
#' @param in_dt data table, contains columns of monthly values.
#' @param fieldex character, example column name containing monthly value.
#' It must include the month number in %m format (two digits e.g. 02 for February)
#' @param mstart integer; position of the first digit of the month in the column
#' of monthly values.
#' @param outcol character; name of new column to be created containing monthly mean
#'
#' @return Modifies in_dt in place. Add column of weighted mean across months
#' (based on number of days in months)
#'
#' @export

comp_ymean<- function(in_dt, fieldex, mstart, outcol) {
  refd <- data.table(month=c('01', '02', '03', '04', '05', '06',
                             '07', '08', '09', '10', '11', '12'),
                     dimonth = c(31, 28.25, 31, 30, 31, 30,
                                 31, 31, 30, 31, 30, 31))
  refd[, mcols := paste0(substr(fieldex, 1, mstart-1),
                         month)]
  in_dt2 <- copy(in_dt[, refd$mcols, with=F])
  in_dt2[, comp_ymeanID := .I]
  in_dt[, outcol] <- melt(in_dt2,
                          id.vars='comp_ymeanID',
                          variable.name = 'mcols') %>%
    .[refd, on='mcols'] %>%
    .[, weighted.mean(value, dimonth), by=comp_ymeanID] %>%
    .[,'V1',with=F]
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
#'   relevant for the intermittency analysis. \cr
#'   \cr
#'   Steps include: \cr
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


  #---- Inspect and correct -9999 and NA values ----
  print('Inspect and correct -9999 values')
  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  in_dt2[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  in_dt2[snw_pc_cmx == -9999, snw_pc_cmx:=0]

  #Places with very high slope value (in Greenland) have -9999 - replace with high value
  #in_dt2[, max(slp_dg_uav)]
  in_dt2[slp_dg_cav == -9999, `:=`(slp_dg_cav = 750,
                                   slp_dg_uav = 650)]

  #bio5 is -9999 in a few places in Greenland. Set at lowest existing value
  # in_dt2[bio5_dc_cav>-9999, min(bio5_dc_cav)]
  # in_dt2[bio5_dc_uav != -9999, min(bio5_dc_uav)]
  in_dt2[bio5_dc_cav == -9999, bio5_dc_cav := -950]
  in_dt2[bio5_dc_uav == -9999, bio5_dc_uav:= -9500]

  print('Number of NA values per column')
  colNAs<- in_dt2[, lapply(.SD, function(x) sum(is.na(x)))]
  print(colNAs)

  print('Number of -9999 values per column')
  col9999<- in_dt2[, lapply(.SD, function(x) sum(x==-9999))]
  print(col9999)

  #-9999 in cly_pc_cav, slt, and snd are places with no soil mask (urban areas, lakes, glaciers, etc.)


  #Define column groups
  gladcols <- unlist(lapply(c('cav', 'uav'),function(s) {
    lapply(c('wloss', 'wdryp', 'wwetp', 'whfrq', 'wseas', 'wperm', 'wfresh',
             'wwpix', 'wdpix'),
           function(v) {
             paste0(v, '_pc_', s)
           })
  })
  )

  sgcols <- c('cly_pc_cav','slt_pc_cav', 'snd_pc_cav',
              'cly_pc_uav','slt_pc_uav', 'snd_pc_uav')

  #Bioclim columns in celsius degrees
  biocolsdc <- unlist(lapply(c('cav', 'uav'),
                             function(s) paste0('bio', 1:11, '_dc_', s)
  ))
  biocolsnegative <- grep('bio[189][01]*_dc_uav', biocolsdc, value=T)

  cmi_cmcols <- paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))
  cmi_umcols <- paste0('cmi_ix_u', str_pad(1:12, width=2, side='left', pad=0))

  #Convert -9999 to NAs
  for (j in which(sapply(in_dt2,is.numeric))) { #Iterate through numeric column indices
    if (!(j %in% which(names(in_dt2) %in% c(gladcols, sgcols, 'wet_cl_cmj')))) {
      set(in_dt2,which(in_dt2[[j]]==-9999),j, NA) #Set those to NA if -9999
    }
  }

  #Scale variables based on HydroATLAS v1.0 documentation and v1.0.9 processing
  in_dt2[, `:=`(
    ari_ix_cav = ari_ix_cav/1000,
    ari_ix_uav = ari_ix_uav/1000,
    dor_pc_pva = dor_pc_pva/100,
    lka_pc_cse = lka_pc_cse/10,
    lka_pc_use = lka_pc_use/10,
    gwt_m_cav = gwt_cm_cav/100
  )]

  in_dt2[, (gladcols) := lapply(.SD, function(x) x/100), .SDcols = gladcols]
  in_dt2[, (biocolsdc) := lapply(.SD, function(x) (x/100)), .SDcols = biocolsdc]
  in_dt2[, (biocolsnegative) := lapply(.SD, function(x) x-100),
         .SDcols=biocolsnegative]

  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  comp_ymean(in_dt=in_dt2, fieldex = 'cmi_ix_c01', mstart=9, outcol='cmi_ix_cyr')
  comp_ymean(in_dt=in_dt2, fieldex = 'cmi_ix_u01', mstart=9, outcol='cmi_ix_uyr')
  comp_ymean(in_dt=in_dt2, fieldex = 'pet_mm_c01', mstart=9, outcol='pet_mm_cyr')
  #comp_ymean(in_dt=in_dt2, fieldex = 'pet_mm_u01', mstart=9, outcol='pet_mm_uyr')

  in_dt2[,
         `:=`(cmi_ix_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
              cmi_ix_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
         .SDcols= cmi_cmcols] %>%
    .[, `:=`(cmi_ix_umn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
             cmi_ix_umx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
      .SDcols= cmi_umcols] %>%
    .[, swc_pc_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
      .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly swc
    .[, `:=`(#min/max monthly CMI
      cmi_ix_cvar= fifelse(cmi_ix_cmx==0, 0, cmi_ix_cmn/cmi_ix_cmx),
      cmi_ix_uvar= fifelse(cmi_ix_umx==0, 0, cmi_ix_umn/cmi_ix_umx),
      #min/max monthly watershed discharge
      dis_m3_pvar=fifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx),
      #min monthly/average yearly watershed discharge
      dis_m3_pvaryr=fifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),
      #catchment average elv - watershec average elev
      ele_pc_rel = fifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav),
      #runoff coefficient (runoff/precipitation)
      runc_ix_cyr = fifelse(bio12_mm_cav==0, 0, run_mm_cyr/bio12_mm_cav),
      #Specific discharge
      sdis_ms_uyr = dis_m3_pyr/UPLAND_SKM,
      sdis_ms_umn = dis_m3_pmn/UPLAND_SKM
    )]
  return(in_dt2)
}

#------ threshold_misclass -----------------
#' Threshold misclassification
#'
#' Compute the misclassification rate, sensitivity, and specificity of probability
#' predictions to a binary classification problem from a random forest learner
#' based on a given threshold.
#'
#' @param i (numeric) Threshold to be used for converting probabilities
#'   predictions to binary predictions, limits=\[0,1\].
#' @param in_preds Either:
#' 1. a \link[mlr3]{PredictionClassif} or
#' 2. a data.table of predictions for a set of CV repetitions as formatted by
#' \code{\link{analyze_benchmark}}. The column of the predicted probability of
#' getting a 'positive' must be called \code{prob.1}.
#'
#' @return A single-row data.table with four columns:
#' \itemize{
#'   \item i - the threshold used to compute the classification statistics
#'   \item misclass - misclassification rate \[0-1\]
#'   \item sens - sensitivity \[0-1\]
#'   \item spec - specificity \[0-1\]
#' }
#'
#' @details
#' \itemize{
#' \item Misclassification rate is the proportion of misclassified
#' records
#' \item Sensitivity is the proportion of actual positives (1) that are
#' correctly identified as such (the proportion of intermittent rivers that
#' are identified as intermittent)
#' \item Specificity is the proportion of actual negatives that are correctly
#' identified as such (the proportion of perennial rivers that are identified as
#' perennial)
#' } \cr
#' See \url{https://en.wikipedia.org/wiki/Sensitivity_and_specificity}
#'
#' @section Warning:
#' This function was only tested for the outputs of a probability random forest
#' (using a classification framework) and a regression forest (using 0 or 1 as
#' dummy continuous variables)
#'
#' @export

threshold_misclass <- function(i=0.5, in_preds) {
  #---- Get confusion matrix ----
  if (inherits(in_preds, 'PredictionClassif')) {
    confu <- as.data.table(in_preds$set_threshold(1-i)$confusion) #Get confusion matrix directly
  }

  if (is.data.table(in_preds)) {
    #If task associated with predictions is a classification and has records
    if (in_preds[task_type == 'classif',.N] > 0) {
      #For each CV repetition:
      #   1. set the probability threshold to compute a confusion matrix to 1-i
      #     (i being the threshold to classify something as 1, set_threshold
      #     being based on prob.0, not prob.1)
      #   2. Compute confusion matrix
      confu <- in_preds[, as.data.table(pred[[1]]$set_threshold(1-i)$confusion),
                        by=outf] %>%
        #Aggregate confusion matrices across repetitions
        .[, .(N=sum(N)), by=.(response, truth)]
    }

    #If task associated with predictions is a regression and has records
    if (in_preds[task_type == 'regr', .N] > 0) {
      #Reclassify continuous predictions into binary response across all records
      confu <- in_preds[, response := fifelse(prob.1>=i, '1', '0')] %>%
        .[, truth := as.character(truth)] %>%
        #Create aggregate confusion matrix
        .[, .N, by=.(response, truth)]
    }
  }

  #---- Compute statistics based on confusion matrix and format into data.table----
  outvec <- data.table(
    i,
    misclas = confu[truth != response, sum(N)] / confu[, sum(N)],
    sens = confu[truth == '1' & response == '1', N] / confu[truth=='1', sum(N)],
    spec  = confu[truth=='0' & response==0, N]/confu[truth=='0', sum(N)]
  )
  return(outvec)
}

#------ threshold_dat ----------------
threshold_dat <- function(bmres) {
  bmres_dt <- as.data.table(bmres) %>%
    .[, learner_id := unlist(lapply(learner, function(x) x$id))] %>%
    .[, task_id := unlist(lapply(task, function(x) x$id))] %>%
    .[, resampling_id := unlist(lapply(resampling, function(x) x$id))] %>%
    unique(by='uhash')

  if (bmres$task_type == 'regr') {
    preds <- lapply(seq_len(bmres_dt[,.N]), function(rsmp_i) {
      preds <- bmres_dt$prediction[[rsmp_i]] %>%
        as.data.table %>%
        .[, `:=`(outf = bmres_dt$iteration[[rsmp_i]],
                 task = bmres_dt$task[[rsmp_i]]$id,
                 task_type = bmres$task_type,
                 learner = bmres_dt$learner[[rsmp_i]]$id)]
      return(preds)
    }) %>%
      do.call(rbind, .)

    if (!('prob.1' %in% names(preds)) & 'response' %in% names(preds)) {
      preds[, prob.1 := response]
    }
  }

  outtab <- lapply(1:bmres$n_resample_results, function(rsmp_i) {
    print(rsmp_i)
    rsmp_preds <- bmres$resample_result(rsmp_i)$prediction()

    if (bmres$task_type == 'regr') {
      rsmp_preds <- preds
    }

    baccthresh <- ldply(seq(0,1,0.01), function(threshold_class) {
      print(threshold_class)
      cbind(
        rsmp_bacc(bmres, rsmp_i=rsmp_i, threshold_class=threshold_class),
        threshold_misclass(i=threshold_class, in_preds=rsmp_preds)
      ) %>%
        .[, i:=NULL]
    }) %>%
      setDT %>%
      cbind(bmres_dt[rsmp_i, .(learner_id, task_id, resampling_id)])
    return(baccthresh)
  }) %>%
    rbindlist


  outtab_melt <- melt(outtab, id.vars=c('threshold_class',
                                        'learner_id', 'task_id', 'resampling_id'))
  return(outtab_melt)
}
#------ get_outerrsmp -----------------
#' Get outer ResampleResult
#'
#' Extract the outer resampling \link[mlr3]{ResampleResult} from the output
#' from \code{\link{selecttrain_rf}} with the option to grab the output from
#' the spatial resampling strategy or not (when there are multiple resampling
#' strategies).
#'
#' @param in_rftuned (list or \link[mlr3]{ResampleResult}) Output from
#' \code{\link{selecttrain_rf}}.
#' @param spatial_rsp (logical) Whether to return the results from the
#' spatial resampling strategy if there is one.
#'
#' @details
#' If \code{in_rftuned} is already a \link[mlr3]{ResampleResult}, then simply
#' return it.
#'
#' @return a \link[mlr3]{ResampleResult}.
#'
#' @section Warning:
#' This will only return the \link[mlr3]{ResampleResult} from a single resampling
#' strategy.
#'
#' @export

get_outerrsmp <- function(in_rftuned, spatial_rsp=FALSE) {
  #Adapt whether return resample result or output from selecttrain_rf
  if (inherits(in_rftuned, 'list')) {
    #If there is more than one type of outer resampling
    if (length(in_rftuned$rf_outer$uhash) > 1) {
      #Check which resampling is spatial — only works if one is spatial
      sp_i <- which(unlist(lapply(in_rftuned$rf_outer, function(x) {
        grepl('.*Resampling.*Sp.*', x$resampling$format())
      })))

      #If user request that spatial resampling be used
      if (spatial_rsp==TRUE) {

        if (length(sp_i)>0) {
          rsmp_res <- in_rftuned$rf_outer[[min(sp_i)]]
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
      print("Only one resampling result, ignoring spatial_rsp argument...")
      rsmp_res <- in_rftuned$rf_outer
    }
    #If in_rftuned is alreadyh a ResampleResult, simply return it
  } else if (inherits(in_rftuned, "ResampleResult")) {
    rsmp_res <- in_rftuned
  }
  return(rsmp_res)
}

#------ weighted_sd -----------------
#' Weighted standard deviation
#'
#' Computes a weighted standard deviation
#'
#' @param x An object containing the values whose weighted mean is to be
#' computed.
#' @param w An integer vector of weights the same length as x giving the
#' weights to use for elements of x (weights are > 1, not decimals).
#' @param na.rm	(logical) Whether NA values in x should be stripped before the
#' computation proceeds.
#'
#' @return A length-one numeric vector.
#'
#' @examples
#' ```{r}
#' wt <- c(5, 4, 3 ,2 , 1)
#' x <- c(1, 2, 3, 4, 5)
#' xm <- weighted.mean(x, wt)
#' xsd <- weighted_sd(x, wt)
#' ```
#'
#' @export

weighted_sd <- function(x, w=NULL, na.rm=FALSE) {
  if (na.rm) {
    x <-  na.omit(x)
    if (length(w) > length(x)) {
      w <- w[-which(is.na(x))]
    }
  }

  if (length(w)==0) {
    w <- rep(1, length(x))
  }

  #Compute weighted standard deviation
  return(sqrt(sum((w) * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)))
}

#------ extract_impperf_nestedrf -----------------
#' Extract variable importance and performance from a trained RF learner
#'
#' Computes the performance, variable importance and associated p-value from
#' either a trained \link[mlr3]{AutoTuner} or a trained \link[mlr3]{Learner}.
#'
#' @param in_rflearner An \link[mlr3]{AutoTuner} or a trained
#' \link[mlr3]{Learner}.
#' @param imp (logical) Whether to compute variable importance.
#' @param perf (logical) Whether to compute the performance measure.
#' @param pvalue (logical) Whether to compute p-values.
#' @param pvalue_permutn (integer) number of permutations to use in p-value
#' calculations
#'
#' @return A data.table with, when all logicals are TRUE, the following columns.
#' \describe{
#'   \item{varnames} - predictor variable name
#'   \item{imp} - variable importance value
#'   \item{pvalue} - variable importance p-value
#'   \item{\[perf\]} - performance value (same for all rows), name changes
#' }
#'
#' @details
#' Accepts both a trained \link[mlr3]{AutoTuner} or a trained \
#' link[mlr3]{Learner} of class \link[mlr3learners]{mlr_learners_classif.ranger}
#' or \link[mlr3learners.partykit]{mlr_learners_classif.cforest}. \cr
#'
#' Also accept learners of class \link[mlr3pipelines]{GraphLearner} \cr
#'
#' If \code{p-value==TRUE} and \link[mlr3learners]{mlr_learners_classif.ranger},
#' then compute p-values with \link[ranger]{importance_pvalues} with Altmann
#' permutation method using \code{pvalue_permutn} permutations. \cr
#'
#' If \link[mlr3learners.partykit]{mlr_learners_classif.cforest} is provided,
#' \code{p-value} is ignored for now.
#'
#' @seealso \link{weighted_vimportance_nestedrf}
#' @section documentation to-do:
#' Can add an example down the line, add source.
#' @export

extract_impperf_nestedrf <- function(in_rflearner, in_task,
                                     imp = TRUE, perf = TRUE,
                                     pvalue = TRUE, pvalue_permutn = 100) {

  in_task <- in_task[[1]]
  in_rflearner <- in_rflearner[[1]]

  if (inherits(in_rflearner, "AutoTuner")) {
    sublrn <- in_rflearner$model$learner
  } else {
    sublrn <- in_rflearner
  }

  print(paste0("Computing variable importance for resampling instance hash #",
               sublrn$hash))

  outobj <- cbind(
    if (imp) {
      ####################### IF GraphLearner ####################################
      if (inherits(sublrn, "GraphLearner")) {

        if ('classif.ranger' %in% names(sublrn$model)) {

          if (pvalue == TRUE) {
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
            warning("p_value calculation is only available for ranger classification rf, ignoring p_value.
                  In addition, default parameters were used in partykit::varimp, adjust as needed.")
          }
          data.table(importance=
                       partykit::varimp(sublrn$model$classif.cforest$model,
                                        nperm = 1,
                                        OOB = TRUE,
                                        risk = "misclassification",
                                        conditional = FALSE,
                                        threshold = .2))
        }
      } else { ####################### IF direct model ####################################
        if (pvalue == TRUE) { #If want pvalue associated with predictor variables
          in_formula <- as.formula(paste0(in_task$target_names, '~.'))

          importance_pvalues(
            in_rflearner$model,
            method = "altmann",
            num.permutations = pvalue_permutn,
            data = in_task$data(),
            formula= in_formula
          )
        } else { #If pvalue == FALSE
          data.table(importance= in_rflearner$model$learner$importance())
        }

      }
    },
    if (perf) {
      perf_id <- in_rflearner$instance_args$measure$id
      outperf <- in_rflearner$tuning_result[, get(perf_id)]
      data.table(outperf) %>% setnames(perf_id)
    }
  )

  return(outobj)
}
#------ weighted_vimportance_nestedrf -----------------
#' Weighted mean of variable importance for resampled RF learner
#'
#' Compute mean and standard deviation of variable importance and mean of
#' p-value across resampling instances (e.g. folds and repetitions) weighted by
#' resampling prediction accuracy.
#'
#' @param rfresamp the \link[mlr3]{ResampleResult} from a classification RF
#' of type \link[mlr3learners]{mlr_learners_classif.ranger}
#' or \link[mlr3learners.partykit]{mlr_learners_classif.cforest}
#' @inheritParams extract_impperf_nestedrf
#'
#' @return A data.table with, when all logicals are TRUE, the following columns.
#' \describe{
#'   \item{varnames} - predictor variable name
#'   \item{imp_wmean} - weighted mean of variable importance across resampling
#'   instances
#'   \item{imp_wsd} - weighted standard deviation of variable importance across
#'   resampling instances
#'   \item{pvalue_wmean} -weighted mean of variable p-value across resampling
#'   instances
#' }
#'
#' @details
#' See \link{extract_impperf_nestedrf} for more details on computations.
#'
#' @seealso \link{extract_impperf_nestedrf} and \link{ggvimp} and
#' \link{benchmark_featsel}
#'
#' @section Warning:
#' Does not accept error measure for weighting
#' - could be amended to be based on 1-error
#'
#' @section documentation to-do:
#' Can add an example down the line, add source.
#' @export

weighted_vimportance_nestedrf <- function(rfresamp,
                                          pvalue = TRUE, pvalue_permutn = 100) {
  varnames <- rfresamp$task$feature_names
  rfresampdt <- as.data.table(rfresamp)

  vimportance_all <- rfresampdt[, extract_impperf_nestedrf(
    in_rflearner = learner, #Extract vimp and perf for each resampling instance
    in_task = task,
    imp=T, perf=T, pvalue=pvalue, pvalue_permutn), by=iteration] %>%
    cbind(., varnames)

  ####!!!!!!!!!!Adapt to allow for other measure than classif.bacc!!!!!!!!######
  out_vimportance <- vimportance_all[
    , list(imp_wmean = weighted.mean(importance, classif.bacc), #Compute weighted mean
           imp_wsd =  weighted_sd(importance, classif.bacc)), #Compute weighted sd
    by=varnames]

  if (pvalue) {
    out_vimportance <- cbind(
      out_vimportance,
      vimportance_all[,
                      list(imp_pvalue = weighted.mean(pvalue, classif.bacc)), #Compute weighted mean of pvalue
                      by=varnames][, !'varnames']
    )
  }

  return(out_vimportance)
}

#------ extract_pd_nestedrf -----------------
#' Extract partial dependence (and performance) from a trained RF learner.
#'
#' Computes the marginal relationship between a subset of the predictors
#' (here, two variables at a time) and the model’s predictions by averaging
#' over the marginal distribution of the compliment of this subset of
#' the predictors, taking in account the interaction between the chosen
#' predictors.
#'
#' @param learner_id (integer) Index of the outer resampling instance to be
#' analyzed.
#' @param in_rftuned \link[mlr3]{ResampleResult} from a classification RF.
#' @param datdf Data from the task that was used to train RF.
#' @param selcols Character vector of the predictor variables to analyze.
#' @param ngrid (integer) Number of values of the  predictor variables over
#' which to compute the marginal relationship.
#'
#' @return A data.table with the following columns.
#' \describe{
#'   \item{value1} - value of the first predictor variable in the pair
#'   \item{value2} - value of the second predictor variable in the pair
#'   \item{0} - predicted probability of 0 (e.g. probability that the river
#'   is perennial) at value1 and value2
#'   \item{1} - predicted probability of 1 (e.g. probability that the river
#'   is intermittent) at value1 and value2
#'   \item{\[perf\]} - performance value (same for all rows), name changes
#'   \item{var1} - name of the first predictor variable
#'   \item{var2} - name of the second predictor variable
#' }
#'
#'
#' @details
#' Also accept learners of class \link[mlr3pipelines]{GraphLearner}. \cr
#' Uses \link[edarf]{partial_dependence} for computing.
#'
#' @seealso \link{weighted_vimportance_nestedrf},
#'  \link{ggpd_bivariate}
#'
#' @section Warning:
#' Has only been tested on \link[mlr3learners]{mlr_learners_classif.ranger}
#'
#' @section documentation to-do:
#' Can add an example down the line, add source.
#'
#' @export

extract_pd_nestedrf <- function(learner_id=1, in_rftuned, datdf,
                                selcols, nvariate, ngrid) {
  in_mod <- as.data.table(in_rftuned)[eval(learner_id),] #Go through data.table format to have access to both tasks and learners

  #Get fold-specific performance measure
  foldperf <- extract_impperf_nestedrf(in_rflearner = in_mod$learner,
                                       in_task = in_mod$task,
                                       imp=F, perf=T, pvalue=F)

  # selcols <- in_vimp_plot$data %>% #Can use that if extracting from tunredrf is expensive
  #   setorder(-imp_wmean) %>%
  #   .[colnums, variable]


  if (inherits(in_mod$learner[[1]]$learner, "GraphLearner")) {
    in_fit <- in_mod$learner[[1]]$learner$model$classif.ranger$model
  } else {
    in_fit <- in_mod$learner[[1]]$learner$model
  }

  ngridvec <- c(ngrid, ngrid)

  #Make dataset of all combinations of selected column names, two at a time
  if (nvariate == 1) {
    pdcomb <- lapply(selcols, function(i) {
      print(i)
      pdout <- edarf::partial_dependence(fit = in_fit, vars = c(i),
                                         n = ngridvec, data = datdf) %>% #Warning: does not work with data_table
        setDT %>%
        .[,(names(foldperf)) := foldperf] %>%
        .[, `:=`(var1=i)] %>%
        setnames(i, 'value1')
    }
    ) %>%
      do.call(rbind, .)

  } else if (nvariate == 2) {
    vargrid <- combn(selcols, 2, simplify=F) %>%
      do.call(rbind, .)

    #Get marginal distribution of the effect of two columns at a time
    pdcomb <- mapply(function(i, j) {
      pdout <- edarf::partial_dependence(fit = in_fit, vars = c(i, j),
                                         n = ngridvec,
                                         interaction = TRUE, data = datdf) %>% #Warning: does not work with data_table
        setDT %>%
        .[,(names(foldperf)) := foldperf] %>%
        .[, `:=`(var1=i, var2=j)] %>%
        setnames(c(i,j), c('value1', 'value2'))


      return(pdout)
    }, vargrid[,1], vargrid[,2], SIMPLIFY = FALSE) %>%
      do.call(rbind, .)
  } else {
    print('Warning: function cannot yet work with more than two variables at a time')
  }

  return(pdcomb)
}

#------ fread_cols -----------------
#' fread columns
#'
#' Fast data.table-based reading of a subset of columns from a table
#'
#' @param file_name path of the table to be read
#' @param cols_tokeep character vector, names of columns to read
#'
#' @return a data.table with the \code{cols_tokeep} that are found in the table
#'
#' @examples
#' fread_cols(iris, c('Sepal.Length', 'Sepal.Width')
#'
#' @export

fread_cols <- function(file_name, cols_tokeep) {
  #Only read the first row from the file
  header <- fread(file_name, nrows = 1, header = FALSE)
  #Check which columns are in the table
  keptcols <- cols_tokeep[cols_tokeep %chin% unlist(header)]
  missingcols <- cols_tokeep[!(cols_tokeep %chin% unlist(header))]
  paste('Importing', file_name, 'with ', length(keptcols),
        'columns out of ', length(cols_tokeep), 'supplied column names')
  #Reading in table
  paste(missingcols, 'columns are not in the file...')

  dt <- fread(input=file_name, header=TRUE,
              select=keptcols, verbose=TRUE)
  return(dt)
}

#------ get_oversamp_ratio -----------------
#' Get oversample ratio
#'
#' Identify minority class and compute ratio between the number of observations
#' in the majority and the minority classes of a binary
#' \link[mlr3]{TaskClassif }.
#'
#' @param in_task binary \link[mlr3]{TaskClassif }
#'
#' @return named list with following items:
#' \describe{
#' \item{minoclass} - Value of minority class (e.g. '1' for intermittent rivers)
#' \item{ratio} - numeric ratio between number of items in majority and minority
#' class
#' }
#'
#' @examples
#' ```{r}
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "in_dt", backend = in_dt, target ='intermittent')
#' get_oversamp_ratio(task)
#' ```
#'
#' @export

#When given an mlr3 task with binary classification target, gets which class is minority and the ratio

get_oversamp_ratio <- function(in_task) {
  return(
    in_task$data()[, .N, by=get(in_task$target_names)] %>%
      setorder(N) %>%
      .[, list(minoclass=get[1], ratio=N[2]/N[1])]
  )
}

# get_oversamp_ratio <- function(in_task, classcol=NULL) {
#   if (inherits(in_task, 'Task')) {
#     if (is.null(classcol)) {
#       classcol <- in_task$target_names
#     }
#     dat <- in_task$data()[, .N, by = classcol] %>%
#       setnames(old=classcol, new='classcol')
#   } else if (inherits(in_task, 'data.table')) {
#     dat <- in_task[, .N, by=classcol]%>%
#       setnames(old=classcol, new='classcol')
#   }
#   return(
#     dat %>%
#       setorder(N) %>%
#       .[, list(minoclass=classcol[1], ratio=N[2]/N[1])]
#   )
# }

#------ convert_clastoregrtask -----------------
#' Convert classification task to regression task
#'
#' Convert a binary \link[mlr3]{TaskClassif} or
#' \link[mlr3spatiotempcv]{TaskClassifST}
#' to a \link[mlr3]{TaskRegr}, optionally oversampling the minority class first.
#'
#' @param in_task \link[mlr3]{TaskClassif} or \link[mlr3spatiotempcv]{TaskClassifST}
#' @param in_id the id of the output task
#' @param oversample (logical) whether to oversample minority class beforehand
#' (this is not exactly equivalent to a PipeOp) but oversample PipeOp is not yet
#' available for regression tasks. If TRUE, will oversample so that classes are
#' equally represented in task.
#'
#' @return a \link[mlr3]{TaskRegr}
#'
#' @examples
#' ```{r}
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "task_classif", backend = in_dt, target ='intermittent')
#' convert_clastoregrtask(task, id="task_regr")
#' ```
#'
#' @export


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
#------ durfreq_parallel -----------------
#' Parallel wrapper for \link{comp_durfreq}
#'
#' Wrapper to run \link{comp_durfreq} (which computes the annual mean number of zero
#' flow days and frequency of zero-flow periods for a streamflow time series)
#' in parallel across a list of time series (paths).
#'
#' @param pathlist vector or list of paths to time series formatted GRDC-style
#' @inheritParams comp_durfreq
#' @param monthsel_list  named list with names being the gauge ID (e.g. GRDC_NO)
#' and values the selected months to compute the statistics over.
#' @param reverse (logical) whether to use monthsel_list as excluding months
#' rather than subsetting
#'
#' @return data.table with each row a gauging station and columns inherited
#' from \link{comp_durfreq}
#'
#' @export
durfreq_parallel <- function(pathlist, maxgap, monthsel_list=NULL,
                             reverse=FALSE) {

  cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)

  return(as.data.table(rbindlist(
    foreach(j=pathlist,
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

#------ rsmp_bbrier ----------------
rsmp_bbrier <- function(bmres, rsmp_i) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  if (inherits(rsmp_pred, 'PredictionClassif')) {
    bbrier <- sum((as.numeric(as.character(rsmp_pred$truth)) -
                     rsmp_pred$prob[,'1'])^2)/length(rsmp_pred$row_ids)
  }

  if (inherits(rsmp_pred, 'PredictionRegr')) {
    bbrier <- sum((as.numeric(as.character(rsmp_pred$truth)) -
                     rsmp_pred$response)^2)/length(rsmp_pred$row_ids)
  }
  return(bbrier)
}

#------ rsmp_auc ----------------
rsmp_auc <- function(bmres, rsmp_i) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  auc <- Metrics::auc(rsmp_pred$truth, rsmp_pred$response)
  return(auc)
}

#------ rsmp_bacc ----------------
rsmp_bacc <- function(bmres, rsmp_i, threshold_class) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred <- rsmp$prediction()

  if (inherits(rsmp_pred, 'PredictionClassif')) {
    prob_resp <- rsmp_pred$prob[,'1']
  }

  if (inherits(rsmp_pred, 'PredictionRegr')) {
    prob_resp <- rsmp_pred$response
  }

  response <- as.data.table(prob_resp)[
    , c(fifelse(prob_resp>=threshold_class, '1', '0'))]

  bacc <- mlr3measures::bacc(factor(as.character(rsmp_pred$truth),
                                    levels=c('0', '1')),
                             factor(response, levels=c('0','1')))
  return(data.table(bacc=bacc, threshold_class=threshold_class))
}

#------ bm_paramstime ----------------
bm_paramstime <- function(bmres) {

  lrns_dt <- as.data.table(bmres) %>%
    unique(by='uhash')

  params_format <- lapply(lrns_dt$learner, function(lrn) {
    print(lrn$id)
    #Get general parameters-----------------------------------------------------
    if (inherits(lrn ,'AutoTuner')) {
      lrn_in <- lrn$learner
    } else {
      lrn_in <- lrn
    }

    #If GraphLearner, get parameters for each Pipe Operator
    if (inherits(lrn_in, 'GraphLearner')) {
      paramtab <- lapply(lrn_in$graph$pipeops, function(po) {
        #Get default params
        def <- as.data.table(po$param_set$default) %>%
          .[, po_id := po$id] %>%
          melt(id.var='po_id')

        #Get set param
        value <- as.data.table(po$param_set$values) %>%
          .[, po_id := po$id]%>%
          melt(id.var='po_id')

        #Merge default and set params, keeping set ones when present
        if (ncol(def)==1) {
          tabmerge <- value
        } else if (ncol(value)==1) {
          tabmerge <- def
        } else {
          tabmerge <- merge(def, value[, .(variable, value, po_id)],
                            by=c('po_id', 'variable'),
                            all.x=T, all.y=T) %>%
            .[is.na(value.y), value.y := value.x] %>%
            setnames('value.y', 'value') %>%
            .[, value.x := NULL]
        }

        graphparams <- tabmerge %>%
          .[, variable := gsub('[.]', '', variable)] %>%
          unique(by='variable')

        return(graphparams)
      }) %>%
        rbindlist(use.names=TRUE)

      #If direct learner
    } else {
      #Get default params
      def <- as.data.table(lrn_in$param_set$default) %>%
        melt(id.vars=integer())

      #Get set param
      value <- as.data.table(lrn_in$param_set$values) %>%
        melt(id.vars=integer())

      if (def[,.N]>0 & value[,.N]>0) {
        tabmerge <- merge(def, value, by='variable', all.x=T, all.y=T)%>%
          .[is.na(value.y), value.y := value.x] %>%
          setnames('value.y', 'value') %>%
          .[, value.x := NULL]
      } else {
        tabmerge <- rbind(def, value)
      }

      #Merge default and set params, keeping set ones when present
      paramtab <- tabmerge %>%
        .[, variable := gsub('[.]', '', variable)] %>%
        unique(by='variable') #Models are not consistent in their use of dots
    }

    #Get inner sampling parameters----------------------------------------------
    if ('instance_args' %in% names(lrn)) {
      param_ranges <- as.list(
        paste0(lrn$instance_args$search_space$lower,
               '-',
               lrn$instance_args$search_space$upper))
      names(param_ranges) <- lrn$instance_args$search_space$ids()


      inner_smp  <- cbind(
        as.data.table(
          lrn$instance_args$resampling$param_set$get_values()),
        as.data.table(
          lrn$instance_args$terminator$param_set$get_values())
      ) %>%
        .[, nruns := prod(.SD)] %>%
        setnames(names(.), paste0('inner_', names(.))) %>%
        cbind(as.data.table(param_ranges))

      #If no inner sampling
    } else {
      inner_smp <- data.table(inner_nruns = 1)
    }


    #Remove pipe operator names and dot from names (inconsistent use across learners)
    setnames(inner_smp,
             old=names(inner_smp),
             new=gsub('[.]', '',
                      gsub(paste(paste0('(', unique(paramtab$po_id), '[.]*)'),
                                 collapse = '|'),
                           '',
                           names(inner_smp))))

    #---------------------------------------------------------------------------
    #Merge inner sampling parameters and general parameters, then format
    return(
      paramtab %>%
        .[!(variable %in% names(inner_smp)),] %>%
        dcast(.~variable, value.var='value') %>%
        cbind(inner_smp)
    )
  }) %>%
    rbindlist(., fill=T)
  ##############################################################################

  #Get task, learner and resampling ids
  lrns_time <- cbind(lrns_dt, params_format) %>%
    .[, learner_id := unlist(lapply(learner, function(x) x$id))] %>%
    .[, task_id := unlist(lapply(task, function(x) x$id))] %>%
    .[, resampling_id := unlist(lapply(resampling, function(x) x$id))]

  #Get number of features
  lrns_time[, npredictors :=
              lapply(lrns_dt$task,function(x) length(x$feature_names)) %>%
              unlist]

  #Compute train and predict time by learner, task and resampling
  bm_ttot <- as.data.table(bmres$aggregate(list(msr('time_train'),
                                                msr('time_predict')))) %>%
    .[lrns_time, on=c('learner_id', 'task_id', 'resampling_id')] %>%
    unique(by='nr') %>%
    .[, time_train := time_train/inner_nruns]

  return(bm_ttot)
}

#------ bm_msrtab ----------------
bm_msrtab <- function(bmres) {
  print('Getting bbrier')
  bbrier_vec <- lapply(1:bmres$n_resample_results, function(i) {
    print(i)
    rsmp_bbrier(bmres=bmres, rsmp_i=i)
  }) %>%
    unlist

  print('Getting auc')
  auc_vec <- lapply(1:bmres$n_resample_results, function(i) {
    print(i)
    rsmp_auc(bmres=bmres, rsmp_i=i)
  }) %>%
    unlist

  print('Getting smp design')
  outer_smp <- lapply(bmres$resamplings$resampling, function(rsmp_design) {
    as.data.table(rsmp_design$param_set$get_values()) %>%
      setnames(names(.), paste0('outer_', names(.))) %>%
      .[, resampling_id := rsmp_design$id]
  }) %>%
    rbindlist(use.names=T)

  print('Getting bacc')
  bacc_vec <- lapply(1:bmres$n_resample_results, function(i) {
    print(paste0('For fold #', i))
    baccthresh <- ldply(seq(0.43,0.55,0.01), function(threshold_class) {
      print(threshold_class)
      rsmp_bacc(bmres, rsmp_i=i, threshold_class=threshold_class)}) %>%
      setorder(bacc) %>%
      setDT %>%
      .[, .SD[.N, .(bacc, threshold_class)]]
  }) %>%
    rbindlist

  print('Getting inner sampling params, general params, and time &  add to rest')
  moddt <- bm_paramstime(bmres) %>%
    cbind(., outer_smp) %>%
    cbind(., bacc_vec) %>%
    .[,`:=`(bbrier = bbrier_vec,
            auc = auc_vec)]

  return(moddt)
}

#------ reset_tuning ----------------------
reset_tuning <- function(in_autotuner, in_task, in_lrnid = NULL) {
  if (inherits(in_autotuner, 'list') & !is.null(in_lrnid)) {
    in_autotuner <- in_autotuner[[
      which(unlist(lapply(in_autotuner, function(lrn) {lrn$id == in_lrnid})))
    ]]
  }

  tuneargs_ini <- in_autotuner$instance_args

  autotuner_new <- set_tuning(in_learner = tuneargs_ini$learner,
                              in_measure = tuneargs_ini$measure,
                              nfeatures = length(in_task$feature_names),
                              insamp_nfolds= tuneargs_ini$resampling$param_set$values$folds,
                              insamp_neval= tuneargs_ini$terminator$param_set$values$n_evals,
                              insamp_nbatch= in_autotuner$tuner$param_set$values$batch_size
  )

  return(autotuner_new)
}
#------ format_modelcompdat --------------
format_modelcompdat <- function(bmres, typecomp=c('classif1', 'regr1', 'classif2')) {
  if (typecomp == 'classif1') {
    bmres[, `:=`(selection = 'Algorithm',
                 type = 'Classif.')]  %>%
      .[, learner_format := dplyr::case_when(
        learner_id == 'classif.ranger'~'default RF',
        learner_id == 'oversample.classif.ranger'~'default RF-oversampled',
        learner_id == 'classweights.classif.ranger'~'default RF-weighted classes',
        learner_id == 'classif.cforest'~'CIF',
        learner_id == 'oversample.classif.cforest'~'CIF-oversampled',
        learner_id == 'classweights.classif.cforest'~'CIF-weighted classes',
      )]

  } else if (typecomp == 'regr1') {
    bmres[, `:=`(selection = 'Algorithm',
                 type = 'Regr.')] %>%
      .[, learner_format := dplyr::case_when(
        task_id == 'inter_regr' ~'MAXSTAT',
        task_id == 'inter_regrover' ~'MAXSTAT-oversampled'
      )]

  } else if (typecomp == 'classif2') {
    bmres[, `:=`(selection = 'Predictors',
                 type = 'Classif.')] %>%
      .[, learner_format := dplyr::case_when(
        task_id == 'inter_class' ~ 'default RF-oversampled-all variables',
        task_id == 'inter_class_featsel' ~ 'default RF-oversampled-selected variables'
      )]
  } else {
    stop('typecomp is not recognized')
  }
}

#------ sfformat_wintri ----------------------
#Convert st to sf and transform to Winkel Trippel projection
sfformat_wintri <- function(in_sp) {
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"

  return(st_as_sf(in_sp) %>%
           st_transform(crs_wintri, use_gdal = FALSE)
  )
}



#------ bin_dt -----------------
#Bin a data_table on a given variable - see tabulate global summary
bin_dt <- function(in_dt, binvar, valuevar, binfunc, binarg,
                   bintrans=NULL, na.rm=FALSE) {
  #Inspired from rbin, adapted to dt and simplified
  in_dt <- copy(in_dt)

  el_freq <- function(byd, bins) {

    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    append(min(byd, na.rm = TRUE), min(byd, na.rm = TRUE) + (bin_length * seq_len(bins)))[1:bins]

  }

  eu_freq <- function(byd, bins) {

    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    ufreq      <- min(byd, na.rm = TRUE) + (bin_length * seq_len(bins))
    n          <- length(ufreq)
    ufreq[n]   <- max(byd, na.rm = TRUE) + 1
    return(ufreq)

  }

  binvar_orig <- copy(binvar)

  #Remove NAs
  if (na.rm) {
    in_dt <- in_dt[!is.na(get(eval(valuevar))),]
  }

  #Transform data if trans
  if (!is.null(bintrans)) {
    transvar <- paste0(binvar, '_bintrans')
    if (bintrans == 'log') {
      nneg <- in_dt[get(eval(binvar)) <= 0, .N]
      warning(paste0('There are ', nneg, ' records with', binvar, ' <= 0...',
                     'removing them for log transformation'))
      in_dt[, eval(transvar) :=
              log(get(eval(binvar)))]

    } else if (is.numeric(bintrans)) {
      in_dt[, eval(transvar) :=
              get(eval(binvar))^eval(bintrans)]

    }
    binvar = transvar
  }

  byd <- in_dt[, get(eval(binvar))]

  if (binfunc == 'manual') {
    l_freq    <- append(min(byd), binarg)
    u_freq    <- c(binarg, (max(byd, na.rm = TRUE) + 1))
    bins      <- length(binarg) + 1
  }

  if (binfunc == 'equal_length') {
    bins = round(binarg)
    l_freq    <- el_freq(byd, bins)
    u_freq    <- eu_freq(byd, bins)
  }

  if (binfunc == 'equal_freq') {
    bins = round(binarg)
    bin_prop     <- 1 / bins
    bin_length   <- in_dt[, round(.N/bins)]
    first_bins   <- (bins - 1) * bin_length
    residual     <- in_dt[, .N - first_bins]
    bin_rep      <- c(rep(seq_len((bins - 1)), each = bin_length),
                      rep(residual, residual))
    l_freq        <- c(1, (bin_length * seq_len((bins - 1)) + 1))
    u_freq       <- c(bin_length * seq_len((bins - 1)), in_dt[,.N])
    setorderv(in_dt, cols= eval(binvar))
    in_dt[, binid := .I]
    binvar = 'binid'
  }

  for (i in seq_len(bins)) {
    in_dt[get(eval(binvar)) >= l_freq[i] & get(eval(binvar)) < u_freq[i],
          bin := i]
    in_dt[bin == i, `:=`(bin_lmin = min(get(eval(binvar_orig)), na.rm=T),
                         bin_lmax = max(get(eval(binvar_orig)), na.rm=T))] %>%
      .[bin == i, bin_lformat := paste(bin_lmin, bin_lmax, sep='-')]

    if (i == bins) {
      in_dt[get(eval(binvar)) == u_freq[i],  bin := i]
    }
  }

  if (binfunc == 'equal_freq') {in_dt[, binid := NULL]}

  return(in_dt)
}

#------ label_manualbins ------------
label_manualbins <- function(binarg, minval) {
  minlabel <- paste(minval, binarg[1], sep=" - ")
  otherlabels <- mapply(function(x, y) {paste(x, y-1, sep=" - ")},
                        binarg[1:(length(binarg)-1)], binarg[2:length(binarg)])
  return(c(minlabel, otherlabels))
}

#------ formathistab -----------------
#Format summary table
formathistab <- function(in_dt, castvar, valuevar, valuevarsub,
                         weightvar, binfunc, binarg, binlabels,
                         datname) {
  rivbin <- bin_dt(in_dt = as.data.table(in_dt),
                   binvar = castvar,
                   valuevar = valuevar,
                   binfunc = binfunc,
                   binarg = binarg)

  netstat <- as.data.table(rivbin)[,sum(get(weightvar)),
                                   by=c(eval(valuevar), "bin")]
  tidyperc_riv <- netstat[, list(perc = 100*.SD[get(valuevar)==valuevarsub,
                                                sum(V1)]/sum(V1),
                                 binsumlength = sum(V1)),
                          by=c('bin')] %>%
    setorder(bin) %>%
    .[, `:=`(dat=datname,
             binformat = binlabels[bin])]
  return(tidyperc_riv)
}

#------ scientific_10 ------------
scientific_10 <- function(x) {
  parse(text=gsub(".*e", "10^", scales::scientific_format()(x)))
}
#------ ggcompare -------------------
#Plot comparing intermittency prevalence and network length by drainage area for two datasets
ggcompare <- function(datmerge, binarg) {
  x_tick <- c(0, unique(datmerge$bin)) + 0.5
  binarg_tick <- c(0, binarg)
  len <- length(x_tick)

  plot_size <- ggplot(datmerge, aes(x=bin, y=binsumlength, fill=dat)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(name = 'Dataset', values=c('#a6cee3', '#1f78b4')) +
    scale_x_continuous(breaks = c(sort(unique(datmerge$bin)), x_tick)[seq(1, 2*len, 2)],
                       labels = c(rep(c(""), len), binarg)[seq(1, 2*len, 2)]) +
    scale_y_continuous(trans='log1p',
                       breaks = c(1, 1000, 100000),
                       labels = scientific_10) +
    labs(x= '',
         y='River length (km)') +
    coord_cartesian(expand=FALSE, clip="off") +
    theme_classic() +
    theme(legend.position = 'none',
          text = element_text(size=9),
          plot.background = element_blank(),
          axis.text.x = element_text(),
          axis.ticks.x = element_line(color = c(rep(NA, len/2 - 1),
                                                rep("black", len/2))))

  plot_inter <- ggplot(datmerge, aes(x=bin, y=perc, fill=dat)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(name = 'Dataset', values=c('#a6cee3', '#1f78b4')) +
    scale_x_continuous(breaks = c(sort(unique(datmerge$bin)), x_tick),
                       labels = c(rep(c(""), len), binarg[1:length(unique(datmerge$bin))])) +
    coord_cartesian(expand=FALSE, ylim=c(0, 1.4*max(datmerge$perc)), clip="off") +
    labs(x= bquote('Discharge'~(m^3~s^-1)),
         y='Prevalence of intermittency (% river length)') +
    theme_classic() +
    theme(legend.position = c(0.15, 0.90),
          legend.background = element_blank(),
          axis.ticks.x = element_line(color = c(rep(NA, len - 1), rep("black", len))))


  xmin_inset <- datmerge[, 0.4*(max(bin)-min(bin))]
  xmax_inset <- datmerge[, 1*max(bin)]
  ymin_inset <- datmerge[,round((0.8 * max(perc, na.rm=T) - min(perc, na.rm=T)) +
                                  min(perc, na.rm=T))]
  plot_join <- plot_inter +
    annotation_custom(grob = ggplotGrob(plot_size),
                      xmin = xmin_inset,
                      xmax = xmax_inset,
                      ymin = ymin_inset,
                      ymax = Inf)
  return(plot_join)
}

#------ bacc_manual ------

##### -------------------- Workflow functions ---------------------------------
#------ def_filestructure -----------------
#' Define file structure
#'
#' Defines file structure for running RF prediction of global intermittency
#'
#' @return Named list with paths to directories and files as well as paths for
#' analysis output
#'
#' @details all paths are relative to file in which project/package is contained. \cr
#' Project must be contained in ~/project_name/src \cr
#' Data must be contained in ~/project_name/data \cr
#' A directory called ~/project_name/results must exist to write outputs
#'
#' @export
#'
#'

def_filestructure <- function() {
  # Get main directory for project
  rootdir <- find_root(has_dir("src"))
  #Directory where raw data are located
  datdir <-  file.path(rootdir, 'data')
  # Directory where results and figures are written
  resdir <- file.path(rootdir, 'results')
  # File geodatabase to write outputs
  outgdb <- file.path(resdir, 'spatialoutputs.gdb')


  # GRDC hydrometric stations that have been joined to RiverATLAS
  in_GRDCgaugep <- file.path(outgdb, 'grdcstations_cleanjoin')
  # Directory containing GRDC hydrometric data
  in_GRDCgaugedir <-  file_in(!!file.path(datdir, 'GRDCdat_day'))
  # GSIM hydrometric stations that have been joined to RiverATLAS
  in_GSIMgaugep <- file.path(resdir, 'GSIM', 'GSIM.gdb',
                             'GSIMstations_cleanriverjoin')
  # Directory containing GSIM monthly hydrometric data
  GSIMgaugedir <- file_in(!!file.path(
    datdir, 'GSIM', 'GSIM_indices', 'TIMESERIES', 'monthly'))

  # River atlas formatted variables
  in_riveratlas_meta <- file_in(!!file.path(datdir, 'HydroATLAS',
                                            'HydroATLAS_metadata_MLMv11.xlsx'))
  #Average monthly discharge for HydroSHEDS network
  in_monthlynetdischarge <- file.path(datdir, 'HydroSHEDS',
                                      'HS_discharge_monthly.gdb',
                                      'Hydrosheds_discharge_monthly')
  #Rasters of dissolved buffers around gauge stations
  in_bufrasdir <- file.path(resdir, 'bufrasdir')
  # River atlas attribute data1
  in_riveratlas <- file_in(!!file.path(resdir, 'RiverATLAS_v10tab.csv'))

  # River atlas attribute data2
  in_riveratlas2 <- file_in(!!file.path(resdir, 'RiverATLAS_v11tab.csv'))

  # French river network for comparison
  compresdir <- file.path(resdir, 'Comparison_databases')
  in_netfr <- file.path(compresdir, 'france.gdb', 'network')
  in_basfr <- file.path(compresdir, 'france.gdb', 'hydrobasins12')

  # Output geopackage of hydrometric stations with appended predicted intermittency class
  out_gauge <- file.path(resdir, 'GRDCstations_predbasic800.gpkg')
  # River atlas predictions table
  out_riveratlaspred <- file.path(resdir, 'RiverATLAS_predbasic800.csv')

  return(c(rootdir=rootdir, datdir=datdir, resdir=resdir, outgdb=outgdb,
           in_gaugep=in_gaugep, in_gaugedir=in_gaugedir,
           in_riveratlas_meta=in_riveratlas_meta,
           in_monthlynetdischarge = in_monthlynetdischarge,
           in_bufrasdir = in_bufrasdir,
           in_netfr = in_netfr, in_basfr = in_basfr,
           out_gauge=out_gauge,
           in_riveratlas=in_riveratlas, in_riveratlas2=in_riveratlas2,
           out_riveratlaspred = out_riveratlaspred))
}

#------ read_monthlydis -------------
#' Read network monthly discharge
#'
#' Import Naturalized mean monthly discharge from WaterGAP v2.2 downscaled to
#' the HydroSHEDS river network (12 monthly values for each river reach)
#'
#' @param in_filestructure named list containing path to RiverATLAS spatial data
#' with monthly discharge, named \code{in_monthlynetdischarge}.
#'
#' @return Data frame of monthly discharge for every reach in RiverATLAS
#'
#' @export
#'
read_monthlydis <- function(in_path) {
  monthlydischarge <- as.data.frame(st_read(
    dsn = dirname(in_path),
    layer = basename(in_path)
  )) %>%
    .[, c(paste0('DIS_', str_pad(seq(1, 12), 2, pad=0), '_CMS'),
          'REACH_ID')]
  return(monthlydischarge)
}

#------ read_gaugep -----------------
#' Read gauge points
#'
#' Import streamgauging station spatial data and attributes, only keeping those
#' whose point representation is within a given distance from the river network.
#'
#' @param in_filestructure named list containing path to point data for gauging
#' stations, named \code{in_gaugep}
#' @param dist maximum distance from the river network beyond which gauging
#' stations are excluded.
#' @param in_monthlydischarge data frame of naturalized monthly discharge for
#' every river reach in HydroSHEDS
#'
#' @return object of class sf
#'
#' @details The distance from the river network must have been determined
#' beforehand and be an attribute of the gauge stations point data called
#' \code{station_river_distance}
#'
#' @export

read_gaugep <- function(inp_GRDCgaugep, inp_GSIMgaugep,
                        inp_riveratlas2, in_monthlydischarge=NULL) {
  #Import gauge stations
  GRDCgaugep <- st_read(dsn = dirname(inp_GRDCgaugep),
                        layer = basename(inp_GRDCgaugep))

  GSIMgaugep <- st_read(dsn = dirname(inp_GSIMgaugep),
                        layer = basename(inp_GSIMgaugep))

  #Rename
  GRDCgaugep$GAUGE_NO <- GRDCgaugep$GRDC_NO
  GSIMgaugep$GAUGE_NO <- GSIMgaugep$gsim_no
  GSIMgaugep$manualsnap_mathis <- GSIMgaugep$manualsnap
  GSIMgaugep$snap_comment_mathis <- GSIMgaugep$snap_comment

  #Row bind GSIM and GRDC stations
  GRDCgaugep$area_correct <- GRDCgaugep$GRDC_AREA
  GRDCgaugep$gsim_no <- NA
  GSIMgaugep$GRDC_NO <- NA
  GRDCgaugep[, c("gsim_no", "reference_db", "reference_no", "grdb_merge",
                 "grdb_no", "paired_db", "paired_db_no")] <- NA
  keepcols_bind <- intersect(names(GSIMgaugep), names(GRDCgaugep))

  gaugep <- rbind(GRDCgaugep[, keepcols_bind],
                  GSIMgaugep[, keepcols_bind])

  #Get new and updated environmental predictors from RiverATLAS v1.0.9
  riveratlas2 <- fread(inp_riveratlas2)
  setnames(riveratlas2,
           names(riveratlas2),
           gsub('_11$', '', names(riveratlas2)))

  #Replace variables from RiverATLAS v1.0 by those that have been re-calculated
  #in v1.0.9
  keepcols <- names(gaugep)[!(names(gaugep) %in% names(riveratlas2))]
  gaugep_attriall <- merge(gaugep[,keepcols, with=F], riveratlas2,
                           by.x = 'HYRIV_ID', by.y = 'REACH_ID',
                           all.x=TRUE, all.y=FALSE)

  #Merge with WaterGAP downscaled monthly naturalized discharge
  if (!is.null(in_monthlydischarge)) {
    gaugep_monthlydischarge <- merge(gaugep_attriall, in_monthlydischarge,
                                     by.x = 'HYRIV_ID', by.y = 'REACH_ID',
                                     all.x=TRUE, all.y=FALSE)
  } else {
    gaugep_monthlydischarge <- gaugep_attriall
  }


  return(gaugep_monthlydischarge)
}

#------ read_GRDCgauged_paths -----------------
#' Read file paths to gauge flow data
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @inheritParams read_gaugep
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables
#'
#' @export

read_GRDCgauged_paths <- function(inp_GRDCgaugedir, in_gaugep) { #, gaugeid = 'GRDC_NO' down the line
  #Get data paths of daily records for gauge stations
  fileNames <- file.path(inp_GRDCgaugedir,
                         paste(
                           in_gaugep[!is.na(in_gaugep$GRDC_NO),]$GRDC_NO,
                           ".txt", sep=""))
  #Check whether any GRDC record does not exist
  print(paste(length(which(!do.call(rbind, lapply(fileNames, file.exists)))),
              'GRDC records do not exist...'))
  return(fileNames)
}

#------ read_GSIMgauged_paths -----------------
#' Read file paths to gauge flow data
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @inheritParams read_gaugep
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables
#'
#' @export

read_GSIMgauged_paths <- function(inp_GSIMindicesdir, in_gaugep, timestep) {
  gaugeno_vec <- in_gaugep[!is.na(in_gaugep$gsim_no),]$gsim_no

  #Get data paths of daily records for gauge stations
  if (timestep == 'month') {
    fileNames <- file.path(inp_GSIMindicesdir, 'TIMESERIES', 'monthly',
                           paste0(gaugeno_vec, ".mon"))
  } else if (timestep == 'season') {
    fileNames <- file.path(inp_GSIMindicesdir, 'TIMESERIES', 'seasonal',
                           paste0(gaugeno_vec, ".seas"))
  }

  #Check whether any GRDC record does not exist
  print(paste(length(which(!do.call(rbind, lapply(fileNames, file.exists)))),
              'GSIM records do not exist...'))
  return(fileNames)
}

#------ comp_GRDCdurfreq -------------------------------
#' Compute intermittency statistics for GRDC gauging stations
#'
#' Determine general characteristics of the whole time series and of the subset
#' of years that have less than a given threshold of missing data as well as
#' intermittency statistics. The intermittency statistics can be computed for a
#' subset of months of the year (e.g. only winter months)
#'
#' @param path file path to a GRDC-formatted streamflow time series table
#' @param maxgap maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param mdurthresh threshold of mean annual number of zero-flow days beyond
#' which to classify gauge as intermittent.
#' @param monthsel selected months to compute the statistics over
#'
#' @return One row data.table with the following columns: \cr
#' \describe{
#' \item{GRDC_NO} - (char) unique identifier for the gauge
#' \item{firstYear} - (num) first year on full record
#' \item{lastYear} - (num) last year on full record
#' \item{totalYears} - (int) total number of years on full record
#' \item{firstYear_kept} - (num) first year on record with < maxgap missing days
#' \item{lastYear_kept} - (num) first year on record with < maxgap missing days
#' \item{totalYears_kept} - (int) total number of years with < maxgap missing days
#' \item{totaldays} - (num) total number of days with discharge data
#' \item{sumDur} - (int) total number of days with discharge = 0
#' \item{mDur} - (num) mean number of days/year with discharge = 0
#' \item{mFreq} - (num) mean number of periods with discharge = 0
#' (with at least one day of flow between periods)
#' }
#'
#' @export

comp_GRDCdurfreq <- function(path, in_gaugep, maxgap, mdurthresh = 1,
                             windowsize = 100, fullwindow = FALSE,
                             monthsel = NULL, verbose = FALSE) {
  if (verbose) {
    print(path)
  }

  #Read and format discharge records
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  gaugetab <- readformatGRDC(path)

  #Function to compute mean zero-flow duration and event frequency by month
  #(average number of drying events in each month over record length)
  #and check whether intermittency only happens in cold months
  comp_monthlyirtemp <- function(gaugetab, gaugeno, in_gaugep, maxgap,
                                 mdurthresh, tempthresh, yearthresh) {
    if (gaugetab[(missingdays < maxgap) & (year >= yearthresh), .N>0]) { #Make sure that there are years with sufficient data
      monthlyfreq <- gaugetab[
        (missingdays < maxgap)  & (year >= yearthresh),
        .(GRDC_NO = unique(GRDC_NO),
          monthrelfreq = length(unique(na.omit(prevflowdate)))/length(unique(year)),
          monthmdur = length(na.omit(prevflowdate))/length(unique(year))),
        by='month'] #Count proportion of years with zero flow occurrence per month
    } else {
      monthlyfreq <- data.table(GRDC_NO = gaugetab[, unique(GRDC_NO)],
                                month=1:12,
                                monthrelfreq=rep(NA,12),
                                monthmdur=rep(NA,12))
    }
    abbrev_months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    monthlyfreq_format <- monthlyfreq %>%
      dcast(GRDC_NO~month, value.var=c('monthrelfreq', 'monthmdur')) %>%
      setnames(c('GRDC_NO',
                 paste0(rep(abbrev_months, times=2), '_',
                        rep(c('mfreq', 'mdur'), each=12)
                 )
      ))

    #Check whether average yearly number of zero-flow days is only above threshold
    # if includes zero-flow days during months with average temperature < 5C
    gaugeirtemp <- as.data.table(in_gaugep)[GRDC_NO == gaugeno,
                                            c('GRDC_NO',
                                              grep('^tmp.*_c[0-9]{2}',
                                                   colnames(in_gaugep),
                                                   value = T)),
                                            with=F] %>%
      melt(id.vars='GRDC_NO') %>%
      .[, month := as.numeric(gsub('^tmp_dc_c', '', variable))] %>%
      merge(monthlyfreq, ., by='month')

    mdur_otempthresh <- gaugeirtemp[value >= 10*tempthresh, sum(monthmdur)]
    mdur_utempthresh <- gaugeirtemp[value < 10*tempthresh, sum(monthmdur)]

    monthlyfreq_format[, winteronlyir := as.numeric(
      mdur_otempthresh < mdurthresh &
        mdur_utempthresh >= mdurthresh)] %>%
      .[, GRDC_NO := NULL]

    setnames(monthlyfreq_format,
             paste0(names(monthlyfreq_format), '_o', yearthresh))

    return(monthlyfreq_format)
  }


  monthlyirtemp_all <- comp_monthlyirtemp(gaugetab, gaugeno, in_gaugep,
                                          maxgap, mdurthresh,
                                          yearthresh=1800, tempthresh=10)
  monthlyirtemp_o1961 <- comp_monthlyirtemp(gaugetab, gaugeno,in_gaugep,
                                            maxgap, mdurthresh,
                                            yearthresh=1961, tempthresh=10)
  monthlyirtemp_o1971 <- comp_monthlyirtemp(gaugetab, gaugeno,in_gaugep,
                                            maxgap, mdurthresh,
                                            yearthresh=1971, tempthresh=10)

  #If analysis is only performed on a subset of months
  if (!is.null(monthsel)) {
    gaugetab <- gaugetab[month %in% monthsel, ]
  }

  #Compute number of days of zero flow for years with number of gap days under threshold
  gaugetab_yearly <- merge(gaugetab[, .(missingdays=max(missingdays, na.rm=T),
                                        datadays=max(datadays, na.rm=T),
                                        integerperc=sum(integervalue,na.rm=T)
  ), by='year'],
  gaugetab[Original == 0, .(dur=.N,
                            freq=length(unique(prevflowdate))), by='year'],
  by = 'year', all.x = T) %>%
    .[!is.finite(missingdays), missingdays := diny(year)] %>%
    .[!is.finite(datadays), datadays := diny(year)] %>%
    .[dur==diny(year), freq:=1] %>%
    .[is.na(dur), `:=`(dur=0, freq=0)]


  #Combine all statistics (and determine which stations are labeled as
  #intermittent based on mdurthresh)
  gaugetab_all <- gaugetab_yearly[, .(GRDC_NO = gaugeno,
                                      firstYear=min(year),
                                      lastYear=max(year),
                                      totalYears=length(unique(year))
  )]

  comp_irstats <- function(tabyearly, maxgap, mdurthresh, yearthresh,
                           windowsize, fullwindow) {
    checkpos <- function(x) {any(x>0)}
    if ((windowsize %% 2)==0) {
      windowsize <- windowsize + 1
    }

    if (!fullwindow) {
      movinginter <- all(
        tabyearly[missingdays <= maxgap & year >= yearthresh,
                  (frollapply(dur, n=windowsize, FUN=checkpos, align="center") >= mdurthresh)],
        na.rm=T)
    } else {
      movinginter <- all(
        tabyearly[missingdays <= maxgap & year >= yearthresh,
                  c((frollapply(dur, n=round(windowsize/2), FUN=checkpos, align="left") >= mdurthresh),
                    (frollapply(dur, n=round(windowsize/2), FUN=checkpos, align="right") >= mdurthresh))],
        na.rm=T)
    }


    irstats <- tabyearly[missingdays <= maxgap & year >= yearthresh,
                         .(firstYear_kept=min(year),
                           lastYear_kept=max(year),
                           totalYears_kept=length(unique(year)),
                           totaldays = sum(datadays),
                           integerperc = sum(integerperc)/(sum(datadays)+sum(missingdays)),
                           sumDur = sum(dur),
                           mDur = mean(dur),
                           mFreq = mean(freq),
                           intermittent =
                             factor(fifelse(mean(dur)>=mdurthresh, 1, 0),
                                    levels=c('0','1')),
                           movinginter = movinginter
                         )]
    setnames(irstats, new = paste0(names(irstats), '_o', yearthresh))
    return(irstats)
  }

  irstats_all <- comp_irstats(tabyearly = gaugetab_yearly, maxgap=maxgap,
                              mdurthresh = mdurthresh,
                              yearthresh = 1800,
                              windowsize = windowsize,
                              fullwindow = fullwindow)
  irstats_1961 <- comp_irstats(tabyearly = gaugetab_yearly, maxgap=maxgap,
                               mdurthresh = mdurthresh,
                               yearthresh = 1961,
                               windowsize = windowsize,
                               fullwindow = fullwindow)
  irstats_1971 <- comp_irstats(tabyearly = gaugetab_yearly, maxgap=maxgap,
                               mdurthresh = mdurthresh,
                               yearthresh = 1971,
                               windowsize = windowsize,
                               fullwindow = fullwindow)

  statsout <- cbind(gaugetab_all,
                    irstats_all, irstats_1961, irstats_1971,
                    monthlyirtemp_all, monthlyirtemp_o1961, monthlyirtemp_o1971)

  #Include local path to discharge records in table
  statsout[, path := path]

  return(statsout)
}


#------ comp_GSIMdurfreq -------------------------------
comp_GSIMdurfreq <- function(path_mo, path_sea,
                             in_gaugep, maxgap, mdurthresh = 1,
                             windowsize = 100, fullwindow = FALSE,
                             monthsel = NULL, verbose = FALSE) {

  #Read and format discharge records and join monthly and seasonal records
  gaugeno <- strsplit(basename(path_mo), '[.]')[[1]][1]
  if (verbose) {
    print(gaugeno)
  }

  gaugetab_mo <- readformatGSIMmon(path_mo)
  gaugetab_sea <- readformatGSIMsea(path_sea)
  gaugetab <- merge(gaugetab_mo,
                    gaugetab_sea[, .(year, month, mDur_minsea)],
                    by=c('year', 'month'), all.x=T) %>%
    setorder(date) %>%
    .[, mDur_minsea := nafill(mDur_minsea, type='nocb')]

  #Function to compute mean zero-flow duration and event frequency by month
  #(average number of drying events in each month over record length)
  #and check whether intermittency only happens in cold months
  comp_monthlyirtemp <- function(gaugetab, gaugeno, in_gaugep, maxgap,
                                 mdurthresh, tempthresh, yearthresh) {
    if (gaugetab[(missingdays < maxgap) & (year >= yearthresh), .N>0]) { #Make sure that there are years with sufficient data
      monthlyfreq <- gaugetab[
        (missingdays < maxgap)  & (year >= yearthresh),
        .(gsim_no = unique(gsim_no),
          monthmdur = mean(mDur_minmo)),
        by='month'] #Count proportion of years with zero flow occurrence per month
    } else {
      monthlyfreq <- data.table(gsim_no = gaugetab[, unique(gsim_no)],
                                month=1:12,
                                monthmdur=rep(NA,12))
    }
    abbrev_months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    monthlyfreq_format <- monthlyfreq %>%
      dcast(gsim_no~month, value.var='monthmdur') %>%
      setnames(c('gsim_no',
                 paste0(abbrev_months, '_',
                        rep('mdur', each=12)
                 )
      ))

    #Check whether average yearly number of zero-flow days is only above threshold
    # if includes zero-flow days during months with average temperature < 5C
    gaugeirtemp <- as.data.table(in_gaugep)[gsim_no == gaugeno,
                                            c('gsim_no',
                                              grep('^tmp.*_c[0-9]{2}',
                                                   colnames(in_gaugep),
                                                   value = T)),
                                            with=F] %>%
      melt(id.vars='gsim_no') %>%
      .[, month := as.numeric(gsub('^tmp_dc_c', '', variable))] %>%
      merge(monthlyfreq, ., by='month')

    mdur_otempthresh <- gaugeirtemp[value >= 10*tempthresh, sum(monthmdur)]
    mdur_utempthresh <- gaugeirtemp[value < 10*tempthresh, sum(monthmdur)]

    monthlyfreq_format[, winteronlyir := as.numeric(
      mdur_otempthresh < mdurthresh &
        mdur_utempthresh >= mdurthresh)] %>%
      .[, gsim_no := NULL]

    setnames(monthlyfreq_format,
             paste0(names(monthlyfreq_format), '_o', yearthresh))

    return(monthlyfreq_format)
  }


  monthlyirtemp_all <- comp_monthlyirtemp(gaugetab, gaugeno, in_gaugep,
                                          maxgap, mdurthresh,
                                          yearthresh=1800, tempthresh=10)
  monthlyirtemp_o1961 <- comp_monthlyirtemp(gaugetab, gaugeno,in_gaugep,
                                            maxgap, mdurthresh,
                                            yearthresh=1961, tempthresh=10)
  monthlyirtemp_o1971 <- comp_monthlyirtemp(gaugetab, gaugeno,in_gaugep,
                                            maxgap, mdurthresh,
                                            yearthresh=1971, tempthresh=10)

  #If analysis is only performed on a subset of months
  if (!is.null(monthsel)) {
    gaugetab <- gaugetab[month %in% monthsel, ]
  }


  #Combine all statistics (and determine which stations are labeled as
  #intermittent based on mdurthresh)
  gaugetab_all <- gaugetab[, .(gsim_no = gaugeno,
                               firstYear=min(year),
                               lastYear=max(year),
                               totalYears=length(unique(year))
  )]

  comp_irstats <- function(tab, maxgap, mdurthresh, yearthresh,
                           windowsize, fullwindow) {
    #Compute the best estimate of minimum number of zero flow days per season
    mDur_final <- tab[, .(mDur_minsea_final = max(mDur_minsea,
                                                  sum(mDur_minmo, na.rm=T),
                                                  na.rm = T)
    ), by=c('year', 'season')] %>%
      .[, .(mDur_minsea_final = sum(mDur_minsea_final, na.rm=T)), by=year]

    yearsel <- tab[missingdays <= maxgap & year >= yearthresh, unique(year)]

    checkpos <- function(x) {any(x>0)}
    if ((windowsize %% 2)==0) {
      windowsize <- windowsize + 1
    }

    if (!fullwindow) {
      movinginter <- all(
        mDur_final[year %in% yearsel, checkpos(mDur_minsea_final), by=year] %>%
          .[, frollapply(V1, n=windowsize, FUN=checkpos, align="center") >= mdurthresh],
        na.rm = T
      )
    } else {
      movinginter <- all(
        mDur_final[year %in% yearsel, checkpos(mDur_minsea_final), by=year] %>%
          .[, c(
            frollapply(V1, n=round(windowsize/2), FUN=checkpos, align="left") >= mdurthresh,
            frollapply(V1, n=round(windowsize/2), FUN=checkpos, align="right") >= mdurthresh
          )],
        na.rm = T
      )
    }

    irstats <- tab[year %in% yearsel,
                   .(firstYear_kept=min(year),
                     lastYear_kept=max(year),
                     totalYears_kept=length(unique(year)),
                     totaldays = sum(n.available),
                     sumDur = mDur_final[year %in% yearsel,
                                         sum(mDur_minsea_final)],
                     mDur = mDur_final[year %in% yearsel,
                                       mean(mDur_minsea_final),],
                     movinginter = movinginter
                   )] %>%
      .[, intermittent := factor(fifelse(mDur>=mdurthresh, 1, 0),
                                 levels=c('0','1'))]

    setnames(irstats, new = paste0(names(irstats), '_o', yearthresh))
    return(irstats)
  }

  irstats_all <- comp_irstats(tab = gaugetab, maxgap=maxgap,
                              mdurthresh = mdurthresh,
                              yearthresh = 1800,
                              windowsize = windowsize,
                              fullwindow = fullwindow)
  irstats_1961 <- comp_irstats(tab = gaugetab, maxgap=maxgap,
                               mdurthresh = mdurthresh,
                               yearthresh = 1961,
                               windowsize = windowsize,
                               fullwindow = fullwindow)
  irstats_1971 <- comp_irstats(tab = gaugetab, maxgap=maxgap,
                               mdurthresh = mdurthresh,
                               yearthresh = 1971,
                               windowsize = windowsize,
                               fullwindow = fullwindow)

  statsout <- cbind(gaugetab_all,
                    irstats_all, irstats_1961, irstats_1971,
                    monthlyirtemp_all, monthlyirtemp_o1961, monthlyirtemp_o1971)

  #Include local path to discharge records in table
  statsout[, path := path_mo]

  return(statsout)
}

#------ plot_GRDCflags ------
plot_GRDCflags <- function(in_GRDCgaugestats, yearthresh,
                           inp_resdir, maxgap, showmissing = FALSE) {
  if (inherits(in_GRDCgaugestats, 'data.table')) {
    GRDCstatsdt <- in_GRDCgaugestats
  }  else if (is.list(in_GRDCgaugestats)) {
    GRDCstatsdt <- rbindlist(in_GRDCgaugestats)
  }

  #Create output directory for IRs
  resdir_GRDCirplots <- file.path(inp_resdir,
                                  paste0('GRDCir_rawplots_', yearthresh, '_',
                                         format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_GRDCirplots))) {
    print(paste0('Creating ', resdir_GRDCirplots ))
    dir.create(resdir_GRDCirplots )
  }

  #Plot
  GRDCstatsdt[(get(paste0('totalYears_kept_o', yearthresh)) >= 10) &
                (get(paste0('intermittent_o', yearthresh)) == 1),
              plotGRDCtimeseries(.SD,
                                 outpath = file.path(resdir_GRDCirplots,
                                                     paste0(GRDC_NO, '.png')),
                                 maxgap=maxgap,
                                 showmissing = showmissing
              ), by=GRDC_NO]

  #Create output directory for non IRs
  resdir_GRDCperplots <- file.path(inp_resdir,
                                   paste0('GRDCper_rawplots_',yearthresh, '_',
                                          format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_GRDCperplots))) {
    print(paste0('Creating ', resdir_GRDCperplots ))
    dir.create(resdir_GRDCperplots )
  }

  GRDCstatsdt[(get(paste0('totalYears_kept_o', yearthresh)) >= 10) &
                (get(paste0('intermittent_o', yearthresh)) == 0),
              plotGRDCtimeseries(.SD,
                                 outpath = file.path(resdir_GRDCperplots,
                                                     paste0(GRDC_NO, '.png')),
                                 maxgap=maxgap,
                                 showmissing = showmissing
              ), by=GRDC_NO]


}

#------ plot_GSIM -------------
plot_GSIM <- function(in_GSIMgaugestats, yearthresh,
                      inp_resdir, maxgap, showmissing) {

  if (inherits(in_GSIMgaugestats, 'data.table')) {
    GSIMstatsdt <- in_GSIMgaugestats
  }  else if (is.list(in_GSIMgaugestats)) {
    GSIMstatsdt <- rbindlist(in_GSIMgaugestats)
  }

  #Create output directory for IRs
  resdir_GSIMirplots <- file.path(inp_resdir,
                                  paste0('GSIMir_rawplots_',yearthresh, '_',
                                         format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_GSIMirplots))) {
    print(paste0('Creating ', resdir_GSIMirplots ))
    dir.create(resdir_GSIMirplots )
  }

  #Plot
  GSIMstatsdt[(get(paste0('totalYears_kept_o', yearthresh)) >= 10) &
                (get(paste0('intermittent_o', yearthresh)) == 1),
              plotGSIMtimeseries(.SD,
                                 outpath = file.path(resdir_GSIMirplots,
                                                     paste0(gsim_no, '.png')),
                                 maxgap=maxgap,
                                 showmissing=showmissing
              ), by=gsim_no]


  #Create output directory for non IRs
  resdir_GSIMperplots <- file.path(inp_resdir,
                                   paste0('GSIMper_rawplots_',yearthresh, '_',
                                          format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_GSIMperplots))) {
    print(paste0('Creating ', resdir_GSIMperplots ))
    dir.create(resdir_GSIMperplots )
  }

  GSIMstatsdt[(get(paste0('totalYears_kept_o', yearthresh)) >= 10) &
                (get(paste0('intermittent_o', yearthresh)) == 0),
              plotGSIMtimeseries(.SD,
                                 outpath = file.path(resdir_GSIMperplots,
                                                     paste0(gsim_no, '.png')),
                                 maxgap=maxgap,
                                 showmissing=showmissing
              ), by=gsim_no]
}

#------ analyze_gaugeir ----------------------------
#Check winter and coastal gauges in particular
#Check availability of data
analyzemerge_gaugeir <- function(in_GRDCgaugestats, in_GSIMgaugestats, yearthresh,
                                 in_gaugep, inp_resdir, plotseries = FALSE) {
  ### Analyze GSIM data ####################################
  GSIMstatsdt <- rbindlist(in_GSIMgaugestats)

  #Outliers from examining plots of ir time series (those that were commented out were initially considered)
  GSIMtoremove_o1800_irartifacts <- c(
    'AT_0000014', #remove - changed flow permanence
    'BR_0000286', #remove - Tapajos river. Impossible that it dries out.
    'BR_0000581', #remove - 0 flow happened only once at the end of record
    'BR_0000786', #remove - probably perennial
    'BR_0000862', #remove - probably perennial
    'BR_0001011', #remove - probably perennial
    'BR_0001116', #remove - probably naturally perennial
    'BR_0001133', #remove - probably naturally perennial
    'CA_0001057', #remove - only one 0 flow event in 28 years
    'CA_0003488', #remove - only one 0 flow event
    'CA_0003526', #remove - only one 0 flow event
    'ES_0000078', #remove - first 20 years without 0 flow
    'ES_0000084', #remove - probably perennial, first 20 years without 0 flow
    'ES_0000388', #remove - record is too discontinued to assess permanence, probably perennial; besides, it looks regulated
    'ES_0000444', #remove - only one 0 valid flow occurrence
    'ES_0000525', #remove - very discontinued record. Seems to have changed flow permanence
    'ES_0000729', #remove - seems to have shifted
    'ES_0000785', #remove - erroneous
    'ES_0000786', #remove - looks regulated
    'ES_0000816', #remove - doesn't look reliable
    'ES_0000856', #remove - 0 flow occured only once in record
    'ES_0000892', #remove - probably perennial
    'ES_0000806', #remove - probably perennial, 0 flow occurred only once
    'ES_0000910', #remove - changed flow permanence
    'ES_0000986', #remove - 0 flow occured only once
    'ES_0000996', #remove - changed flow permanence
    'ES_0001020', #remove - hard to tell flow permanence
    'ES_0001052', #remove - regulated
    'ES_0001082', #remove - changed flow permanence, either regulated or...
    'ES_0001085', #remove - changed flow permanence
    'ES_0001116', #remove - occurred only the last year not the first 22 years
    'ES_0001162', #remove - most occurrences seem buggy
    'FI_0000156', #remove - no occurrences in first 60% of record
    'FR_0000052', #remove - no occurrences in first 20 years
    'HU_0000017', #remove - only one zero flow occurence in 25 years
    'IE_0000014', #remove - no 0 flow occurrence in first 25 years
    'IN_0000023', #remove - no 0 flow occurrence in first 20 years
    'IN_0000045', #remove - no 0 flow occurrence in first 20 years
    'IN_0000046', #remove - changed flow permanence
    'IN_0000050', #remove - changed flow permanence
    'IN_0000063', #remove - changed flow permanence
    'IN_0000064', #remove - changed flow permanence
    'IN_0000113', #remove - only one 0 flow occurrence
    'IN_0000124', #remove - changed flow permanence
    'IN_0000125', #remove - changed flow permanence
    'IN_0000134', #remove - changed flow permanence
    'IN_0000170', #remove - only one 0 flow occurrence
    #'IN_0000190', #checked- looked regulated
    'IN_0000255', #remove - changed flow permanence
    'IN_0000283', #remove - changed flow permanence
    'IN_0000312', #remove - changed flow permanence
    'MZ_0000010', #remove - changed flow permanence
    'NO_0000028', #remove - changed flow permanence
    'NO_0000030', #remove - occurred only once
    'RU_0000189', #remove - occurred only once
    'RU_0000265', #remove - changed flow permanence
    'RU_0000358', #remove - changed flow permanence
    'TZ_0000018', #remove - hard to tell but probably changed flow permanence
    'US_0005106', #remove - hard to tell, probably changed flow permanence
    'US_0005874', #remove - changed flow permanence
    'US_0005876', #remove - changed flow permanence
    'US_0006103', #remove -  regulated
    #'US_0006156', #checked -  regulated, barrier made it perennial?
    'US_0006396', #remove - buggy looking pre-1960
    'US_0006440', #remove - not enough data to tell, seemed like changed flow permanence
    'US_0006537', #remove - changed flow permanence
    'ZA_0000074', #remove - missing data
    'ZA_0000268', #remove - seemingly changed flow permanence
    'ZA_0000270' #remove - changed flow permanence
  )

  GSIMtoremove_o1961_irartifacts <- c(
    #'AR_0000014', #maybe regulated --- checked
    'AT_0000021', #remove
    'AT_0000026', #remove
    'AT_0000038', #remove
    'AT_0000059', #remove
    'BR_0000557', #inspected - confirmed dry channel visually
    'BR_0000662', #is regulated -- remove
    #'BR_0000664', #change of regime, maybe regulated --- but was intermittent before
    'BR_0000706', #remove
    'BR_0000717', #remove --regulated intermittency
    'BR_0000778', #remove, too much missing data
    'BR_0001104', #remove
    #'CA_0003473', downstream of natural lake, all good
    #'CA_0003526', #downstream of lakes, maybe low lake levels
    'CA_0003544', #already in GRDC and a bit buggy looking
    'CN_0000004', #remove
    'CN_0000009', #remove
    'CN_0000010', #remove
    'CN_0000012', #remove
    'CN_0000013', #remove
    'CN_0000021', #remove
    'CN_0000022', #remove
    'CN_0000026', #remove
    'CN_0000029', #remove
    'CN_0000038', #remove - got dewatered by dam
    'CN_0000043', #remove
    'CN_0000062', #remove
    #'ES_0000525', #maybe regulated --- check
    'ES_0000581', #remove
    'ES_0000660', #remove -- regulated
    #'ES_0000676', #looks regulated --- but probably intermittent before
    'ES_0000794', #remove - same as 0000832
    'ES_0000784', #remove -- same as 'ES_0000785'
    'ES_0000818', #remove - same as ES_0000785 and ES_0000784 which are not intermittent
    'ES_0000841', #remove -- regulated
    'ES_0000958', #remove -- 0s seem erroneous
    'FI_0000107', #remove -- lake inlet, maybe just become standing water when high water level
    #'IE_0000014', #downstream of lake -- must be drop in lake levels
    'IN_0000105', #remove -- downstream of major dams -- no records pre-1968 time of building
    'IN_0000159', #remove -- regulated, no records before dam
    'IN_0000280', #remove -- regulated, no records before dam. intermittency started after building of Parambikulam dam
    #'IN_0000309', #intermittent before building of reservoir in 1988
    'NA_0000050', #remove, unreliable, 0s and interpolations
    #'NO_0000018', #checked -- downstream of natural lake. probably water level decrease
    'NO_0000044', #remove
    'NO_0000090', #remove
    'SE_0000058', #remove -- downstream of dam, didn't use to be intermittent. and experiences ice
    #'US_0001855', #maybe rounded values --- check
    #'US_0001861', #maybe rounded values --- check
    #'US_0001868', #maybe rounded avlues --- check
    'US_0002247', #remove -- regulated -- GRDC 4149415 upstream not intermittent
    'US_0002248', #remove -- regulated just downstream of 2247
    'US_0002791', #on usgs website: "close proximity to the Ohio River. During periods of high water on the Ohio River, the computed discharge at this site may be incorrectly displayed due to the backwater effect created."
    #'US_0003591', #maybe rounded values --- check
    # 'US_0003774', #maybe rounded values --- check
    # 'US_0003836', #maybe rounded values --- check
    # 'US_0004023', #maybe rounded values --- check
    # 'US_0004216', #maybe rounded values --- check
    # 'US_0004232', #maybe rounded values --- check
    # 'US_0004658', #maybe rounded values --- check
    'US_0004773', #remove - regulated. no records prior to reservoir buiding
    #US_005099, #confirmed dried bed on imagery
    # 'US_0005161', #maybe rounded values --- check
    # 'US_0005177', #maybe rounded values --- check
    #'US_0005303', #looks fine on usgs website and imagery. just small
    #'US_0005596', #looks fine, small
    #'US_0005597', #looks fine, small
    # 'US_0005622', #maybe rounded values --- check
    # 'US_0005623', #maybe rounded values --- check
    # 'US_0005684', #maybe rounded values --- check
    # 'US_0005687', #maybe rounded values --- check
    'US_0005732', #remove -- regime shift because of Lake Arcadia (reservoir)
    #'US_0005859', #maybe rounded values --- check
    #'US_0005879', #maybe rounded values --- check
    # 'US_0006073', #maybe rounded values --- check
    # 'US_0006109', #maybe rounded values --- check
    # 'US_0006154', #maybe rounded values --- check
    'US_0006155', #currently regulated but used to be intermittent but same as 6156 so remove
    #'US_0006156', #currently regulated but used to be intermittent so keep
    #'US_0006301', #maybe rounded values --- check
    #'US_0006327', #maybe rounded values --- check
    #'US_0006387', #maybe rounded values --- check
    #'US_0006206', #all good, confirmed dry with imagery
    #'US_0006574', #maybe regulated? ---checked all good
    #'US_0006975', #maybe rounded values --- check
    #'US_0006984', #maybe rounded values --- check
    #'US_0006985', #maybe rounded values --- check
    #'US_0006986', #maybe rounded values --- check
    'US_0008607', #regulated. 1978 intermittency started after building reservoir
    'US_0008687', #regulated. intermittency started after bui;lding reserovir 1935
    # 'US_0008726', #maybe rounded values --- check
    # 'US_0008779', #maybe rounded values --- check
    'ZA_0000008', #remove
    'ZA_0000084' #remove
  )

  # readformatGSIMmon(
  #   GSIMstatsdt[gsim_no == 'US_0000546', path]) %>%
  #   .[MIN != 0, min(MIN)]


  #-----  Check flags in winter IR for GSIM
  wintergaugesall_GSIM <- plot_winterir(
    dt = GSIMstatsdt, dbname = 'gsim', inp_resdir = inp_resdir,
    yearthresh = 1800, plotseries = plotseries)

  #Check suspicious canadian ones
  GSIMwintermeta <- in_gaugep[in_gaugep$gsim_no %in% wintergaugesall_GSIM$gsim_no,]
  canadians_toinspect <- in_gaugep[in_gaugep$gsim_no %in%
                                     paste0('CA_000', c(3469, 3473, 3526, 3544, 6082, 6122)),]$reference_no

  # if (!dir.exists(hy_dir())) download_hydat()
  # cancheck <- lapply(canadians_toinspect, function(refno) {
  #   merge(hy_daily(station_number = refno),
  #         hy_stn_regulation(station_number = refno),
  #         by='STATION_NUMBER') %>%
  #     setDT
  # }) %>%
  #   rbindlist
  # cancheck[REGULATED==T, unique(STATION_NUMBER)] #No regulated station
  # cancheck[Value==0, .N, by=.(STATION_NUMBER, Symbol)]

  #E - Estimate:  no measured data available for the day or missing period,
  #    and the water level or streamflow value was estimated by an indirect method
  #A - Partial Day:  daily mean value was estimated despite gaps of more than
  #    120 minutes in the data string or missing data not significant enough to
  #    warrant the use of the E symbol.
  #B - Ice conditions: value was estimated with consideration for the presence
  #    of ice in the stream. Ice conditions alter the open water relationship
  #    between water levels and streamflow.
  #D - Dry: stream or lake is "dry" or that there is no water at the gauge.
  #    This symbol is used for water level data only.
  #R - Revised: The symbol R indicates that a revision, correction or addition
  # `  has been made to the historical discharge database after January 1, 1989.

  # ggplot(cancheck[Value > 0, ], aes(x=Date, y=Value, color=Symbol)) +
  #   geom_vline(data=cancheck[is.na(Value),], aes(xintercept = Date), color='grey', alpha=1/4) +
  #   geom_point(alpha=1/6) +
  #   geom_point(data=cancheck[Value==0,]) +
  #   facet_wrap(~STATION_NUMBER, scales='free') +
  #   theme_classic()

  #Remove 06NB002, 06AF001 — CA_0003544, CA_0003473

  #Check others 'CN_0000047', 'NO_0000018', 'RU_0000089',
  #'RU_0000391', 'RU_0000393', 'RU_00000395', 'RU_0000436',
  #'RU_0000470', 'US_0008687')
  check <- readformatGSIMmon(GSIMstatsdt[gsim_no == 'US_0008687',path])

  #Remove CN_
  GSIMtoremove_winterIR <- c('CA_0003544', #erroneous patterns (abnormally high values)
                             'CA_0003473', #sudden peak — unsure about estimated discharge under ice conditions
                             'CN_0000047', #Anomalous change from near 0 discharge to 150 m3/s, no explanation)
                             'NO_0000018', #Record for 18 years without a 0, 0s every month for last two years before discontinuation
                             'RU_0000391', #Stopped recording during the winter the last ~10 years. Maybe questionable winter data
                             'RU_0000393', #Didn't record during the winter for the first 20 years. Maybe questionable winter data
                             'RU_0000436', #Same
                             'RU_0000470', #Same
                             'US_0008687' #Just downstream of reservoir which appears to have caused intermittence
  )

  #-----  Check flags in coastal IR for GSIM
  GSImcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GSIMstatsdt,
                                     dbname = 'gsim', inp_resdir = inp_resdir,
                                     yearthresh = 1800, plotseries = plotseries)


  #Already checked CA_0006122
  GSIMtoremove_coastalIR <- c('NO_0000044',
                              'NO_0000090')


  #------ Remove stations with unstable intermittent flow regime
  #Remove those which have at least one day per year of zero-flow day but instances
  #of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
  GSIMtoremove_unstableIR <- GSIMstatsdt[
    (mDur_o1800 >= 1) & (!movinginter_o1800), gsim_no]

  ### Analyze GRDC data ####################################
  GRDCstatsdt <- rbindlist(in_GRDCgaugestats)

  #Remove all gauges with 0 values that have at least 99% of integer values as not reliable (see GRDC_NO 6140700 as example)
  GRDCtoremove_allinteger <- GRDCstatsdt[integerperc_o1800 >= 0.99 &
                                           intermittent_o1800 == 1, GRDC_NO]

  #Outliers from examining plots of ir time series (those that were commented out were initially considered)
  GRDCtoremove_o1800_irartifacts <- c(
    1134300, #remove - Appears to have been a change in flow permanence, but too much missing data
    1159830, #check missing values - probably perennial
    1160785, #remove - changed flow permanence
    1160800, #remove - changed flow permanence
    1160850, #remove - changed flow permanence
    1160881, #remove - changed flow permanence
    1196141, #remove - doesn't look reliable, hard to assess long term flow permanence
    1199410, #remove - changed flow permanence
    1259500, #remove - hard to assess long term flow permanence
    1428500, #remove - changed flow permanence
    1434200, #remove - almost all integers, impossible to tell flow permanence
    1434300, #remove - almost all integers
    1491815, #remove - changed flow permanence
    1591110, #remove - doesn't look reliable, changed flow permanence
    3652050, #remove - looks naturally perennial
    3652200, #remove - looks natural perennial (first 15 years)
    4208372, #remove
    4208655, #remove - insufficient data to tell flow permanence
    4208855, #remove - insufficient data to tell flow permanence
    4208857, #remove - probably not intermittent only occured once - missing years suggest perennial too
    4213566, #remove - probably perennial
    4214297, #remove - insufficient data to tell flow permanence
    4214298, #remove - probably perennial
    4769200, ##remove - insufficient data to tell flow permanence
    5204170, #remove - probably changed flow permanence
    5405095, ##remove - hard to tell whether actually intermittent (change at beginning, maybe regulated?)
    5708200 #remove - probably perennial
  )

  GRDCtoremove_o1961_irartifacts <- c(#1104800, #Keep
    1134500, #Goes from 1 cm3/s to 0. Occurs only one time in entire record.
    1159302, #Anomalous zero values
    1159303, #weird patterns, isolated 0s, sudden jumps and capped at 77
    #1159320, #0 values for the first 14 years but doesn't affect post-1961
    #1159325, #0 values for most record. probably episodic and due to series of agricultural ponds
    1159510, #0 values for most record, on same segment as 1159511 but weird record, remove.
    #1159512, #0 values for first 4 years. Otherwise no real reason to doubt.
    #1159520, #values seem capped after 1968, otherwise seem fine. Could just be rating curve
    1160101, #0 values are erroneous (goes from 1-10 m3/s to 0 then right back up )
    1160310, #Only at the beginning of record. Now regulated. Probably dried because of dam building or natural regime
    1160340, #pretty definitely erroneous. Missing data must have been labeled as 0.
    #1160378, #just downstream of dam. reservoir fully dry on satellite imagery, so keep as intermittent even when not regulated
    #1160420, #decreases to 0 appear a bit sudden but ok
    1160427, #0s at beginning of record, must be missing. Wouldn't be considered intermittent without that.
    1160435, #discard. 0s are mostly missing data. unreliable
    1160440, #only 0 values
    1160470, #Maybe discard. change of rating curve in 1947, mostly missing data until 1980 but truly intermittent based on imagery
    #1160510, #sudden shift in regime 1974. maybe change in bed morphology or dam building
    1160511, #unreliable, remove
    #1160540, #only 0 - nodata for first 15 years. seemingly good data post 1979. Keep.
    #1160580, #change of regime post-1974 from intermittent to perennial. not sure why. on same river as 1160510
    #1160635, #0 values post-1985 seem abrupt but enough values otherwise to make it intermittent
    1160650, #0 values are outliers
    1160670, #0 values are outliers due to dam just upstream
    #1160675, #Some outliers but otherwise most 0 values seem believable
    1160701, #0 values look like outliers
    1160756, #most 0 values look like outliers
    1160770, #only happened once. maybe outlier zeros
    1160780, #unreliable record. remove
    #1160785, #some 0 values look like outliers but most seem believable
    #1160790, #weird record pre 1963 but does not affect much. rating curve must have changed
    1160795, #most 0 values look like outliers
    1160825, #first 2 years of data are 0s, some other outlying 0s but enough believable ones to be intermittent
    1160840, #only 2 zero values are believable, others are outliers
    #1160850, #unsure
    1160880, #outlying 0 values. Tugela river. Perennial
    1160900, #most 0 values look like outliers
    1160911, #Remove. only one zero-flow vlaues appears right
    1160971, #check record. most zero-flow values look outliers.
    1160975, #Remove. unreliable
    1196100, #Half of the zero-flow values are outliers
    1196102, #Remove, unreliable
    #1196160, #Some outlying 0 values but most are good
    1197500, #Remove, unreliable
    1197540, #Check, abrupt 0 values. probably remove
    1197591, #Check, abrupt 0 values. probably remove
    1197700, #Check, abrupt 0 values. probably remove
    #1197740, #several outlying 0 flow values, but most are believable
    1199100, #Most 0s are outliers
    1199200, #Remove
    #1199410, #a few outliers pre-1990 but otherwise good
    1259800, #Remove, 0 values come from integer-based record
    1286690, #Remove, it seems that values were rounded to second decimal
    1289230, #Change in regime past 1963, only one outlier zero flow value after, remove
    1428400, #0 values come from integer-based record
    1434810, #0 values come from integer-based record
    #1491870, #some outliers but most 0s are believable
    1494100, #outlier 0s
    #1496351, #regulated but naturally intermittent
    1591720, #buggy, some values also look interpolated
    1591730, #anomalous drops
    1733600, #0 values come from integer-based record and outlier
    1837410, #Check, abrupt 0 values. probably remove
    #1837430, #Essentially the same as 1837410. Naturally intermittent before dam
    #1896501, #Faulty values but intermittent
    1897550, #outlier 0s
    1898501, #outlier 0s
    1992400, #nearly half of 0s are seemingly outliers
    1992600, #Outlier 0 values
    2549230, #Remove
    2588500, #Remove
    2588551, #Remove
    2588630, #Remove
    2588640, #Remove
    2588708, #Remove
    2588820, #Remove
    2589230, #Remove
    2589370, #Remove
    2591801, #Remove
    2694450, #Remove
    2917100, #Remove
    2969081, #Remove
    2999920, #0s come from integer-vased values
    3627900, #Remove
    #3650460, #some issues but enough believable intermittent ones
    3650470, #0s come from integer-based values
    #3650475, #earlier values pre-1973 are erroneous but post-2002 are valid 0s
    #3650610, #integers pre-1960s but still intermittent after
    3650640, #lots of outliers. remove
    3650649, #change of flow regime due to dam building but intermittent before
    3650690, #weird values but intermittent
    3650928, #weird. remove
    3652135, #some outliers, but enough 0s otherwise
    3844460, #Remove
    #4101451, #station downstream also has 0s
    4103700, #zero values are outliers
    #4113304, #Not super representative of current regime
    4122150, #just downstream of dam.
    4148955, #just downstream of dam which made it intermittent. May even be tidally influenced
    4149411, #Remove
    #4150605, #just downstream of lwesville dam in Dallas. previously intermittent as well but will be removed anyways as >50% dor
    #4151513, #Looks regulated but will be removed as > 50% regulated #################### good example for showing effect of regulation
    4152120, #regulated. intermittency occured right at the building of Alamo Lake on the Bill William
    4208195, #0s stem from interpolation
    #4213540, #one outlier, rest is good although occurred only once in 1965..
    #4213675, #looks funny but just downstream of natural lake
    4213905, #Regulated. hence the intermittency
    4213911, #Is regulated, remove
    4214075, ##Remove
    4234300, #remove. regulated
    4243610, #remove, regulated ############################ Good example of regulation
    4351710, #Remove. 0s due to outlier and integer-based values
    4355500, #Remove. outlier values. regulated
    4357510, #Remove
    4773050, #look erroneous (going from 140L/s to 0 in one day. must be a rating curve issue)
    5101020, #Remove
    5101101, #Remove
    5101130, #Remove
    5101201, #Remove
    #5101290, #valkues before 2000 aren't great. But 0 values post 2000 are believable
    5101305, #most 0 values are wrong
    5101380, #Remove
    5101381, #Remove
    5109110, #Remove
    5109200, #Remove, 0 values at beginning, interpolations used, not reliable
    5109230, #Remove
    5109251, #Remove
    5202140, #Remove
    5202145, #Most 0s look erroneous
    5202228, #Regulated? Remove
    5302229, #Remove
    5302251, #Remove
    5302261, #Remove

    #5405046, #canal through adelaide - sturt river
    #5405095, #looks wierd but intermittent it is
    5606130, #Remove
    5608100, #Remove
    5803160, #Remove
    5864500, #Remove. rounded to 10L/s
    5870100, #remove, rounded to 100L/S
    6119100, #Remove, rounded to 10L/s
    #6128220, #Weird. occurred only once but seems believable
    6140700, #Remove, perfect example of what an integer-based record involves
    6401800, #Remove, became regulated
    6442300, ##Remove, perfect example of what an integer-based record involves
    6444250, #Remove
    6444350, #Remove
    6444400, #Check, abrupt 0 values
    6935570 #Remove, rounded to 10L/s
  )

  #plotGRDCtimeseries(GRDCstatsdt[GRDC_NO == 4101450,])

  #### Check intermittent record
  # checkno <- 6444400 #GRDC_NO
  # check <- checkGRDCzeroes( #Check area around 0 values
  #   GRDCstatsdt, in_GRDC_NO=checkno, period=15, yearthresh=1800,
  #   maxgap=20, in_scales='free', labelvals = F)
  # checkno %in% GRDCtoremove_allinteger #Check whether all integers
  # in_gaugep[in_gaugep$GRDC_NO==checkno & !is.na(in_gaugep$GRDC_NO), "dor_pc_pva"] #check DOR
  # GRDCstatsdt[GRDC_NO == checkno, integerperc_o1800] #Check % integers

  #Outliers from examining plots of perennial time series (those that were commented out were initially considered)
  #Try to find those:
  # whose low flow plateaus could be 0s
  # whose perennial character is dam-driven or maybe irrigation driven (changed from IR to perennial but hard to find)
  # whose missing data are actually 0s
  # whose quality is too low to be reliable
  GRDCtoremove_o1800_pereartifacts <- c(
    1159800, #look regulated
    #1160302, #looked regulated -- but nothing obvious
    4101200, #low flows may be zeros as they plateau unless all integers
    4118850, #low flows plateau
    4125903, #too regulated
    4126351, #looks too regulated (as far as the record goes)
    4148850 #unsure, missing data have lots of zeros
  )

  GRDCtoremove_o1961_pereartifacts <- c(
    1160331, #remove- Lower plateaus are likely overestimated 0 values
    #1593100 #checked - clearly not intermittent but really bad quality
    1593751, #remove - missing values seem to contain intermittency
    1599111, #remove - maybe going down to 0 in missing data
    1160788, #remove - looks like plateauing at 0 but rating curve is off
    1899100, #remove-maybe intermittent now. intermittent last year of record (missing gap)
    3628200, #remove- lower values may be 0s but rating curve is off
    3652030, #remove- 0s in missing years and low flows in other years may also be 0s
    4115225, #remove -looks like it became regulated and may have otherwised become intermittent
    4146610, #remove- low flows may be 0s but rating curve is off
    4151801, #remove - regulated and may otherwise go dry Rio Grande
    4152651, #remove - regulated by blue mesa reservoir, may have been intermittent otherwise ######### good example
    4208610, #remove - too much missing data but if not would be intermittent ###################good example of that
    4213055, #remove- too much missing data but if not would be intermittent
    4213802, #remove - identical to 4213801
    4214320, #remove - lower values may be 0s and in missing years
    #4231620, #check -- maybe regulated but would probably otherwise be perennial
    4362100, #remove = lower values may be 0s and in missing years
    5606090, #remove - lower values may be 0s
    56064140, #remove - lower values may be 0s
    6123630, #remove- lower values may be 0s
    6233410, #remove -  looks erroneous
    6335020, #remove - looks identical to 633060
    6335050, #remove - looks identical to 6335060
    6337503, #remove -looks heavily regulated. cannot tell whether may have been intermittent before
    6442100, #remove - identical to 6442600
    6935146, #remove - looks identical to 6935145
    6972350, #remove
    6935600  #remove - identical to 6935145
  )

  #---------- Check flags in winter IR
  plot_winterir(dt = GRDCstatsdt, dbname = 'grdc', inp_resdir = inp_resdir,
                yearthresh = 1800, plotseries = plotseries)

  #Checked for seemingly anomalous 0s. Sudden decreases.
  #Check for flags, check satellite imagery, station name, check for construction of reservoir
  #If no way to explain, remove or if caused by reservoir/dam that is not in GranD
  #4220310 and 4243610, just downstream of dams that are in GranD — should be taken in account
  GRDCtoremove_winterIR <- c('2588640', #Sudden shift
                             '2589230',#Sudden shift
                             '4213540',#Sudden shift
                             '4214075',#Sudden shift
                             '6401800' #Just downstream of a reservoir that is not in GranD
  )

  #------ Check time series of stations within 3 km of seawater
  GRDCcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GRDCstatsdt,
                                     dbname = 'grdc', inp_resdir = inp_resdir,
                                     yearthresh = 1800, plotseries = plotseries)
  GRDCcoastalirall[, unique(readformatGRDC(path)$Flag), by=GRDC_NO]
  #Nothing obviously suspect beyond those that ad already been flagged

  #Inspect statistics for 4208857, 4213531 as no flow days occurred only one year
  # ID = '6976300'
  # GRDCstatsdt[GRDC_NO == ID,]
  # check <- readformatGRDC(GRDCstatsdt[GRDC_NO == ID,path])
  # unique(check$Flag)
  #
  # plotGRDCtimeseries(GRDCstatsdt[GRDC_NO == ID,], outpath=NULL)


  #------ Remove stations with unstable intermittent flow regime
  #Remove those which have at least one day per year of zero-flow day but instances
  #of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
  GRDCtoremove_unstableIR <- GRDCstatsdt[
    (mDur_o1800 >= 1) & (!movinginter_o1800) &
      !(GRDC_NO %in% c(1160115, 1160245, 4146400)), GRDC_NO]

  #Before cleaning
  GRDCtoremove_all <- unique(c(GRDCtoremove_allinteger,
                               GRDCtoremove_o1800_irartifacts,
                               GRDCtoremove_o1800_pereartifacts,
                               GRDCtoremove_o1961_irartifacts,
                               GRDCtoremove_o1961_pereartifacts,
                               GRDCtoremove_winterIR,
                               GRDCtoremove_unstableIR))

  GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10, .N]
  GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10 &
                !(GRDC_NO %in% GRDCtoremove_all), .N]

  ### Check changes in GSIM discharge data availability and flow regime over time ####
  GSIMstatsdt_clean <- GSIMstatsdt[!(gsim_no %in%  c(GSIMtoremove_o1961_irartifacts,
                                                     GSIMtoremove_o1800_irartifacts,
                                                     GSIMtoremove_coastalIR,
                                                     GSIMtoremove_winterIR,
                                                     GSIMtoremove_unstableIR)),]
  mvars <- c('intermittent_o1800',
             'intermittent_o1961',
             'intermittent_o1971')
  alluv_formatGSIM <- melt(GSIMstatsdt_clean,
                           id.vars = c('gsim_no',
                                       paste0('totalYears_kept_o',
                                              c(1800, 1961, 1971))),
                           measure.vars = mvars) %>%
    .[totalYears_kept_o1800 < 10 & variable %in% mvars, value := NA] %>%
    .[totalYears_kept_o1961 < 10 & variable %in% mvars[2:3], value := NA] %>%
    .[totalYears_kept_o1971 < 10 & variable %in% mvars[3], value := NA] %>%
    .[, count := .N, by=.(variable, value)]


  ### Check changes in GRDC discharge data availability and flow regime over time ####
  GRDCstatsdt_clean <- GRDCstatsdt[!(GRDC_NO %in% GRDCtoremove_all),]

  alluv_formatGRDC <- melt(GRDCstatsdt_clean,
                           id.vars = c('GRDC_NO',
                                       paste0('totalYears_kept_o',
                                              c(1800,1961, 1971))),
                           measure.vars = mvars) %>%
    .[totalYears_kept_o1800 < 10 & variable %in% mvars, value := NA] %>%
    .[totalYears_kept_o1961 < 10 & variable %in% mvars[2:3], value := NA] %>%
    .[totalYears_kept_o1971 < 10 & variable %in% mvars[3], value := NA] %>%
    .[, count := .N, by=.(variable, value)]


  ###Analyze change in number of gauges with different intermittency criterion
  irsensi_format <- melt(rbind(GRDCstatsdt_clean, GSIMstatsdt_clean,
                               use.names=TRUE, fill=T)[totalYears_kept_o1961 >= 10,],
                         id.vars = c('GRDC_NO', 'gsim_no'),
                         measure.vars = paste0('mDur_o', c(1800, 1961, 1971))) %>%
    .[!is.na(value) & value >0,] %>%
    setorder(variable, -value) %>%
    .[, cumcount := seq(.N), by=.(variable, is.na(GRDC_NO))]

  ggirsensi <- ggplot(irsensi_format, aes(x=value, y=cumcount,
                                          color=variable, linetype=is.na(GRDC_NO))) +
    geom_line(size=1.1) +
    coord_cartesian(expand=0, clip='off') +
    scale_x_sqrt(breaks=c(1, 5, 10, 30, 90, 180, 365),
                 labels=c(1, 5, 10, 30, 90, 180, 365)) +
    geom_vline(xintercept=c(1, 5)) +
    annotate(geom='text', x=c(1.7,6.5), y=150, angle=90,
             label=c(sum(irsensi_format[value==1 & variable=='mDur_o1800',
                                        max(cumcount), by=is.na(GRDC_NO)]$V1),
                     sum(irsensi_format[value==5 & variable=='mDur_o1800',
                                        max(cumcount), by=is.na(GRDC_NO)]$V1))) +
    theme_classic()

  plots <- grid.arrange(
    ggalluvium_gaugecount(dtformat = alluv_formatGRDC, alluvvar = 'GRDC_NO'),
    ggalluvium_gaugecount(dtformat = alluv_formatGSIM, alluvvar = 'gsim_no'),
    ggirsensi
  )

  ### Bind GRDC and GSIM records ####################################
  databound <- rbind(GRDCstatsdt_clean,
                     GSIMstatsdt_clean,
                     use.names=TRUE, fill=T)

  return(list(plots=plots, data=databound))
}

#------ format_gaugestats --------------------------------------------------------
#' Format gauge statistics
#'
#' Format gauge attributes
#'
#' @param path file path to a GRDC-formatted streamflow time series table
#' @param maxgap maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param monthsel selected months to compute the statistics over
#'
#' @return One row data.table with the following columns: \cr
#' \describe{
#' \item{GRDC_NO} - (char) unique identifier for the gauge
#' \item{firstYear} - (num) first year on full record
#' }
#'
#' @export

format_gaugestats <- function(in_gaugestats, in_gaugep, yearthresh) {
  #Join intermittency statistics to predictor variables and subset to only
  #include those gauges with at least 10 years of data
  gaugestats_join <- in_gaugestats[
    , GAUGE_NO := fifelse(is.na(gsim_no), GRDC_NO, gsim_no)] %>%
    .[!is.na(get(paste0('totalYears_kept_o', yearthresh))) &
        get(paste0('totalYears_kept_o', yearthresh))>=10,] %>%  # Only keep stations with at least 10 years of data pas yearthresh
    merge(as.data.table(in_gaugep)[, -c('GRDC_NO', 'gsim_no'), with=F],
          by='GAUGE_NO', all.x=T, all.y=F) %>%
    .[, c('X', 'Y') := as.data.table(sf::st_coordinates(geometry))] %>%
    .[, DApercdiff := (area_correct-UPLAND_SKM)/UPLAND_SKM] %>%
    .[, DApercdiffabs := abs(DApercdiff)]

  #Check for multiple stations on same HydroSHEDS segment
  dupliseg <- gaugestats_join[duplicated(HYRIV_ID) |
                                duplicated(HYRIV_ID, fromLast = T),] %>%
    .[, interdiff := length(unique(
      eval(paste0('intermittent_o', yearthresh))))-1, by=HYRIV_ID]

  #Only two stations have differing status. Inspect their record
  # dupliseg[interdiff==1,
  #          c('HYRIV_ID', 'GAUGE_NO', 'DApercdiff',
  #            'mDur_o1961', paste0('intermittent_o', yearthresh)), with=F]
  # plotGRDCtimeseries(dupliseg[GAUGE_NO == '1197560',])
  # plotGRDCtimeseries(dupliseg[GAUGE_NO == '1197591',]) #Remove because lots of erroneous data
  # plotGRDCtimeseries(dupliseg[GAUGE_NO == '4150605',])
  # plotGSIMtimeseries(dupliseg[GAUGE_NO == 'US_0006104',]) #Remove because identical record but without precise zero flow day count

  gaugestats_joinsel <- gaugestats_join[!GAUGE_NO %in% c('1197591', 'US_0006104')] %>%
    setorder(HYRIV_ID, DApercdiffabs) %>%
    .[!duplicated(HYRIV_ID),]

  print(paste0('Removing ', gaugestats_joinsel[dor_pc_pva >= 5000, .N],
               ' stations with >=50% flow regulation'))
  gaugestats_derivedvar <- gaugestats_joinsel[dor_pc_pva < 5000, ] %>% #Only keep stations that have less than 50% of their discharge regulated by reservoir
    comp_derivedvar #Compute derived variables, rescale some variables, remove -9999

  return(gaugestats_derivedvar)
}

#------ selectformat_predvars -----------------
selectformat_predvars <- function(inp_riveratlas_meta, in_gaugestats) {
  #---- List predictor variables ----
  monthlydischarge_preds <- paste0('DIS_',
                                   str_pad(seq(1, 12), 2, pad=0),
                                   '_CMS')

  predcols<- c(
    #monthlydischarge_preds,
    'UPLAND_SKM',
    'dis_m3_pyr',
    'dis_m3_pmn',
    'dis_m3_pmx',
    'dis_m3_pvar',
    'dis_m3_pvaryr',
    'run_mm_cyr',
    'runc_ix_cyr', #runoff coefficient (runoff/precipitation)
    'sdis_ms_uyr', #specific discharge
    'sdis_ms_umn',
    'inu_pc_umn',
    'inu_pc_umx',
    'inu_pc_cmn',
    'lka_pc_cse',
    'lka_pc_use',
    #'dor_pc_pva', #anthropogenic - degree of regulation
    'gwt_m_cav',
    'ele_pc_rel',
    'slp_dg_cav',
    'slp_dg_uav',
    'clz_cl_cmj',
    'snw_pc_uyr',
    'snw_pc_cyr',
    'snw_pc_cmx',
    'glc_cl_cmj',
    'glc_pc_c16',
    'glc_pc_u16',


    'pnv_cl_cmj',
    'wet_pc_cg1',
    'wet_pc_cg2',
    'wet_pc_ug1',
    'wet_pc_ug2',
    'wet_pc_c07',
    'wet_pc_c09',
    'wet_pc_u07',
    'wet_pc_u09',
    'for_pc_use',
    'for_pc_cse',
    #'ire_pc_use', #anthropogenic - irrigated area extent
    #'ire_pc_cse', #anthropogenic - irrigated area extent
    'gla_pc_use',
    'gla_pc_cse',
    'prm_pc_use',
    'prm_pc_cse',
    'swc_pc_uyr',
    'swc_pc_cyr',
    'swc_pc_cmn',
    'lit_cl_cmj',
    'kar_pc_use',
    'kar_pc_cse',
    #'ppd_pk_cav', #anthropogenic - pop density
    #'ppd_pk_uav', #anthropogenic - pop density
    #'urb_pc_cse', #anthropogenic - urban cover
    #'urb_pc_use', #anthropogenic - urban cover
    #'hft_ix_c93', #anthropogenic - human footprint
    #'hft_ix_u93', #anthropogenic - human footprint
    #'hft_ix_c09', #anthropogenic - human footprint
    #'hft_ix_u09', #anthropogenic - human footprint
    #'hdi_ix_cav', #anthropogenic - dev index

    'cly_pc_cav',
    'cly_pc_uav',
    'slt_pc_cav',
    'slt_pc_uav',
    'snd_pc_cav',
    'snd_pc_uav',
    # 'pre_mm_c01',
    # 'pre_mm_c02',
    # 'pre_mm_c03',
    # 'pre_mm_c04',
    # 'pre_mm_c05',
    # 'pre_mm_c06',
    # 'pre_mm_c07',
    # 'pre_mm_c08',
    # 'pre_mm_c09',
    # 'pre_mm_c10',
    # 'pre_mm_c11',
    # 'pre_mm_c12',
    # 'pre_mm_u01',
    # 'pre_mm_u02',
    # 'pre_mm_u03',
    # 'pre_mm_u04',
    # 'pre_mm_u05',
    # 'pre_mm_u06',
    # 'pre_mm_u07',
    # 'pre_mm_u08',
    # 'pre_mm_u09',
    # 'pre_mm_u10',
    # 'pre_mm_u11',
    # 'pre_mm_u12',
    'aet_mm_cyr',
    'aet_mm_uyr',
    'pet_mm_cyr',
    'pet_mm_uyr',
    # 'cmi_ix_cyr',
    # 'cmi_ix_uyr',
    'cmi_ix_cmn',
    'cmi_ix_umn',
    # 'cmi_ix_cvar',
    # 'cmi_ix_uvar',
    # 'cmi_ix_c01',
    # 'cmi_ix_c02',
    # 'cmi_ix_c03',
    # 'cmi_ix_c04',
    # 'cmi_ix_c05',
    # 'cmi_ix_c06',
    # 'cmi_ix_c07',
    # 'cmi_ix_c08',
    # 'cmi_ix_c09',
    # 'cmi_ix_c10',
    # 'cmi_ix_c11',
    # 'cmi_ix_c12',
    # 'cmi_ix_u01',
    # 'cmi_ix_u02',
    # 'cmi_ix_u03',
    # 'cmi_ix_u04',
    # 'cmi_ix_u05',
    # 'cmi_ix_u06',
    # 'cmi_ix_u07',
    # 'cmi_ix_u08',
    # 'cmi_ix_u09',
    # 'cmi_ix_u10',
    # 'cmi_ix_u11',
    # 'cmi_ix_u12',
    'ari_ix_cav',
    'ari_ix_uav',
    'wloss_pc_cav',
    'wdryp_pc_cav',
    'wwetp_pc_cav',
    'whfrq_pc_cav',
    'wseas_pc_cav',
    'wperm_pc_cav',
    'wfresh_pc_cav',
    'wloss_pc_uav',
    'wdryp_pc_uav',
    'wwetp_pc_uav',
    'whfrq_pc_uav',
    'wseas_pc_uav',
    'wperm_pc_uav',
    'wfresh_pc_uav',
    'bio1_dc_cav',
    'bio2_dc_cav',
    'bio3_dc_cav',
    'bio4_dc_cav',
    'bio5_dc_cav',
    'bio6_dc_cav',
    'bio7_dc_cav',
    'bio8_dc_cav',
    'bio9_dc_cav',
    'bio10_dc_cav',
    'bio11_dc_cav',
    'bio12_mm_cav',
    'bio13_mm_cav',
    'bio14_mm_cav',
    'bio15_mm_cav',
    'bio16_mm_cav',
    'bio17_mm_cav',
    'bio18_mm_cav',
    'bio19_mm_cav',
    'bio1_dc_uav',
    'bio2_dc_uav',
    'bio3_dc_uav',
    'bio4_dc_uav',
    'bio5_dc_uav',
    'bio6_dc_uav',
    'bio7_dc_uav',
    'bio8_dc_uav',
    'bio9_dc_uav',
    'bio10_dc_uav',
    'bio11_dc_uav',
    'bio12_mm_uav',
    'bio13_mm_uav',
    'bio14_mm_uav',
    'bio15_mm_uav',
    'bio16_mm_uav',
    'bio17_mm_uav',
    'bio18_mm_uav',
    'bio19_mm_uav'
  )

  #Check that all columns are in dt
  message(paste(length(predcols[!(predcols %in% names(in_gaugestats))]),
                'variables are missing from formatted gauge dataset'))

  #---- Associate HydroATLAS column names with variables names ----

  #Get predictor variable names
  metaall <- readxl::read_xlsx(inp_riveratlas_meta,
                               sheet='Overall') %>%
    setDT

  metascale <- readxl::read_xlsx(inp_riveratlas_meta,
                                 sheet='scale') %>%
    setDT %>%
    setnames(c('Key','Spatial representation'),
             c('Keyscale', 'Spatial.representation'))

  metastat <- readxl::read_xlsx(inp_riveratlas_meta,
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
                     fifelse(grepl("[0-9]",Keystat),
                             str_pad(Keystat, 2, side='left', pad='0'),
                             Keystat)),
    varname = paste(Attribute,
                    Spatial.representation,
                    Temporal.or.statistical.aggregation.or.other.association))]

  #Add newly generated variables to meta_format (variable labels)
  addedvars <- data.table(varname=c('Precipitation catchment Annual min/max',
                                    'Discharge watershed Annual min/max',
                                    'Discharge watershed Annual min/average',
                                    'Elevation catchment average - watershed average',
                                    'Runoff coefficient catchment Annual average',
                                    'Specific discharge watershed Annual average',
                                    'Specific discharge watershed Annual min',
                                    paste0('Discharge watershed ', month.name),
                                    'Drainage area',
                                    'Groundwater table depth catchment average'),
                          varcode=c('pre_mm_cvar',
                                    'dis_m3_pvar', 'dis_m3_pvaryr',
                                    'ele_pc_rel',
                                    'runc_ix_cyr',
                                    'sdis_ms_uyr', 'sdis_ms_umn',
                                    monthlydischarge_preds,
                                    'UPLAND_SKM',
                                    'gwt_m_cav'
                          )
  )

  oldcolnames <- c('Spatial.representation',
                   'Temporal.or.statistical.aggregation.or.other.association',
                   'Source Data')
  newcolnames <- c('Spatial representation',
                   'Temporal/Statistical aggreg.',
                   'Source')

  predcols_dt <- merge(data.table(varcode=predcols),
                       rbind(meta_format, addedvars, fill=T),
                       by='varcode', all.x=T, all.y=F)   %>%
    setnames(oldcolnames, newcolnames) %>%
    setorder(Category, Attribute,
             `Spatial representation`, `Temporal/Statistical aggreg.`)

  #Format table
  predcols_dt[varcode=='UPLAND_SKM', `:=`(
    Category = 'Physiography',
    Attribute= 'Drainage Area',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='',
    Source = 'HydroSHEDS',
    Citation = 'Lehner & Grill 2013'
  )]

  predcols_dt[varcode=='dis_m3_pvar', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='dis_m3_pvaryr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='runc_ix_cyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Runoff coefficient',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2, WorldClim v2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='sdis_ms_uyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='sdis_ms_umn', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  # predcols_dt[grepl('DIS_[0-9]{2}.*', varcode), `:=`(
  #   Category = 'Hydrology',
  #   Attribute= 'Natural Discharge',
  #   `Spatial representation`='p',
  #   `Temporal/Statistical aggreg.`= gsub('[A-Z_]', '', varcode),
  #   Source = 'WaterGAP v2.2',
  #   Citation = 'Döll et al. 2003'
  # )]

  predcols_dt[varcode=='ele_pc_rel', `:=`(
    Category = 'Physiography',
    Attribute= 'Elevation',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='(cav-uav)/uav',
    Source = 'EarthEnv-DEM90',
    Citation = 'Robinson et al. 2014'
  )]

  predcols_dt[varcode=='cmi_ix_cvar', `:=`(
    Category = 'Climate',
    Attribute= 'Climate Moisture Index',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WorldClim v2 & Global-PET v2',
    Citation = 'Fick et al. 2017',
    varname = 'Climate moisture index catchment monthly mn/mx'
  )]

  predcols_dt[varcode=='cmi_ix_uvar', `:=`(
    Category = 'Climate',
    Attribute= 'Climate Moisture Index',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WorldClim v2 & Global-PET v2',
    Citation = 'Fick et al. 2017',
    varname = 'Climate moisture index watershed monthly mn/mx'
  )]

  predcols_dt[varcode=='gwt_m_cav', `:=`(
    Category = 'Hydrology',
    Attribute= 'Groundwater table depth',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='av',
    Source = 'Global Groundwater Map',
    Citation = 'Fan et al. 2013'
  )]

  #Remove duplicates (which were created with keystat meaning different things e.g; 09 meaning september, 2009, class 9)
  predcols_dtnodupli<- predcols_dt[!(
    (Category == "Climate" &
       grepl('Class.*', `Temporal/Statistical aggreg.`)) |
      (Category == "Landcover" &
         !grepl('(Class|Spatial).*', `Temporal/Statistical aggreg.`))
  ),]  %>%
    .[grepl('hft_ix_[cu]09', varcode),
      `:=`(`Temporal/Statistical aggreg.`='2009',
           varname = gsub('(?<=Human\\sFootprint\\s(watershed|catchment)).*',
                          ' 2009',
                          varname,
                          perl=T)
      )] %>%
    .[grepl('hft_ix_[cu]93', varcode),
      `:=`(`Temporal/Statistical aggreg.`='1993',
           varname = gsub('(?<=Human\\sFootprint\\s(watershed|catchment)).*',
                          ' 1993',
                          varname,
                          perl=T)
      )] %>%
    unique(by='varcode')

  return(predcols_dtnodupli)
}

#------ create_tasks  -----------------
create_tasks <- function(in_gaugestats, in_predvars,
                         id_suffix=NULL, include_discharge = TRUE) {
  #Create subset of gauge data for analysis (in this case, remove records with missing soil data)
  datsel <- in_gaugestats[, c('intermittent_o1800',in_predvars$varcode, 'X', 'Y'),
                          with=F] %>%
    na.omit

  #Remove WaterGAP variables if include_discharge == TRUE
  if (!include_discharge) {
    watergapcols <- c('dis_m3_pyr',
                      'dis_m3_pmn',
                      'dis_m3_pmx',
                      'dis_m3_pvar',
                      'dis_m3_pvaryr',
                      'run_mm_cyr',
                      'runc_ix_cyr', #runoff coefficient (runoff/precipitation)
                      'sdis_ms_uyr', #specific discharge
                      'sdis_ms_umn')
    datsel <- datsel[, (watergapcols) := NULL]
  }

  #Basic task for classification
  task_classif <- mlr3spatiotempcv::TaskClassifST$new(
    id=  paste0("inter_class", id_suffix),
    backend = datsel,
    target = "intermittent_o1800",
    coordinate_names = c("X", "Y"))

  #Basic task for regression without oversampling
  task_regr <- convert_clastoregrtask(in_task = task_classif,
                                      in_id = paste0('inter_regr', id_suffix),
                                      oversample=FALSE)

  #Basic task for regression with oversampling to have the same number of minority and majority class
  task_regrover <- convert_clastoregrtask(in_task = task_classif,
                                          in_id = paste0('inter_regrover',
                                                         id_suffix),
                                          oversample=TRUE)
  return(list(classif=task_classif, regr=task_regr, regover=task_regrover))
}

#------ create_baselearners -----------------
create_baselearners <- function(in_task) {
  #---------- Create learners --------------------------------------------------
  lrns <- list()

  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }

  if (inherits(in_task, 'TaskClassif')) {
    #Compute ratio of intermittent to perennial observations
    imbalance_ratio <- get_oversamp_ratio(in_task)$ratio

    #Create basic learner
    lrns[['lrn_ranger']] <- mlr3::lrn('classif.ranger',
                                      num.trees = 800,
                                      sample.fraction = 0.632,
                                      replace = FALSE,
                                      splitrule = 'gini',
                                      predict_type = 'prob',
                                      importance = 'impurity_corrected',
                                      respect.unordered.factors = 'order')

    #print(lrn_ranger$param_set)

    #Create a conditional inference forest learner with default parameters
    #mtry = sqrt(nvar), fraction = 0.632
    lrns[['lrn_cforest']] <- mlr3::lrn('classif.cforest',
                                       ntree = 800,
                                       fraction = 0.632,
                                       replace = FALSE,
                                       alpha = 0.05,
                                       mtry = round(sqrt(length(in_task$feature_names))),
                                       predict_type = "prob")

    #Create mlr3 pipe operator to oversample minority class based on major/minor ratio
    #https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
    #https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
    #Sampling happens only during training phase.
    po_over <- mlr3pipelines::po("classbalancing", id = "oversample", adjust = "minor",
                                 reference = "minor", shuffle = TRUE,
                                 ratio = imbalance_ratio)
    #table(po_over$train(list(in_task))$output$truth()) #Make sure that oversampling worked

    #Create mlr3 pipe operator to put a higher class weight on minority class
    po_classweights <- mlr3pipelines::po("classweights", minor_weight = imbalance_ratio)

    #Create graph learners so that oversampling happens systematically upstream of all training
    lrns[['lrn_ranger_overp']] <- mlr3pipelines::GraphLearner$new(
      po_over %>>% lrns[['lrn_ranger']])
    lrns[['lrn_cforest_overp']] <- mlr3pipelines::GraphLearner$new(
      po_over %>>% lrns[['lrn_cforest']])

    #Create graph learners so that class weighin happens systematically upstream of all training
    lrns[['lrn_ranger_weight']] <- mlr3pipelines::GraphLearner$new(
      po_classweights %>>% lrns[['lrn_ranger']])
    lrns[['lrn_cforest_weight']] <- mlr3pipelines::GraphLearner$new(
      po_classweights  %>>% lrns[['lrn_cforest']])
  }


  if (inherits(in_task, 'TaskRegr')) {
    #Create regression learner with maxstat
    lrns[['lrn_ranger_maxstat']] <- mlr3::lrn('regr.ranger',
                                              num.trees=800,
                                              sample.fraction = 0.632,
                                              min.node.size = 10,
                                              replace=FALSE,
                                              splitrule = 'maxstat',
                                              importance = 'impurity_corrected',
                                              respect.unordered.factors = 'order')
  }

  return(lrns)
}

#------ set_tuning -----------------
set_tuning <- function(in_learner, in_measure, nfeatures,
                       insamp_nfolds, insamp_neval, insamp_nbatch) {

  if (is.list(in_learner)) {
    in_learner <- in_learner[[1]]
  }

  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)
    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T),
                   lower = floor(0.1*nfeatures),
                   upper = floor(0.5*nfeatures)), #Half number of features
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

  #Define termination rule
  evalsn = mlr3tuning::trm("evals", n_evals = insamp_neval) #termine tuning after insamp_neval rounds

  if (in_learner$task_type == 'classif') {
    if (inherits(in_measure, 'list')) {
      in_measure <- in_measure$classif
    }

    if (grepl('classif[.]cforest$', in_learner$id)) {
      learnertune <- in_learner
    } else if (grepl('classif[.]ranger$', in_learner$id)) {
      learnertune <- AutoTuner$new(learner= in_learner,
                                   resampling = rcv_rf,
                                   measure = in_measure,
                                   search_space = regex_tuneset(in_learner),
                                   terminator = evalsn,
                                   tuner =  tnr("random_search",
                                                batch_size = insamp_nbatch)) #batch_size determines level of parallelism
    } else{
      stop('The classification learner provided is not configurable with this workflow yet...')
    }
  } else if (in_learner$task_type == 'regr') {
    if (inherits(in_measure, 'list')) {
      in_measure <- in_measure$regr
    }

    learnertune <- AutoTuner$new(learner= in_learner,
                                 resampling = rcv_rf,
                                 measure = in_measure,
                                 search_space = regex_tuneset(in_learner),
                                 terminator = evalsn,
                                 tuner =  tnr("random_search",
                                              batch_size = insamp_nbatch))
  }

  learnertune$id <- in_learner$id

  return(learnertune)
}

#------ set_cvresampling ---------------
set_cvresampling <- function(rsmp_id, in_task, outsamp_nrep, outsamp_nfolds) {
  #repeated_cv or repeated-spcv-coords
  outer_resampling = rsmp(rsmp_id,
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)
  outer_resampling$instantiate(in_task)

  return(outer_resampling)
}


#------ dynamic_resample ------------------
#Run resample on in_task, selected learner (in_lrnid) from in_bm, in_resampling
dynamic_resample <- function(in_task, in_learner, in_resampling, type,
                             store_models = TRUE) {
  if (is.list(in_learner)) {
    in_learner <- in_learner[[1]]
  }

  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }

  if (inherits(in_learner, 'BenchmarkResult')) {
    print(('BenchmarkResults was provided, getting the learner...'))
    in_learner <- in_learner$learners$learner[[1]]
  }

  #Make sure autotuner matches task (adjust mtry)
  if (inherits(in_learner, 'AutoTuner')) {
    in_learner <- reset_tuning(in_autotuner = in_learner,
                               in_task = in_task)
  }

  if ((in_learner$task_type == 'classif' & type=='classif') |
      (in_learner$task_type == 'regr' & type=='regr')) {
    resmp_rs <- mlr3::resample(learner = in_learner, task = in_task,
                               resampling = in_resampling, store_models = store_models
    )
    return(resmp_rs)
  }
}

#------ dynamic_resamplebm ------------------
#Run resample on in_task, selected learner (in_lrnid) from in_bm, in_resampling
dynamic_resamplebm <- function(in_task, in_bm, in_lrnid, in_resampling, type,
                               store_models = TRUE) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  #get desired resampled_results/learner
  in_rf <- in_bm$filter(learner_ids = in_lrnid)

  rsmp_out <- dynamic_resample(in_task = in_task,
                               in_learner = in_rf,
                               in_resampling = in_resampling,
                               type = type,
                               store_models = store_models)

  return(rsmp_out)
}

#------ combined_bm: resample results into benchmark results -------------
combine_bm <- function(in_resampleresults, write_qs = NULL, inp_resdir = NULL) {
  #When tried as_benchmark_result.ResampleResult, got "Error in setcolorder(data, slots) :
  # x has some duplicated column name(s): uhash. Please remove or rename the
  # duplicate(s) and try again.". SO use this instead
  print('Converting to benchmark results...')
  if (length(in_resampleresults) > 1) {
    bmres_list <- lapply(
      in_resampleresults[!sapply(in_resampleresults, is.null)],
      function(rsmpres) {
        print(rsmpres)
        if (!is.null(rsmpres)) {
          as_benchmark_result(rsmpres)
        }
      })
    #BenchmarkResult$new(rsmpres$data)})

    print('Combining...')
    bmrbase = bmres_list[[1]]
    for (i in 2:length(bmres_list)) {
      if (in_resampleresults[[i]]$task$task_type ==
          in_resampleresults[[1]]$task$task_type) {
        print(i)
        bmrbase$combine(bmres_list[[i]])
      } else {
        warning('ResampleResult #', i,
                'is not of the same task type as the first ResampleResult you provided, skipping...')
      }
    }
  } else {
    warning('You provided only one resample result to combine_bm,
            simply returning output from as_benchmark_result...')
    bmrbase = as_benchmark_result(in_resampleresults[[1]])
  }
  print('Done combining, now writing to qs...')
  if (write_qs) {
    out_qs <- file.path(inp_resdir,
                        paste0('combine_bm',
                               format(Sys.time(), '%Y%m%d%H%M%s'),
                               '.qs')
    )
  }

  qs::qsave(bmrbase, out_qs)

  return(out_qs)
}

#------ select_features ------------------
select_features <- function(in_bm, in_lrnid, in_task, pcutoff) {

  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  #get desired resampled_results/learner
  in_rf <- in_bm$filter(learner_ids = in_lrnid)

  #Apply feature/variable selection
  vimp <- weighted_vimportance_nestedrf(
    rfresamp = in_rf$resample_result(uhash=unique(as.data.table(in_rf)$uhash)),
    pvalue = TRUE) %>%
    .[,imp_wmeanper := imp_wmean/sum(imp_wmean)]

  task_featsel <- in_task$clone()$select(
    vimp[imp_pvalue <= pcutoff, as.character(varnames)])
  task_featsel$id <- paste0(in_task$id, '_featsel')

  return(list(in_task, task_featsel))
}

#------ selecttrain_rf -----------------
selecttrain_rf <- function(in_rf, in_learnerid, in_taskid = NULL,
                           insamp_nfolds =  NULL, insamp_nevals = NULL) {
  #Prepare autotuner for full training
  if (inherits(in_rf, 'ResampleResult')) {
    in_bmsel <- in_rf$clone()
    lrn_autotuner <- in_bmsel$learners[[1]]
    in_task <- in_bmsel$task
    outer_resampling_output <- in_rf

  } else {
    in_bmsel <- in_rf$clone()$filter(learner_ids = in_learnerid,
                                     task_id = in_taskid)

    lrn_autotuner <- in_bmsel$clone()$learners$learner[[1]]
    in_task <-in_bmsel$tasks$task[[1]]

    #Return outer sampling object for selected model (or list of outer sampling objects)
    uhashes <- unique(as.data.table(in_bmsel)$uhash)
    if (length(uhashes) == 1) {
      outer_resampling_output <- in_bmsel$resample_result(uhash=uhashes)
    } else {
      outer_resampling_output <- lapply(uhashes, function(x) {
        in_bmsel$resample_result(uhash=x)
      })
    }
  }

  if (!is.null(insamp_nfolds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- insamp_nfolds
  }

  if (!is.null(insamp_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- insamp_nevals
  }

  #Train learners
  lrn_autotuner$param_set$values = mlr3misc::insert_named(
    lrn_autotuner$param_set$values,
    list(classif.ranger.importance = 'permutation')
  )
  lrn_autotuner$train(in_task)

  return(list(rf_outer = outer_resampling_output, #Resampling results
              rf_inner = lrn_autotuner, #Core learner (with hyperparameter tuning)
              task = in_task)) #Task
}


#------ rformat_network ------------------
rformat_network <- function(in_predvars, in_monthlydischarge=NULL,
                            inp_riveratlasmeta, inp_riveratlas, inp_riveratlas2) {
  cols_toread <-  unique(
    c("HYRIV_ID", "HYBAS_L12", "HYBAS_ID03", "LENGTH_KM",'INLAKEPERC',
      in_predvars[, varcode], 'pop_ct_csu', 'pop_ct_usu',
      'ele_mt_cav','ele_mt_uav', 'gwt_cm_cav', 'dor_pc_pva', 'ORD_STRA',
      #paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
      #paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
      paste0('pet_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
      paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0)))
  )

  riveratlas <- fread_cols(file_name=inp_riveratlas,
                           cols_tokeep = cols_toread)

  if (!is.null(in_monthlydischarge)) {
    riveratlas <- merge(riveratlas, as.data.table(in_monthlydischarge),
                        by.x='HYRIV_ID', by.y='REACH_ID')
  }

  setorder(riveratlas, HYRIV_ID)

  cols_v11 <-  names(fread(inp_riveratlas2, nrows=1)) %>%
    gsub('_11$', '', .)

  cols_tokeep <- names(riveratlas)[!(names(riveratlas) %in% cols_v11)]

  riveratlas_format <- fread(inp_riveratlas2) %>%
    setorder(REACH_ID) %>%
    setnames(gsub('_11$', '', names(.))) %>%
    cbind(riveratlas[, cols_tokeep, with=F], .) %>%
    comp_derivedvar

  return(riveratlas_format)
}

#------ make_gaugepreds -------------
make_gaugepreds <- function(in_rftuned, in_res_spcv, in_gaugestats,
                            in_predvars, interthresh) {
  #Get fully trained predictions
  gpreds_full <- in_rftuned$rf_inner$predict(in_rftuned$task)

  #Get unique ID for each spatial fold-repetition combination
  spcv_clusters <- in_res_spcv$resampling$instance
  spcv_clustersformat <- dcast(spcv_clusters, row_id~rep, value.var = 'fold') %>%
    setnames(old=unique(spcv_clusters$rep)+1, paste0('fold', unique(spcv_clusters$rep))) %>%
    .[, spfold_unique := do.call(paste0, .SD),
      .SDcols = paste0('fold', unique(spcv_clusters$rep))] %>%
    setorder(row_id)

  #Get average predictions for oversampled rows and across repetitions - simple CV
  rsmp_res_nosp <- get_outerrsmp(in_rftuned, spatial_rsp=FALSE)
  gpreds_CV_nosp <-  rsmp_res_nosp$prediction() %>%
    as.data.table %>%
    .[, list(truth=first(truth), prob.1=mean(prob.1)), by=row_id]%>%
    setorder(row_id)

  #Get average predictions for oversampled rows and across repetitions - spatial CV
  # rsmp_res_sp <- get_outerrsmp(in_rftuned, spatial_rsp=TRUE)
  # gpreds_CV_sp <-  rsmp_res_sp$prediction()$set_threshold(1-interthresh) %>%
  #   as.data.table %>%
  #   .[, list(truth=first(truth), prob.1=mean(prob.1)), by=row_id] %>%
  #   .[, response := fifelse(prob.1 >= interthresh, 1, 0)] %>%
  #   setorder(row_id)

  #Format gauge data.table prior to merging
  datsel <- na.omit(in_gaugestats, c('intermittent_o1800',
                                     in_predvars$varcode, 'X', 'Y'))

  #Merge all predictions
  datsel[, `:=`(
    IRpredprob_full = gpreds_full$prob[,2],
    # IRpredprob_CVsp = gpreds_CV_sp$prob.1,
    # IRpredcat_CVsp = gpreds_CV_sp$response,
    IRpredprob_CVnosp = gpreds_CV_nosp$prob.1,
    spfold_unique =  spcv_clustersformat$spfold_unique
  )] %>%
    setnames(old='spfold_unique', new=paste0(in_rftuned$task$id, '_spfold'))

  return(datsel)
}

#------ bind_gaugepreds -----------------
bind_gaugepreds <- function(in_gpredsdt, interthresh) {
  #Get binary classification threshold
  if (inherits(interthresh, 'data.table')) {
    interthresh <- interthresh[learner == in_learnerid, thresh]
  }

  if (is.list(in_gpredsdt) & !inherits(in_gpredsdt, 'data.table')) {
    in_gpredsdt <- rbindlist(in_gpredsdt,
                            use.names = TRUE, fill = TRUE, idcol = "modelgroup")
  }

  if (in_gpredsdt[duplicated(GAUGE_NO), .N] > 0 ) {
    out_gpredsdt <- in_gpredsdt[, `:=`(
      IRpredprob_full=mean(IRpredprob_full, na.rm=T),
      IRpredprob_CVnosp=mean(IRpredprob_CVnosp, na.rm=T)
    ),
    by=.(GAUGE_NO)] %>% #Compute mean predicted probability if overlapping models
      .[!duplicated(GAUGE_NO),]     #Remove duplicates (if overlapping model)
  } else {
    out_gpredsdt <- in_gpredsdt
  }

  out_gpredsdt[,`:=`(IRpredcat_full = fifelse(IRpredprob_full >= interthresh,
                                              '1', '0'),
                     IRpredcat_CVnosp =fifelse(IRpredprob_CVnosp >= interthresh,
                                              '1', '0'),
                     preduncert_full = IRpredprob_full -
                       as.numeric(as.character(intermittent_o1800)),
                     preduncert_CVnosp = IRpredprob_CVnosp -
                       as.numeric(as.character(intermittent_o1800))
  )]

  return(out_gpredsdt)
}



#------ write_gaugepreds -----------------
write_gaugepreds <- function(in_gaugep, in_gpredsdt, outp_gaugep) {
  cols_toditch<- colnames(in_gpredsdt)[
    !(colnames(in_gpredsdt) %in% c('GAUGE_NO', 'geometry'))]

  out_gaugep <- base::merge(
    in_gaugep[,- which(names(in_gaugep) %in% cols_toditch)],
    in_gpredsdt[, -'geometry', with=F],
    by='GAUGE_NO',
    all.x=F)

  st_write(obj=out_gaugep,
           dsn=outp_gaugep,
           driver = 'gpkg',
           delete_dsn=T)

  return(out_gaugep)
}

#------ write_netpreds -------------------
write_netpreds <- function(in_network, in_rftuned, in_predvars,
                           discharge_interval = c(-Inf, Inf),
                           interthresh = 0.5, outp_riveratlaspred) {

  write_netpred_util <- function(in_network, in_rftuned, in_predvars,
                                 discharge_interval, interthresh) {
    # ----- Make predictions across river network -----
    #Get subset of river network
    in_network <- in_network[(dis_m3_pyr >= discharge_interval[1]) &
                               (dis_m3_pyr < discharge_interval[2]),]

    #Get rows for which no predictor variable is NA (seecomp_derivedvar for formatting/determining variables)
    netnoNArows <- in_network[, c(!(.I %in% unique(unlist(
      lapply(.SD, function(x) which(is.na(x))))))),
      .SDcols = in_predvars$varcode]

    #Predict model — chunk it up by climate zone to avoid memory errors
    for (clz in unique(in_network$clz_cl_cmj)) {
      print(clz)
      tic()
      in_network[netnoNArows & clz_cl_cmj == clz,
                 predbasic800 := in_rftuned$rf_inner$predict_newdata(
                   as.data.frame(.SD))$prob[,2],]
      toc()
    }

    #Label each reach categorically based on threshold
    in_network[, predbasic800cat := fifelse(predbasic800>=interthresh, 1, 0)]

    return(in_network[, c('HYRIV_ID', 'HYBAS_L12',
                          'predbasic800', 'predbasic800cat'), with=F])
  }


  networkpreds <- mapply(FUN = write_netpred_util,
                         in_rftuned = in_rftuned,
                         discharge_interval = discharge_interval,
                         interthresh = interthresh,
                         MoreArgs = list(
                           in_network = in_network,
                           in_predvars = in_predvars
                         ),
                         USE.NAMES = TRUE,
                         SIMPLIFY = FALSE
  ) %>%
    rbindlist

  if (networkpreds[duplicated(HYRIV_ID), .N] > 0 ) {
    networkpreds <- networkpreds[, predbasic800 := mean(predbasic800, na.rm=T),
                                 by=.(HYRIV_ID)] %>% #Compute mean predicted probability if overlapping models
      .[!duplicated(HYRIV_ID),] %>%     #Remove duplicates (if overlapping model)
      .[, predbasic800cat := fifelse(predbasic800 >= interthresh, '1', '0')]
  }

  outp_riveratlaspred <- paste0(gsub('[.][a-z]*$', '', outp_riveratlaspred), '_',
                                format(Sys.Date(), '%Y%m%d'), '.csv')
  fwrite(networkpreds, outp_riveratlaspred)

  # --------- Return data for plotting ------------------------
  return(outp_riveratlaspred)
}

# list(out_gaugep = out_gaugep,
#      rivpredpath =


##### -------------------- Diagnostics functions -------------------------------

#------ netpredformat ------
netpredformat <- function(outp_riveratlaspred, in_rivernetwork) {
  fread(file_in(outp_riveratlaspred))%>%
    .[in_rivernetwork[, c('HYRIV_ID', 'HYBAS_L12', 'LENGTH_KM', 'dis_m3_pyr',
                          'UPLAND_SKM', 'clz_cl_cmj'),
                      with=F], on='HYRIV_ID']
}

#------ ggmisclass_single -----------------

ggmisclass_single <- function(in_predictions=NULL, in_rftuned=NULL, spatial_rsp=FALSE) {
  #Get predicted probabilities of intermittency for each gauge
  # in_gaugestats[!is.na(cly_pc_cav), intermittent_predprob :=
  #                 as.data.table(in_predictions)[order(row_id), mean(prob.1), by=row_id]$V1]
  #Get misclassification error, sensitivity, and specificity for different classification thresholds
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class
  if (!is.null(in_rftuned)) {
    rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
    in_predictions <- rsmp_res$prediction()
  }

  threshold_confu_dt <- ldply(seq(0,1,0.01), threshold_misclass, in_predictions) %>%
    setDT

  #Get classification threshold at which sensitivity and specificity are the most similar
  balanced_thresh <- threshold_confu_dt[which.min(abs(spec-sens)),]
  print(paste('Sensitivity =', round(balanced_thresh$sens,2),
              'and Specificity =', round(balanced_thresh$spec,2),
              'at a classification threshold of', balanced_thresh$i))

  gout <- ggplot(melt(threshold_confu_dt, id.vars='i'),
                 aes(x=i, y=value, color=variable, linetype=variable)) +
    geom_line(size=1.2) +
    geom_vline(xintercept=balanced_thresh$i, alpha=1/2) +
    geom_hline(yintercept=balanced_thresh$spec, alpha=1/2) +
    annotate('text', x=(balanced_thresh$i), y=0.4,
             label=balanced_thresh$i, angle=-90) +
    annotate('text', x=0.9, y=(balanced_thresh$spec),
             label=round(balanced_thresh$sens,2)) +
    scale_x_continuous(expand=c(0,0), name='Threshold') +
    scale_y_continuous(expand=c(0,0), name='Value') +
    scale_color_brewer(palette='Dark2',  #colorblind friendly
                       labels=c('Misclassification rate',
                                'Sensitivity (true positives)',
                                'Specificity (true negatives)')) +
    theme_bw()

  #Plot it
  return(list(plot = gout,
              interthresh = balanced_thresh$i))
}

#------ analyze_benchmark -----------------
analyze_benchmark <- function(in_bm, in_measure) {

  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  print(paste('It took',
              in_bm$aggregate(mlr3::msr('time_both'))$time_both,
              'seconds to train and predict with the',
              in_bm$aggregate(msr('time_both'))$learner_id,
              'model...'))

  bmdt <- as.data.table(in_bm)

  if (in_bm$task_type == 'regr') {
    print(in_bm$aggregate(in_measure$regr))
    boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure$regr)

    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- bmdt$prediction[[rsmp_i]] %>%
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
    print(in_bm$aggregate(measures=in_measure$classif))
    boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure$classif)

    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- data.table(outf = bmdt$iteration[[rsmp_i]],
                          task = bmdt$task[[rsmp_i]]$id,
                          learner = bmdt$learner[[rsmp_i]]$id,
                          task_type = in_bm$task_type,
                          pred = list(bmdt$prediction[[rsmp_i]]))
      return(preds)
    }) %>%
      do.call(rbind, .)
  }

  tasklearner_unique <- preds[, expand.grid(unique(task), unique(learner))] %>%
    setnames(c('task', 'learner')) %>%
    setDT

  tasklearner_unique[, learner_format := dplyr::case_when(
    learner == 'classif.ranger'~'default RF',
    learner == 'oversample.classif.ranger'~'default RF - oversampled',
    learner == 'classweights.classif.ranger'~'default RF - weighted classes',
    learner == 'classif.cforest'~'CIF',
    learner == 'oversample.classif.cforest'~'CIF - oversampled',
    learner == 'classweights.classif.cforest'~'CIF - weighted classes',
  )]

  glist <- lapply(1:nrow(tasklearner_unique), function(tsklrn) {
    print(tasklearner_unique[tsklrn,])
    subpred <- preds[task ==tasklearner_unique$task[tsklrn] &
                       learner == tasklearner_unique$learner[tsklrn],]

    ggmisclass_out <- ggmisclass_single(in_predictions = subpred)

    gout <- ggmisclass_out$plot +
      ggtitle(paste(tasklearner_unique$task[tsklrn],
                    tasklearner_unique$learner_format[tsklrn])) +
      labs(x='Threshold', y='Value')

    if (tsklrn < nrow(tasklearner_unique)) {
      gout <- gout +
        theme(legend.position = 'none')
    }

    return(list(plot = ggplotGrob(gout),
                interthres_dt = data.table(
                  learner = as.character(tasklearner_unique$learner[tsklrn]),
                  thresh = ggmisclass_out$interthresh
                )
    )
    )
  }) %>%
    unlist(recursive=F)


  return(list(
    bm_misclasscomp=do.call("grid.arrange",
                            list(grobs=glist[seq(1, length(glist), 2)])), #Get all plots out of the nested list
    bm_boxcomp = boxcomp,
    interthresh_dt = rbindlist(glist[seq(2, length(glist), 2)]) #Get all threshold data.table rows out of nested list
  ))
}

#------ bin_misclass --------------------------------
bin_misclass <-  function(in_predictions=NULL, in_rftuned=NULL,
                          in_gaugestats, binvar, binfunc, binarg,
                          interthresh=0.5, spatial_rsp=FALSE) {
  #Get misclassification error, sensitivity, and specificity for different classification thresholds
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class
  if (!is.null(in_rftuned)) {
    if (inherits(in_rftuned, 'list')) {
      rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
    }
    else if (inherits(in_rftuned, 'ResampleResult')) {
      rsmp_res <- in_rftuned
    }

    in_predictions <- rsmp_res$prediction() %>%
      as.data.table

    #Get average predicted probability and response across CV repeats
    predsformat <- in_predictions[, list(prob.1 = mean(prob.1)),
                                  by=.(row_id, truth)] %>%
      .[, response := fifelse(prob.1 >= interthresh, 1, 0)]

    #Get bins
    bin_gauges <- bin_dt(in_dt = in_gaugestats, binvar = 'dis_m3_pyr',
                         binfunc = 'manual', binarg=binarg,
                         bintrans=NULL, na.rm=FALSE) %>%
      .[, row_id := .I]

  } else if (is.null(in_rftuned) & !is.null(in_predictions)) {
    #Get bins
    bin_gauges <- bin_dt(in_dt = in_predictions, binvar = 'dis_m3_pyr',
                         binfunc = 'manual', binarg=binarg,
                         bintrans=NULL, na.rm=FALSE) %>%
      .[, row_id := .I]

    predsformat <- in_predictions[, `:=`(truth = intermittent_o1800,
                                         response = IRpredcat_CVnosp,
                                         row_id = .I)]
  }

  #Merge bins with predictions
  gselpreds <- merge(bin_gauges, predsformat, by='row_id', all.y=F)

  #Function to get statistics for each bin
  binconfustats <- function(in_gselpreds) {
    confu <- in_gselpreds[, .N, by=.(truth, response)]

    outvec <- data.table(
      inter_confu =  paste0(confu[truth==1 & response==1, max(0L, N)], '|',
                            confu[truth==1 & response==0, max(0L, N)]),
      pere_confu = paste0(confu[truth==0 & response==1, max(0L, N)], '|',
                          confu[truth==0 & response==0, max(0L, N)]),
      misclas = round(confu[truth != response, sum(N)] / confu[, sum(N)], 2),
      sens = round(
        confu[truth == '1' & response == '1', max(0L, N)]/confu[truth=='1', sum(N)],
        2),
      spec  = round(
        confu[truth=='0' & response==0, max(0L, N)]/confu[truth=='0', sum(N)],
        2),
      N = sum(confu$N),
      predtrue_inter = paste0(round(confu[response==1, 100*sum(N)/sum(confu$N)]),
                              '|',
                              round(confu[truth==1, 100*sum(N)/sum(confu$N)]))
    )
    return(outvec)
  }

  #Run function for each bin
  performtable <- gselpreds[, binconfustats(.SD),
                            by=.(bin, paste0(round(bin_lmin, 2), '-',
                                             round(bin_lmax, 2)))] %>%
    setorder(bin) %>%
    rbind(binconfustats(gselpreds)[, paste0 := 'All'], use.names=T, fill=T)

  return(performtable)
}

#------ ggvimp -----------------
ggvimp <- function(in_rftuned, in_predvars, varnum = 10, spatial_rsp=FALSE) {
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)

  #Get variable importance and format them
  varimp_basic <- weighted_vimportance_nestedrf(rfresamp = rsmp_res,
                                                pvalue = FALSE) %>%
    merge(., in_predvars, by.x='varnames', by.y='varcode') %>%
    .[, `:=`(varname = factor(varname, varname[order(-imp_wmean)]),
             Category = factor(Category,
                               levels = c('Climate', 'Hydrology', 'Landcover',
                                 'Physiography', 'Soils & Geology'))
             )] %>%
    setorder(-imp_wmean)

  #Plot 'em
  outp <- ggplot(varimp_basic[1:varnum,],aes(x=varname, fill=Category)) +
    geom_bar(aes(y=imp_wmean), stat = 'identity') +
    geom_errorbar(aes(ymin=imp_wmean-imp_wsd, ymax=imp_wmean+imp_wsd)) +
    scale_x_discrete(labels = function(x) {
      stringr::str_wrap(tolower(x), width = 30)
    }) +
    scale_fill_manual(values=c('#fdb462','#80b1d3','#b3de69','#bc80bd','#d9d9d9'),
                      drop=FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(size=8),
          axis.title.x  = element_blank(),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.position = c(0.80, 0.5)
    ) +
    scale_y_continuous(expand=c(0,0)) +
    coord_flip(ylim=c(0, min(varimp_basic[, max(imp_wmean+imp_wsd)+1], 100)),
               clip='off')

  return(outp)
}


#------ ggpd_bivariate -----------------
ggpartialdep <- function (in_rftuned, in_predvars, colnums, ngrid, nodupli=T,
                          nvariate = 2,
                          parallel=T, spatial_rsp=FALSE) {

  #Get outer resampling of interest
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)

  #Get partial dependence across all folds and repeats
  nlearners <-with(rsmp_res$resampling$param_set$values, folds*repeats)
  datdf <- as.data.frame(rsmp_res$task$data()) #This may be shortened
  varimp <- weighted_vimportance_nestedrf(rsmp_res, pvalue=FALSE) %>%
    setorder(-imp_wmean)

  if (length(colnums) > nrow(varimp)) {
    colnums <- colnums[1:nrow(varimp)]
    print('colnums argument exceeded the number of variables,
          reduced it to ', nrow(varimp), ' variables')
  }

  if (nodupli) {
    selcols <- as.character(
      varimp$varnames[!duplicated(substr(varimp$varnames, 1,3))][colnums])
  } else {
    selcols <- as.character(
      varimp$varnames[colnums])
  }

  if (parallel) {
    print(paste("Computing partial dependence with future.apply across", nlearners,
                "CV folds"))
    pd <- future.apply::future_lapply(seq_len(nlearners),
                                      extract_pd_nestedrf,
                                      in_rftuned = rsmp_res,
                                      datdf = datdf,
                                      selcols = selcols,
                                      nvariate = nvariate,
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
                 nvariate = nvariate,
                 ngrid = ngrid)
  }

  #Get weighted mean
  varvec <- paste0('var', 1:nvariate)
  valvec <- paste0('value', 1:nvariate)

  pdformat <- do.call(rbind, pd) %>%
    setDT %>%
    .[, list(mean1 = weighted.mean(`1`, classif.bacc)),
      by= c(varvec, valvec)] %>%
    .[, variables := var1] %>%
    merge(in_predvars[, .(varcode, varname)], by.x='var1', by.y='varcode')

  datdf2 <- as.data.table(datdf)[, intermittent_o1800 := as.numeric(as.character(intermittent_o1800))]

  if (nvariate ==1) {
    tileplots_l <- pdformat[,list(list(ggplotGrob(
      ggplot(.SD, aes(x=value1, y=mean1)) +
        geom_line() +
        geom_rug(data=datdf2,
                 aes_string(x=eval(var1),y='intermittent_o1800'),
                 alpha=1/3) +
        scale_y_continuous(name='Partial dependence (probability of intermittency)',
                           limits= c(min(mean1)-0.01, max(mean1)+0.01),  #c(0.25, 0.425),
                           expand=c(0,0))+
        scale_x_continuous(name=stringr::str_wrap(eval(varname), width = 30)) +
        theme_classic() +
        theme(text = element_text(size=12),
              axis.title.y = element_blank())
    ))), by=.(var1)]

  } else if (nvariate == 2) {
    pdformat[, variables := paste(var1, var2)]

    vargrid <- t(combn(1:length(selcols), 2))
    #leglims <- pdformat[, c(min(mean1), max(mean1))]

    #Iterate over every pair of variables

    tileplots_l <- pdformat[,list(list(ggplotGrob(
      ggplot(.SD, aes(x=value1, y=value2)) +
        geom_tile(aes(fill = mean1)) +
        scale_fill_distiller(palette='Grey') +
        geom_jitter(data=datdf,
                    aes_string(color='intermittent_o1800', x=eval(var1),y=eval(var2)),
                    alpha=1/3) +
        scale_color_manual(values=c('#0F9FD6','#ff9b52')) +
        labs(x=stringr::str_wrap(in_predvars[varcode==eval(var1), varname],
                                 width = 20),
             y=stringr::str_wrap(in_predvars[varcode==eval(var2), varname],
                                 width = 20)) +
        theme_bw() +
        theme(text = element_text(size=12))
    )))
    , by=.(var1, var2)]
  }

  pagelayout <-   lapply(1:(nrow(tileplots_l) %/% 9), function(p_i) {
    (p_i-1)*9+(1:9)
  })
  if (nrow(tileplots_l) %% 9 > 0) {
    pagelayout[[nrow(tileplots_l) %/% 9 + 1]] <-
      (nrow(tileplots_l) %/% 9)*9+(1:(nrow(tileplots_l) %% 9))
  }


  tileplots_multipl <- lapply(pagelayout, function(page) {
    print(page)
    return(do.call("grid.arrange",list(
      grobs=(tileplots_l[page,V1]),
      left = 'Partial dependence (probability of intermittency)')))
  })
  return(tileplots_multipl)
}

#------ gggaugeIPR -----------------
gggaugeIPR <- function(in_gpredsdt, in_predvars, spatial_rsp = FALSE,
                       interthresh = 0.5, yearthresh, in_learnerid = NULL) {
  #Plot numeric variables
  in_gpredsdt[, DApercsdiff := abs(
    (area_correct-UPLAND_SKM)/mean(c(area_correct, UPLAND_SKM))), by=GAUGE_NO]

  predmelt_num <- in_gpredsdt[
    , which(as.vector(unlist(lapply(in_gpredsdt, is.numeric)))), with=F] %>%
    cbind(in_gpredsdt[, c('GAUGE_NO', 'intermittent_o1800'), with=F]) %>%
    melt(id.vars=c('GAUGE_NO', 'intermittent_o1800',
                   'IRpredprob_CVnosp', 'preduncert_CVnosp'))

  #Set variable labels
  varlabels <- copy(in_predvars) %>%
    .[, .(varcode, varname)] %>%
    setkey(varcode) %>%
    .[levels(predmelt_num$variable)] %>%
    .[is.na(varname), varname:=varcode] %>%
    .[varcode=='totalYears_kept_o1800', varname := 'Years of record kept'] %>%
    .[varcode=='yearskeptratio', varname := 'Years kept/Years total'] %>%
    .[varcode=='mDur_o1800', varname :='Mean annual # of dry days'] %>%
    .[varcode=='mFreq_o1800', varname :='Mean annual # of dry periods'] %>%
    .[varcode=='DApercsdiff', varname :='Absolute symmetric difference in reported drainage area (%)'] %>%
    .[varcode=='ORD_STRA', varname :='Strahler river order'] %>%
    .[varcode=='UPLAND_SKM', varname :='Drainage area (km2)'] %>%
    .[varcode=='dor_pc_pva', varname :='Degree of regulation'] %>%
    .[, varname := str_wrap(varname, 30)]

  levels(predmelt_num$variable) <-  varlabels$varname
  varstoplot_continuous <- varlabels[varcode %in% c('totalYears_kept_o1800',
                                                    'mDur_o1800',
                                                    'mFreq_o1800',
                                                    'DApercsdiff',
                                                    'ari_ix_uav'),]

  varstoplot_log <- varlabels[varcode %in% c('UPLAND_SKM',
                                             'dis_m3_pyr',
                                             'dor_pc_pva'),]

  colorpal <- c('#1f78b4', '#ff7f00')
  rectdf <- data.table(
    xmin=rep(0, 4),
    xmax=rep(Inf, 4),
    ymin=c(-1, interthresh-1, 0, interthresh),
    ymax=c(interthresh-1, 0, interthresh, 1),
    fillpal = rep(colorpal, 2)
  )

  ggIPR_util <- function(plotdt, logx = FALSE, legend=FALSE) {
    p <- ggplot(plotdt) +
      geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
                                 ymin=ymin, ymax=ymax, fill=fillpal),
                alpha=1/4) +
      scale_y_continuous(limits=c(-1, 1), expand=c(0,0)) +
      scale_fill_manual(values=colorpal,
                        name='Predicted regime',
                        labels = c('Perennial', 'Intermittent')) +
      geom_point(aes(x=value, y=preduncert_CVnosp, color=intermittent_o1800),
                 alpha = 1/8) +
      geom_hline(yintercept=0, alpha=1/2) +
      new_scale_fill() +
      geom_smooth(aes(x=value, y=preduncert_CVnosp, color=intermittent_o1800),
                  method='gam', formula = y ~ s(x, k=3)) +
      # annotate("text", x = Inf-5, y = 0.5, angle = 90,
      #          label = "Pred:Int, Obs:Per",
      #          color = colorpal[[1]]) +
      # annotate("text", x = Inf-5, y = -0.5, angle = 90,
      #          label = "Pred:Per, Obs:Int",
      #          color = colorpal[[2]]) +
      scale_color_manual(values=colorpal,
                         name='Observed regime',
                         labels = c('Perennial', 'Intermittent')) +
      labs(x='Value', y='Intermittency Prediction Residuals (IPR)') +
      #scale_x_sqrt(expand=c(0,0)) +
      coord_cartesian(expand=FALSE, clip='off') +
      facet_wrap(~variable, scales='free', ncol=3,
                 labeller=as_labeller(varlabels)) +
      theme_classic() +
      theme(legend.position = c(0.85, 0.2))

    if (logx) {
      p <- p + scale_x_continuous(
        trans='log1p',
        breaks = c(1, 10, 100, 1000, 10000, 100000,1000000),
        labels = scientific_10)
    }

    if (!legend) {
      p <- p + theme(legend.position = 'none')
    }

    return(p)
  }


  p_continuous <- ggIPR_util(
    plotdt = predmelt_num[variable %in% varstoplot_continuous$varcode |
                            variable %in% varstoplot_continuous$varname,],
    logx=F, legend=T) +
    labs(x="")
  p_log <- ggIPR_util(
    plotdt = predmelt_num[variable %in% varstoplot_log$varcode |
                            variable %in% varstoplot_log$varname,],
    logx=T, legend=F) +
    labs(y="", x='')


  #Plot categorical variables
  predmelt_cat <- in_gpredsdt[, c('GAUGE_NO', 'intermittent_o1800', 'preduncert_CVnosp',
                                  'ENDORHEIC', 'clz_cl_cmj'), with=F] %>%
    .[, `:=`(intermittent_o1800 = factor(intermittent_o1800, levels=c('0','1')))] %>%
    melt(id.vars=c('GAUGE_NO', 'intermittent_o1800', 'preduncert_CVnosp'))
  levels(predmelt_cat$variable) <- c('Endorheic',
                                     'Climate Zone (catchment majority)')

  rectdf_cat <- copy(rectdf)[,xmin:=-Inf]
  gaugeIPR_catplot <-
    ggplot(predmelt_cat) +
    geom_rect(data=rectdf_cat, aes(xmin=xmin, xmax=xmax,
                               ymin=ymin, ymax=ymax, fill=fillpal),
              alpha=1/4) +
    scale_fill_manual(values=colorpal,
                      name='Predicted regime',
                      labels = c('Perennial', 'Intermittent')) +
    new_scale_fill() +
    geom_boxplot(aes(x=as.factor(value), y=preduncert_CVnosp,
                     fill=intermittent_o1800, color=intermittent_o1800),
                 alpha=0.75) +
    facet_wrap(~variable, scales='free', labeller=label_value) +
    scale_fill_manual(values=colorpal,
                      name='Observed regime',
                      labels = c('Perennial', 'Intermittent')) +
    scale_color_manual(values=c('#175885', '#9e3f00'),
                       name='Observed regime',
                       labels = c('Perennial', 'Intermittent')) +
    geom_hline(yintercept=0, alpha=1/2) +
    labs(x='Value', y='') + #'Intermittency Prediction Residuals (IPR)') +
    coord_cartesian(expand=FALSE, clip='off') +
    theme_classic() +
    theme(legend.position='none')


  outp <- grid.arrange(p_continuous, p_log, gaugeIPR_catplot,
               layout_matrix = rbind(c(1,1,1),
                                     c(1,1,1),
                                     c(2,2,2),
                                     c(3,3,3))
  )


  return(outp)
}

#------ krige_spgaugeIPR----
krige_spgaugeIPR <- function(in_rftuned, in_gaugep, in_gaugestats,
                             kcutoff=50000, inp_bufrasdir,
                             overwrite = FALSE) {
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=TRUE)
  predsp <- rsmp_res$prediction() %>%
    as.data.table %>%
    .[, list(IRpredprob_spcv = 100*mean(prob.1),
             IRpredcoefcv_spcv = sd(100*prob.1, na.rm=T)/mean(100*prob.1, na.rm=T)),
      by=.(row_id, truth)] %>%
    .[, IPR := IRpredprob_spcv-100*as.numeric(as.character(truth))] %>%
    .[, IPR_abs := abs(IPR)] %>%
    .[, IRpedcat_spcv := as.character(fifelse(IRpredprob_spcv > 40, 1, 0))] %>%
    setorder(row_id) %>%
    cbind(in_gaugestats[!is.na(cly_pc_cav), list(GAUGE_NO=GAUGE_NO)])

  predsp_gaugep <- merge(in_gaugep, predsp, by='GAUGE_NO') %>%
    .[order(.[,'IPR_abs']$IPR_abs),]

  bufras_vec <- file.path(inp_bufrasdir,
                          grep('bufras_T.*proj[.]tif$',
                               list.files(inp_bufrasdir),
                               value=T)
  )

  #Convert to SpatialPointDataFrame to use in gstats and project to Goode Homolosine
  crs_aeqd <- "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  predsp_gaugep_df <- st_transform(predsp_gaugep, crs=crs_aeqd)  %>%
    as_Spatial()

  # Compute the sample variogram; note that the f.1 trend model is one of the
  # parameters passed to variogram(). This tells the function to create the
  # variogram on the de-trended data.
  var_smpl <- gstat::variogram(IPR ~ 1, predsp_gaugep_df,
                               cutoff=kcutoff, width=500, cressie=TRUE) #cloud=T)

  ggplot(var_smpl, aes(x=dist, y=gamma)) +
    geom_point(alpha=1/3) +
    coord_cartesian(expand=F) +
    scale_y_log10() +
    theme_classic()

  # scale_x_sqrt(breaks=c(0,100,1000,10000,25000,50000, 100000),
  #              labels=c(0,100,1000,10000,25000,50000, 100000)) +

  # Compute the variogram model by passing the nugget, sill and range values
  # to fit.variogram() via the vgm() function.
  varmods <- as.character(vgm()[,'short'])
  dat_fit  <- fit.variogram(var_smpl, fit.ranges = TRUE, fit.sills = FALSE,
                            vgm(varmods[!(varmods %in% c('Pow', 'Int'))]),
                            fit.kappa = TRUE)
  # The following plot allows us to assess the fit
  plot(var_smpl, dat_fit)

  # Import gauge buffer mask


  #Predict error
  kmod <- gstat(formula=IPR~1, locations=predsp_gaugep_df, model=dat_fit,
                maxdist = kcutoff)
  kpl <-  lapply(seq_along(bufras_vec), function(i) {
    bufmask <- raster(bufras_vec[[i]]) %>%
      as('SpatialGrid')
    kp <- round(raster(predict(kmod, bufmask)))

    outras = file.path(dirname(bufras_vec[[i]]),
                       gsub('bufras', 'krigpred', basename(bufras_vec[[i]])))


    print(paste0('Writing ', outras, '...'))
    writeRaster(kp, outras, datatype = "INT2S", overwrite=overwrite)
    return(outras)
  }) %>%
    do.call(rbind, .)

  return(kpl)
}
#------ mosaic_kriging -------------
mosaic_kriging <- function(in_kpathlist, outp_krigingtif, overwrite) {
  kpl <- in_kpathlist

  kplextents <- lapply(kpl, function(ras_path) {
    ras <- raster(ras_path)
    data.table(x=c(xmin(ras),xmax(ras)), y= c(ymin(ras),ymax(ras)))
  }) %>%
    rbindlist %>%
    .[, lapply(.SD, function(i) c(min(i), max(i)))]

  template <- raster(extent(kplextents[, x], kplextents[, y]))

  writeRaster(template, file=outp_krigingtif, datatype='INT2S',
              format="GTiff", overwrite=overwrite)
  gdalUtils::mosaic_rasters(gdalfile=as.vector(kpl),
                            dst_dataset=out_krigingtif,
                            force_ot = "Int16",
                            of="GTiff")
  gdalinfo(out_krigingtif)

  return(out_krigingtif)
}


#------ map_BACC ------
map_basinBACC <- function(in_gaugepred, #rfpreds_gauges,
                          in_rivernetwork, #################################### NEXT RUN - SHOULDN'T NEED THAT - COULD ADD HYBAS_ID03 DIRECTLY TO GAUGEP
                          inp_basin,
                          outp_basinerror) {
  bas03 <- st_read(dsn = dirname(inp_basin),
                   layer = basename(inp_basin))

  gnetjoin <- merge(in_gaugepred,
                    in_rivernetwork[,.(HYRIV_ID, HYBAS_ID03)],
                    by='HYRIV_ID', all.x=T, all.y=F)

  classcol='intermittent_o1800'
  dat <- in_gaugepred[, .N, by=classcol]%>%
    setnames(old=classcol, new='classcol')
  minoratio <- dat %>%
    setorder(N) %>%
    .[, list(minoclass=classcol[1], ratio=N[2]/N[1])]

  gnetjoin[, classweights := fifelse(intermittent_o1800 == minoratio$minoclass,
                                     minoratio$ratio, 1)]

  basbacc <- gnetjoin[, list(
    basbacc =  sum((intermittent_o1800==IRpredcat_CVnosp) *
                     classweights) / sum(classweights),
    acc = sum(intermittent_o1800==IRpredcat_CVnosp)/.N,
    gnum = .N,
    percinter = round((100*sum(intermittent_o1800=='1')/.N), 2),
    predbias = (sum(IRpredcat_CVnosp==1)/.N) - (sum(intermittent_o1800=='1')/.N)
  ), by=HYBAS_ID03]


  basbacc_format <- base::merge(
    bas03, basbacc,
    by.x='HYBAS_ID', by.y='HYBAS_ID03',
    all.x.=T, all.y=T)
  basbacc_format$gdens <- with(basbacc_format, gnum/UP_AREA)

  st_write(obj=basbacc_format,
           dsn=outp_basinerror,
           driver = 'gpkg',
           delete_dsn=F)
}

##### -------------------- Report functions -----------------------------------
#------ get_basemaps ------------
get_basemapswintri <- function() {
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"

  wcountries <- rnaturalearth::ne_countries(
    scale = "medium", returnclass = "sf") %>%
    sfformat_wintri
  wland <- rnaturalearth::ne_download(
    scale = 110, type = 'land', category = 'physical') %>%
    sfformat_wintri
  wlakes <- rnaturalearth::ne_download(
    scale = 110, type = 'lakes', category = 'physical') %>%
    sfformat_wintri

  grat_wintri <-
    st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
    st_transform(crs = crs_wintri)

  # riv_query <- "SELECT HYRIV_ID, ORD_STRA, dis_m3_pyr FROM \"RiverATLAS_v10\" WHERE ORD_STRA > 4"
  # riv_lines <- st_read(dsn = dirname(in_filestructure['in_riveratlas']),
  #                      layer = basename(in_filestructure['in_riveratlas']),
  #                      query = riv_query)
  #
  # riv_wintri <- sfformat_wintri(riv_lines)
  # riv_simple <-  sf::st_simplify(riv_wintri, preserveTopology = TRUE, dTolerance = 1000)

  return(list(wcountries = wcountries,
              wland = wland,
              wlakes = wlakes,
              grat = grat_wintri))
  #,riv_simple = riv_simple
}

#------ ggrivers -----------------------
ggrivers <- function(in_basemaps) {
  p <- ggplot() +
    geom_sf(data = in_basemaps[['grat']],
            color = alpha("black", 1/5), size = 0.25/.pt) +
    geom_sf(data = in_basemaps[['wland']],
            color='#dee0e0') +
    geom_sf(data = in_basemaps[['wlakes']],
            alpha=1/3, fill='black', color=NA) +
    # geom_sf(data=in_basemaps['riv_simple'][riv_simple$ORD_STRA < 6,][1:10000,],
    #         aes(size=ORD_STRA), color=alpha('black', 1/20), show.legend = F) +
    # geom_sf(data=in_basemaps['riv_simple'][riv_simple$ORD_STRA >= 6,][1:10000,],
    #         aes(size=ORD_STRA), color=alpha('black', 1/10), show.legend = F) +
    scale_size_continuous(range=c(0.2, 1.0)) +
    coord_sf(datum = NA, expand=F,
             xlim=c(-14500000, 18000000), ylim=c(-7000000,8700000)) +
    theme_map() +
    theme(legend.position=c(0.1,0),
          legend.direction="horizontal",
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.key.height = unit(0.1,"in"))

  return(p)
}

#------ gggauges --------------------------
gggauges <- function(in_gaugepred, in_basemaps,
                     binarg, binvar) {
  gaugepred <- in_gaugepred %>%
    sfformat_wintri

  #Bin data
  dt <- as.data.table(gaugepred)
  byd <- dt[, binvar, with=F]
  l_freq    <- append(min(byd), binarg)
  u_freq    <- c(binarg, (max(byd, na.rm = TRUE) + 1))
  bins      <- length(binarg) + 1

  for (i in seq_len(bins)) {
    gaugepred[dt[,get(binvar)] >= l_freq[[i]] &
                dt[,get(binvar)] < u_freq[[i]], 'bin'] <- i

    if (i == bins) {
      gaugepred[dt[,get(binvar)]  == u_freq[i], 'bin'] <- i
    }
  }

  gaugepred$bin <- factor(gaugepred$bin, levels=seq_len(bins))


  #Create color scales
  cs_per <- scales::seq_gradient_pal('#45c6f7', '#152540', "Lab")(
    seq(0,1,length.out=bins)) #'#0F9FD6'
  cs_ir <- scales::seq_gradient_pal('#ffb88f','#940404', "Lab")( #'#ff9b52'  '#a32d18'
    seq(0,1,length.out=bins)) #


  #Subset data into perennial and intermittent rivers
  perennial_gauges <- gaugepred[gaugepred$intermittent_o1800=='0',]
  ir_gauges <-  gaugepred[gaugepred$intermittent_o1800=='1',]


  #Make histograms
  gggaugehist <- function(in_gdf, cs) {
    x_tick <- c(0, unique(as.numeric(levels(in_gdf$bin)))) + 0.5
    binarg_tick <- c(0, seq_len(bins))
    len <- length(x_tick)

    permean <- as.data.table(in_gdf)[,  mean(get(binvar))]
    permeanbindiff <- c(permean - l_freq)
    meanbin <- which(permeanbindiff > 0)[which.min(permeanbindiff[permeanbindiff  > 0])]
    meanpos <- meanbin + (permean-l_freq[meanbin])/(u_freq[meanbin]-l_freq[meanbin])

    gaugehist <- ggplot(in_gdf, aes(x=as.numeric(bin), fill=bin)) +
      geom_histogram(stat="count", width=1) +
      geom_vline(xintercept = meanpos - 1) +
      annotate(geom='text', angle=90,
               x=meanpos-0.6,
               y=0.5*max(as.data.table(in_gdf)[, .N, by=bin]$N),
               label=paste(round(permean), 'y')) +
      scale_fill_manual(values=cs) +
      scale_x_continuous(name = 'Years of data',
                         breaks = c(binarg_tick, x_tick),
                         labels = c(rep(c(""), len), c(l_freq, u_freq[bins]))) +
      # scale_y_continuous(name = paste0('Number of gauging stations (total: ',
      #                                  nrow(in_gdf),
      #                                  ')'))+
      coord_cartesian(clip='off', expand=c(0,0)) +
      theme_classic() +
      theme(legend.position = 'none',
            text = element_text(size=12),
            plot.background = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_text(vjust=3),
            axis.ticks.x = element_line(color = c(rep(NA, len - 1),
                                                  rep("black", len))))
    return(gaugehist)
  }

  histper <- gggaugehist(perennial_gauges, cs=cs_per)
  histir <- gggaugehist(ir_gauges, cs=cs_ir)


  #Plot it
  p_pr <- ggrivers(in_basemaps) +
    geom_sf(data=perennial_gauges, aes(color=bin), size=1.2, alpha=0.8) +
    scale_color_manual(values=cs_per) +
    coord_sf(datum = NA, expand=F,
             xlim=c(-16000000, 18000000), ylim=c(-7000000,8700000)) +
    theme(legend.position = 'none') +
    annotation_custom(grob = ggplotGrob(histper),
                      xmin = -16000000,
                      xmax = -8000000,
                      ymin = -9000000,
                      ymax = -1000000)


  p_ir <- ggrivers(in_basemaps) +
    geom_sf(data=ir_gauges, aes(color=bin), size=1) +
    scale_color_manual(values=cs_ir) +
    coord_sf(datum = NA, expand=F,
             xlim=c(-16000000, 18000000), ylim=c(-7000000,8700000)) +
    theme(legend.position = 'none') +
    annotation_custom(grob = ggplotGrob(histir),
                      xmin = -16000000,
                      xmax = -8000000,
                      ymin = -9000000,
                      ymax = -1000000)

  p_patch <- p_pr/p_ir
  outp <- p_patch + plot_annotation(tag_levels = 'A')
  return(outp)
}

#------ get_gaugestats ----------
comp_GRDCqstats <- function(path, maxgap,
                            minyear = 1971 , maxyear = 2000, verbose = FALSE) {
  if (verbose) {
    print(path)
  }

  #Read and format discharge records
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  gaugetab <- readformatGRDC(path)

  #Remove years with too many gaps, no data records, and
  gaugetabsub <- gaugetab[(missingdays <= maxgap) &
                            (year >= minyear & year <= maxyear) &
                            !(Original %in% c(-999, -99, -9999, 99, 999, 9999) |
                                is.na(Original)),]

  gaugestatsyr <- gaugetabsub[, list(
    qminyr = min(Original, na.rm=T),
    q3minyr = min(frollmean(Original, 3, align='center'), na.rm=T)
  ),by=year] %>%
    .[, list(qminyr = mean(qminyr, na.rm=T),
             q3minyr = mean(q3minyr, na.rm=T))]

  gaugestats <- gaugetabsub[, list(
    GRDC_NO = unique(GRDC_NO),
    nyears = length(unique(year)),
    qmean = mean(Original, na.rm=T),
    q10 = quantile(Original, 0.9, na.rm=T),
    q90 = quantile(Original, 0.1, na.rm=T),
    q99 = quantile(Original, 0.01, na.rm=T)
  )] %>%
    cbind(gaugestatsyr)

  return(gaugestats)
}


#------ eval_watergap -------
eval_watergap <- function(in_qstats, in_selgauges, binarg) {
  qsub <- in_qstats[!is.na(GRDC_NO) &
                      nyears >= 20 &
                      GRDC_NO %in% in_selgauges$GRDC_NO,] %>%
    .[!(GRDC_NO %in% c('5204010')),]
  #5204010 - issue with units which does not affect zero flow assessment
  #Checked 1160520 - mixed but not clear enough

  qsubp <- merge(qsub, in_selgauges, by='GRDC_NO', all.x=T, all.y=F) %>%
    .[dor_pc_pva < 1000,]
  qsubp_bin <- bin_dt(in_dt = qsubp, binvar='qmean', binarg = binarg, binfunc = 'manual')

  #Compute a simple set of performance statistics, including rsquare for ols without studentized outliers
  getqstats <- function(dt, x, y, rstudthresh= 3, log=FALSE) {
    if (log) {
      in_form = paste0('log10(', y, '+0.1)~log10(', x, '+0.1)')
    } else {
      in_form = paste0(y,'~',x)
    }

    mod <- lm(as.formula(in_form), data=dt)

    outrows <- setDT(olsrr::ols_prep_rstudlev_data(mod)$`levrstud`)[
      abs(rstudent)>rstudthresh & color == 'outlier & leverage',]
    if (nrow(outrows) >0) {
      dtsub <- dt[-(outrows$obs),]
    } else {
      dtsub <- dt
    }

    mod_nooutliers <- lm(as.formula(in_form), data=dtsub)

    outstats <- dt[, list(pearsonr = round(cor(get(y), get(x)), 3),
              mae = round(Metrics::mae(get(y), get(x)), 2),
              smape = round(Metrics::smape(get(y), get(x)), 2),
              #pbias = round(Metrics::percent_bias(get(y), get(x))),
              rsq = round(summary(mod)$r.squared, 3),
              rsq_nooutliers = round(summary(mod_nooutliers)$r.squared, 3),
              n_total = .N,
              noutliers = nrow(outrows)
              )
       ,]

    return(outstats)
  }

  qsubp_bin[qmean >= 100,
            getqstats(dt=.SD, x='dis_m3_pyr', y='qmean',
                      rstudthresh = 3)]


  #Get stats for mean Q ~ dis_m3_pyr (watergap mean annual)
  qmean_stats <- qsubp_bin[,
                           getqstats(dt=.SD, x='dis_m3_pyr', y='qmean',
                                     rstudthresh = 3),
                           by=.(bin, bin_lformat)] %>%
    rbind(
      getqstats(dt=qsubp_bin, x='dis_m3_pyr', y='qmean', log=T)[
        , `:=`(bin_lformat='all', bin=length(binarg)+1)]
      ) %>%
    .[, comp := 'qmean_dism3pyr']

  #Get stats for Q90 ~ dis_m3_pyr (watergap min monthly)
  q90_stats <- qsubp_bin[,
                           getqstats(dt=.SD, x='dis_m3_pmn', y='q90',
                                     rstudthresh = 3),
                           by=.(bin, bin_lformat)] %>%
    rbind(getqstats(dt=qsubp_bin, x='dis_m3_pmn', y='q90', log=T)[
      , `:=`(bin_lformat='all', bin=length(binarg)+1)]) %>%
    .[, comp := 'q90_dism3mn']

  outstats <- rbind(qmean_stats, q90_stats) %>%
    setorder(-comp, bin)
  return(outstats)
}

#------ formatscales ------------
formatscales <- function(in_df, varstoplot) {
  scales_x <- list(
    ari_ix_uav = scale_x_continuous(expand=c(0,0)),
    bio12_mm_uav  = scale_x_sqrt(expand=c(0,0),
                                 labels=c(0, 1000, 2000, 5000, 10000)),
    bio14_mm_uav  = scale_x_sqrt(expand=c(0,0),
                                 breaks = c(0, 50, 100, 200, 500),
                                 labels=c(0, 50, 100, 200, 500)),
    cly_pc_uav = scale_x_continuous(labels=percent_format(scale=1), expand=c(0,0)),
    cmi_ix_uyr = scale_x_continuous(),
    dis_m3_pyr = scale_x_log10(breaks=c(1, 10^2,
                                        10^(0:log10(max(in_df$dis_m3_pyr)))),
                               labels=c(0, 10^2,
                                        10^(0:log10(max(in_df$dis_m3_pyr)))),
                               expand=c(0,0)),
    dor_pc_pva = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    for_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    gla_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    kar_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    lka_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    pet_mm_uyr = scale_x_continuous(expand=c(0,0)),
    sdis_ms_uyr = scale_x_continuous(expand=c(0,0)),
    snw_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    run_mm_cyr = scale_x_continuous(expand=c(0,0)),
    swc_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    tmp_dc_uyr = scale_x_continuous(expand=c(0,0)),
    hdi_ix_cav = scale_x_continuous(expand=c(0,0)),
    hft_ix_c93 = scale_x_continuous(expand=c(0,0)),
    ORD_STRA = scale_x_continuous(expand=c(0,0)),
    UPLAND_SKM = scale_x_log10(breaks=c(1, 10^2,
                                        10^(0:log10(max(in_df$UPLAND_SKM)))),
                               labels=c(1, 10^2,
                                        10^(0:log10(max(in_df$UPLAND_SKM)))),
                               expand=c(0,0)),
    gwt_m_cav = scale_x_sqrt(expand=c(0,0)),
    ire_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0))
  ) %>%
    .[(names(.) %in% names(in_df)) & names(.) %in% varstoplot]
  #Only keep those variables that are actually in df and that we want to plot

  scales_y <- unlist(rep(list(scale_y_continuous(expand=c(0,0))),
                         labels = scientific_format(),
                         length(scales_x)),
                     recursive=F) %>%
    setNames(names(scales_x))

  scales_y[['dis_m3_pmn']] <- scale_y_sqrt(expand=c(0,0))
  scales_y[['glc_pc_u16']] <- scale_y_continuous(trans='log1p',
                                                 breaks=c(10, 1000, 100000, 10000000))

  coordcart <- lapply(varstoplot, function(var) {
    coord_cartesian(xlim=as.data.table(in_df)[, c(min(get(var), na.rm=T),
                                                  max(get(var), na.rm=T))])
  }) %>%
    setNames(varstoplot)

  coordcart[['clz_cl_cmj']] <-  coord_cartesian(
    xlim=c(1,max(in_df$clz_cl_cmj)))
  coordcart[['kar_pc_use']] <-  coord_cartesian(
    xlim=c(0, 100))
  coordcart[['pet_mm_uyr']] <-  coord_cartesian(
    xlim=c(0, max(in_df$pet_mm_uyr)))
  coordcart[['ORD_STRA']] <-  coord_cartesian(
    xlim=c(1, 10))
  coordcart[['ari_ix_uav']] <-  coord_cartesian(
    xlim=c(0, 100))

  return(list(scales_x=scales_x, scales_y=scales_y, coordcart=coordcart))
}

#------ ggenvhist -------------
ggenvhist <- function(vartoplot, in_gaugedt, in_rivdt, in_predvars,
                      scalesenvhist, intermittent=TRUE) {
  print(vartoplot)
  if (intermittent) {
    vartoplot2 <- c(vartoplot, 'intermittent_o1800')
  }

  varname <- in_predvars[varcode==vartoplot, Attribute]
  #paste0(Attribute, ' ',Keyscale,Keystat,' (',unit,')')]

  if (vartoplot == "clz_cl_cmj") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=as.factor(clz_cl_cmj)]
    gclz <- in_gaugedt[,.N/in_gaugedt[,.N],by=as.factor(clz_cl_cmj)]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))

    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'dodge', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#2b8cbe', '#dd3497'))

  } else if (vartoplot == "glc_pc_u16") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=glc_pc_u16]
    gclz <- in_gaugedt[,.N/in_gaugedt[,.N],by=glc_pc_u16]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))

    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'identity', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#2b8cbe', '#dd3497'))
    #
    #     penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
    #       geom_histogram(data=in_rivdt, aes(weight = LENGTH_KM),
    #                      fill='#2b8cbe', alpha=0.5, bins=101) +
    #       geom_histogram(fill='#dd3497', alpha=0.5, bins=101)

  } else {
    penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
      geom_density(data=in_rivdt, aes(weight = LENGTH_KM),
                   fill='#2b8cbe', alpha=0.5) +
      geom_density(fill='#dd3497', alpha=0.5) +
      ylab('Density')
  }

  penvhist <- penvhist +
    scalesenvhist$scales_x[[vartoplot]] +
    #scalesenvhist$scales_y[[vartoplot]] +
    scalesenvhist$coordcart[[vartoplot]] +
    xlab(varname) +
    theme_classic() +
    theme(strip.background=element_rect(colour="white", fill='lightgray'),
          legend.position = 'none',
          axis.title.y = element_blank(),
          axis.title = element_text(size=12))

  # if (which(vartoplot %in% varstoplot_hist)!=length(varstoplot_hist)) {
  #   penvhist <- penvhist +
  #     theme(legend.position='none')
  # }

  return(ggplotGrob(penvhist))
}

#------ layout_ggenvhist --------------------------
layout_ggenvhist <- function(in_rivernetwork, in_gaugepred, in_predvars) {
  varstoplot_hist <- c(
    "bio1_dc_uav", "bio7_dc_uav", "bio12_mm_uav", "bio14_mm_uav", "clz_cl_cmj",
    "ari_ix_uav", "dis_m3_pyr", "sdis_ms_uyr", "gwt_m_cav", "UPLAND_SKM",
    "lka_pc_use", "snw_pc_uyr", "kar_pc_use", "for_pc_use") #, "glc_pc_u16")

  if ("dis_m3_pyr" %in% varstoplot_hist) {
    setDT(in_rivernetwork)[, dis_m3_pyr := dis_m3_pyr + 1]
    setDT(in_gaugepred)[, dis_m3_pyr := dis_m3_pyr + 1]
  }

  #Get legend
  pleg <- ggplot(in_gaugepred, aes(x=dis_m3_pyr, fill=factor(IRpredcat_full))) +
    geom_density(alpha=1/2) +
    scale_fill_manual(values=c('#2b8cbe', '#dd3497'),
                      name = 'Dataset',
                      labels=c('Global river network',
                               'Training gauges')) +
    theme(text=element_text(size=14))

  tmp <- ggplot_gtable(ggplot_build(pleg))
  leg <- tmp$grobs[[
    which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
  ]]

  #Get scales
  scalesenvhist <- formatscales(in_df=in_rivernetwork, varstoplot=varstoplot_hist)

  #Plot each facet
  penvhist_grobs <- lapply(varstoplot_hist, ggenvhist,
                           in_gaugedt = in_gaugepred,
                           in_rivdt = in_rivernetwork,
                           in_predvars = in_predvars,
                           scalesenvhist = scalesenvhist)
  #Add legend
  penvhist_grobs[[length(penvhist_grobs) + 1]] <- leg

  #Plot
  grid.newpage()
  do.call("grid.arrange", list(grobs=penvhist_grobs, nrow=5))
}

#------ tabulate_benchmarks ------------
tabulate_benchmarks <- function(in_bm, in_bmid) {

  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  print('Getting table content...')
  tbbm <- bm_msrtab(in_bm) %>%
    format_modelcompdat(typecomp=in_bmid)

  metrics_dat <- data.table(selection=character(),
                            type=character(),
                            learner_format=character(),
                            inner_folds=integer(),
                            inner_n_evals=integer(),
                            ntree = integer(),
                            numtrees = integer(),
                            alpha=numeric(),
                            mtry=integer(),
                            minnodesize=integer(),
                            fraction=numeric(),
                            samplefraction = numeric(),
                            `minor_weight|ratio`=numeric(),
                            minor_weight = numeric(),
                            ratio = numeric(),
                            npredictors=integer(),
                            outer_repeats=integer(),
                            outer_folds=integer()) %>%
    rbind(tbbm, use.names=TRUE, fill=TRUE)

  #Continue formatting table
  print('Format table...')
  metrics_dat[is.na(ntree), ntree := numtrees]
  metrics_dat[is.na(fraction), fraction := samplefraction]
  metrics_dat[, `minor_weight|ratio` := fifelse(is.na(minor_weight),
                                                as.numeric(ratio),
                                                as.numeric(minor_weight))]


  #Create model specification tables
  setup_table <- metrics_dat[, .(selection, type, learner_format,
                                 inner_folds, inner_n_evals,
                                 alpha, mtry, minnodesize, fraction,
                                 `minor_weight|ratio`, npredictors,
                                 outer_repeats, outer_folds)] %>%
    .[!grepl('(CIF)|(MAXSTAT)', learner_format), alpha := NA] %>%
    .[grepl('MAXSTAT', learner_format), minnodesize := NA] %>%
    unique(by='learner_format')

  #Create model results table
  results_table <- metrics_dat[, .(selection, learner_format,
                                   resampling_id, outer_repeats, outer_folds,
                                   time_train=round(time_train),
                                   time_predict=round(time_predict),
                                   bacc=round(bacc, 3), threshold_class,
                                   bbrier=round(bbrier, 3), auc=round(auc, 3)
  )]

  return(list(setup=setup_table, results=results_table))
}

#------ formatmisclass_bm -------------
formatmisclass_bm <- function(in_bm, in_bmid) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  print('Getting table content...')
  thresh_dt <- threshold_dat(bmres =in_bm) %>%
    format_modelcompdat(in_bmid)
  return(thresh_dt)
}


#------ ggmisclass_bm -------------
ggmisclass_bm <- function(in_threshdts) {

  #
  threshplot_datall <- rbindlist(in_threshdts) %>%
    .[, list(value_mean = mean(value)),
      by=.(threshold_class, variable, learner_format, selection)]


  ggthresh_benchmark <- ggplot(threshplot_datall[variable %in% c('sens', 'spec'),],
                               aes(x=threshold_class, y=value_mean,
                                   linetype=variable,
                                   gourp=learner_format,
                                   color=learner_format)) +
    geom_line(size=1.2) +
    #scale_color_brewer(palette='Dark2') +
    scale_linetype(labels=c('Sensitivity (true positives)',
                            'Specificity (true negatives)')) +
    scale_x_continuous(expand=c(0,0), name='Threshold') +
    scale_y_continuous(expand=c(0,0), name='Value') +
    facet_wrap(~selection) +
    theme_bw() +
    theme(panel.spacing = unit(0.75, "cm"),
          legend.title = element_blank())

  return(ggthresh_benchmark)
}


#------ tabulate_globalsummary -----
#Example
# gad_tabledis <- tabulate_globalsummary(in_filestructure = filestructure,
#                                        idvars = 'gad_id_cmj',
#                                        castvar = 'dis_m3_pyr',
#                                        castvar_num = TRUE,
#                                        weightvar = 'LENGTH_KM',
#                                        valuevar = 'predbasic800cat',
#                                        valuevar_sub = 1,
#                                        binfunc = 'manual',
#                                        binarg = c(0.1, 1, 10, 100, 1000, 10000, 100000),
#                                        na.rm=T,
#                                        tidy = FALSE)

tabulate_globalsummary <- function(outp_riveratlaspred, inp_riveratlas,
                                   inp_riveratlas_legends,
                                   idvars,
                                   castvar, castvar_num=TRUE,
                                   weightvar,
                                   valuevar, valuevarsub,
                                   binfunc=NULL, binarg=NULL, bintrans=NULL,
                                   na.rm=T, tidy=FALSE,
                                   nolake = TRUE,
                                   nozerodis = FALSE) {

  #Import global predictions
  if (is.character(outp_riveratlaspred)) {
    rivpred <- fread(outp_riveratlaspred)
  }

  #Columns to import from full network
  incols <- c('HYRIV_ID', 'INLAKEPERC', castvar, idvars, valuevar, weightvar)
  #Import global river network and join to predictions
  riveratlas <- fread_cols(file_name=inp_riveratlas,
                           cols_tokeep = incols) %>%
    .[rivpred, on='HYRIV_ID']

  #Exclude either those that intersect lakes and/or those that have zero discharge
  if (nolake) {
    riveratlas <- riveratlas[INLAKEPERC < 100,]
  }

  if (nozerodis) {
    riveratlas <- riveratlas[dis_m3_pyr > 0,]
  }


  #If global administrative boundaries were selected for idvar, get country names
  if ('gad_id_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(inp_riveratlas_legends, sheet='gad_id') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='gad_id_cmj', by.y='Country_ID') %>%
      .[, gad_id_cmj := Country_Name]
  } else if ('fmh_cl_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(inp_riveratlas_legends, sheet='fmh_cl') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='fmh_cl_cmj', by.y='MHT_ID') %>%
      .[, fmh_cl_cmj := MHT_Name]
  } else if ('tbi_cl_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(inp_riveratlas_legends, sheet='tbi_cl') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='tbi_cl_cmj', by.y='Biome_ID') %>%
      .[, tbi_cl_cmj := Biome_Name]
  } else if ('clz_cl_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(inp_riveratlas_legends, sheet='clz_cl') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='clz_cl_cmj', by.y='GEnZ_ID') %>%
      .[, clz_cl_cmj := GEnZ_Name]
  }

  #Bin castvar if needed
  if (!is.null(binfunc) & !is.null(binarg)) {
    riveratlas <- bin_dt(in_dt = riveratlas, binvar = castvar, valuevar = valuevar,
                         binfunc = binfunc, binarg = binarg, bintrans = bintrans,
                         na.rm = na.rm)
    castvar = 'bin_lformat'
    castvar_num = FALSE
  }


  #Compute overall number or weight for all combinations of castvar, valuevar and idvars
  #e.g. total number or river length for each flow state category for each country and river order
  statall <- riveratlas[if (na.rm) {!is.na(eval(valuevar))},
                        if (!is.null(weightvar)) sum(get(weightvar)) else .N,
                        by=c(eval(castvar), eval(valuevar), eval(idvars))] %>%
    rbind(
      .[, list('World', sum(V1)), by=c(eval(castvar), eval(valuevar))] %>%
        setnames(c('V1', 'V2'), c(idvars, 'V1'))
    )


  #Compute for each cast and id var, the percentage of the valuevar that is of the category of interest
  #(e.g. for each country and river order, percentage of river length that is intermittent)
  tidyperc <- statall[, 100*.SD[get(valuevar)==valuevarsub, sum(V1)]/sum(V1),
                      by=c(eval(castvar), eval(idvars))]

  #If cast variable is a numeric (e.g. Strahler Order), sort table in increasing order
  if (castvar_num) {
    tidyperc[, eval(castvar) := as.numeric(as.character(get(eval(castvar))))]
  }
  setorderv(tidyperc, cols=c(castvar, idvars))


  #Compute totals for each valuevar and idvar summed across castvar
  #e.g. (total length of rivers for each country, and of each flow regime)
  tidytotal <- statall[, list('Total intermittency (%)',
                              100*.SD[get(valuevar)==valuevarsub, sum(V1)]/sum(V1)),
                       by=eval(idvars)] %>%
    rbind(  statall[, list('Total stream length (10^3 km)', sum(V1)/1000),
                    by=c(eval(idvars))]) %>%
    setnames(c('V1', 'V2'), c(castvar, 'V1'))
  totalcols <- unique(tidytotal[, get(eval(castvar))]) #Name of total cols

  #Prepare casting formula for formatting
  castformula <- as.formula(paste0(paste(idvars, collapse='+'),
                                   '~',
                                   paste(castvar, collapse='+')))

  #Merge percentages and totals, then cast table if tidy != TRUE
  tidypercformat <- tidyperc[, eval(castvar) := as.factor(get(eval(castvar)))] %>%
    rbind(tidytotal)

  if (!tidy) {
    tidyperc_cast <- dcast(tidypercformat, castformula, value.var = 'V1') %>%
      .[, (totalcols) := lapply(.SD, function(x) fifelse(is.na(x), 0 , x)), #Replace NAs in total cols by 0
        .SDcols = totalcols]

    tidyperc_format <- rbind(tidyperc_cast[get(eval(idvars)) != 'World',],
                             tidyperc_cast[get(eval(idvars)) == 'World',])
  }

  return(tidyperc_format)
}
#------ compare_fr --------------------------------------
compare_fr <- function(inp_frdir, in_rivpred, binarg) {
  in_netpath <- file.path(inp_frdir, 'network')
  in_baspath <- file.path(inp_frdir, 'hydrobasins12')
  valuevarsub <- "1"

  net_all <- st_read(dsn = dirname(in_netpath),
                 layer = basename(in_netpath))
  net <- net_all[(net_all$rhtvs2_all_phi_qclass_SURF_BV>=10) |
        (net_all$rhtvs2_all_phi_qclass_MODULE >= 0.1),]

  bas <- st_read(dsn = dirname(in_baspath),
                 layer = basename(in_baspath)) %>%
    .[, 'HYBAS_ID', with=F]

  rivpredsub <- merge(in_rivpred, bas, by.x="HYBAS_L12", by.y="HYBAS_ID", all.x=F) %>%
    .[, UPLAND_SKM := round(UPLAND_SKM)]

  binlabels <- label_manualbins(binarg=binarg,
                                minval=min(net$rhtvs2_all_phi_qclass_SURF_BV))

  tidyperc_fr <- formathistab(in_dt = net,
                              castvar = "rhtvs2_all_phi_qclass_MODULE",
                              valuevar = "INT_RF_txt_V1",
                              valuevarsub = valuevarsub,
                              weightvar = "rhtvs2_all_phi_qclass_LONG_",
                              binfunc = 'manual',
                              binarg =  binarg,
                              binlabels = binlabels,
                              datname = 'Snelder') %>%
    .[, binsumlength := binsumlength/1000]

  tidyperc_riv  <- formathistab(in_dt = rivpredsub,
                                castvar = 'dis_m3_pyr',
                                valuevar = 'predbasic800cat',
                                valuevarsub = valuevarsub,
                                weightvar = 'LENGTH_KM',
                                binfunc = 'manual',
                                binarg =  binarg,
                                binlabels = binlabels,
                                datname = 'Global')

  datmerge <- rbind(tidyperc_fr, tidyperc_riv) %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]
  print('Percentage intermittence for France')
  print(datmerge[, weighted.mean(perc, binsumlength), by=dat])

  return(
    ggcompare(datmerge, binarg) + theme(axis.title.y = element_blank(),
                                        axis.title.x = element_blank(),
                                        legend.title = element_blank())
  )
}

#------ compare_us ----------------
compare_us <- function(inp_usresdir, inp_usdatdir, in_rivpred, binarg) {

  in_netpath_hr <- file.path(inp_usdatdir,  'NHDhr_attris.csv')
  in_netpath_mr <- file.path(inp_usdatdir, 'NHDmr_attris.csv')
  in_baspath <- file.path(inp_usresdir, 'hydrobasins12')
  valuevarsub <- "1"

  #GEt NHD high and medium resolution
  #55800: Artificial path
  #46000: Stream/River
  #46003: Stream/River: Hydrographic Category = Intermittent
  #46006: Stream/River: Hydrographic Category = Perennial
  #46007: Stream/River: Hydrographic Category = Ephemeral
  keepfcodes <- c(55800, 46000, 46003, 46007, 46006) #Flow line types to keep

  #Read and format NHD high res
  nethr <- fread(in_netpath_hr,
                 colClasses = c('character', 'character', 'numeric', 'character',
                                'integer', 'numeric', 'numeric'))
  nethr[, HUC8 := substr(ReachCode, 1, 8)]

  #Read NHD medium res
  mrcolclasses <- list(character=c('COMID', 'REACHCODE'),
                       numeric=c('LENGTHKM', 'TotDASqKM', 'QE_MA'),
                       integer = c('FTYPE', 'FCODE', 'StreamOrde'))
  netmr <- fread(in_netpath_mr,
                 colClasses = list(character=c('COMID', 'REACHCODE'),
                                   numeric=c('LENGTHKM', 'TotDASqKM', 'QE_MA'),
                                   integer = c('FTYPE', 'FCODE', 'StreamOrde'))) %>%
    .[,unlist(mrcolclasses), with=F]
  netmr[, `:=`(HUC8 = substr(REACHCODE, 1, 8),
               QE_MAm3 = QE_MA*0.028316847)]

  netmrsubhr <- netmr[HUC8 %in% unique(nethr$HUC8),]

  #Get HydroSHEDS basins that overlap with selected NHD HUC8s
  bas <- st_read(dsn = dirname(in_baspath),
                 layer = basename(in_baspath)) %>%
    .[, c('HYBAS_ID', 'HUC8'), with=F]

  #Join HydroSHEDS basins with RiverATLAS network and subselect network to match NHD selection
  rivpredbas <- merge(in_rivpred, bas,
                      by.x="HYBAS_L12", by.y="HYBAS_ID", all.x=F) #To match full US
  rivpredsubhr <- rivpredbas[, UPLAND_SKM := round(UPLAND_SKM)] %>%
    .[HUC8 %in% unique(nethr$HUC8),]
  rivpredsubmr <- rivpredbas[HUC8 %in% unique(netmr$HUC8),] %>%
    .[dis_m3_pyr > 0,] #Remove segments with 0 discharge

  #Compare some statistics
  nethrlen <- nethr[FCode %in% keepfcodes,
                    round(sum(LengthKM, na.rm=T)/(10^6), 1)]

  print(paste0('Total length of streams in NHDplus highres:' ,
               nethrlen,
               ' million km'))
  print(paste0(round(nethrlen/netmrsubhr[FCODE %in% keepfcodes,
                                         sum(LENGTHKM, na.rm=T)/(10^6)],
                     1),
               ' times the length in HDplus res'))
  print(paste0(round(nethrlen/rivpredsubhr[, sum(LENGTH_KM)/(10^6)],
                     1),
               ' times the length in HydroSHEDS'))


  print(paste0('% of stream length in NHDplus highres with DA < 10 km2:' ,
               round(
                 100 * nethr[FCode %in% keepfcodes &
                               !is.na(TotDASqKm) & TotDASqKm < 10,
                             sum(LengthKM, na.rm=T)]/
                   nethr[FCode %in% keepfcodes & !is.na(TotDASqKm),
                         sum(LengthKM, na.rm=T)]
               ), ' %'))

  #Total range in intermittency
  print(paste0(
    'Range of prvalence of intermittency in NHDplus high resolution,',
    'depending on whether unclassified reaches are counted as fully intermittent or perennial: ',
    round(100*nethr[FCode %in% c(46003, 46007), sum(LengthKM)]/
            nethr[FCode %in% keepfcodes, sum(LengthKM)]),
    '-',
    round(100*nethr[FCode %in% c(46003, 46007, 46000, 55800), sum(LengthKM)]/
            nethr[FCode %in% keepfcodes, sum(LengthKM)]),
    '%'))

  print(paste0(
    'Range of prevalence of intermittency in NHDplus medium resolution,',
    'depending on whether unclassified reaches are counted as fully intermittent or perennial: ',
    round(100*netmrsubhr[FCODE %in% c(46003, 46007), sum(LENGTHKM)]/
            netmrsubhr[FCODE %in% keepfcodes, sum(LENGTHKM)]),
    '-',
    round(100*netmrsubhr[FCODE %in% c(46003, 46007, 46000, 55800), sum(LENGTHKM)]/
            netmrsubhr[FCODE %in% keepfcodes, sum(LENGTHKM)]),
    '%'))


  netmr_o10 <- netmr[((TotDASqKM >= 10)| (QE_MA*0.028316847) >= 0.1),]

  print(paste0(
    'Range of prevalence of intermittency in NHDplus medium resolution >= 10 km2,',
    'depending on whether unclassified reaches are counted as fully intermittent or perennial: ',
    round(100*netmr_o10[FCODE %in% c(46003, 46007), sum(LENGTHKM)]/
            netmr_o10[FCODE %in% keepfcodes, sum(LENGTHKM)]),
    '-',
    round(100*netmr_o10[FCODE %in% c(46003, 46007, 46000, 55800), sum(LENGTHKM)]/
            netmr_o10[FCODE %in% keepfcodes, sum(LENGTHKM)]),
    '%'))

  print(paste0(
    'Total estimated prevalence in HydroSHEDS: ',
    rivpredsubmr[predbasic800cat==1, sum(LENGTH_KM)]/
      rivpredsubmr[, sum(LENGTH_KM)]
  ))


  #Check the percentage of each type of line in the NHD by HUC8 and Drainage area
  netmr_fsub <- copy(netmr[FCODE %in% keepfcodes,])
  netmr_fsub[, HUC8sum := sum(LENGTHKM), by=HUC8]
  FCode_HUC8 <- netmr_fsub[FCODE %in% keepfcodes,
                           as.integer(100*.SD[,sum(LENGTHKM)]/max(HUC8sum)),
                           by=.(HUC8, FCODE)]

  netmrbinned <- bin_dt(in_dt=netmr[TotDASqKM > 0,],
                        binvar='QE_MAm3',
                        valuevar='intermittent',
                        binfunc='manual',
                        binarg=binarg)
  # netmrbinned[FCode %in% c(55800, 46000, 46003, 46006, 46007),
  #           DAsum := .N, by=bin_lformat]
  # FCode_DA <- netmrbinned[FCode %in% c(55800, 46000, 46003, 46006, 46007),
  #                       as.integer(100*max(.SD[, .N]/DAsum)), by=.(bin_lformat, FCode)]

  #Create intermittency variable
  netmr_fsub[, intermittent := fifelse(FCODE %in% c(46003, 46007), '1', '0')]
  netmr_fsub[FCODE %in% c(46000, 55800), intermittent := -1]
  netmr_fsub[!(FCODE %in% keepfcodes), intermittent := NA]

  #Get bin labels
  binlabels <- label_manualbins(binarg=binarg,
                                minval=min(netmr_fsub$QE_MAm3))

  # --------------------- Compute HUC8 comparison ----------------------------------
  tidyperc_nhdHUC <- netmr_fsub[!is.na(intermittent),][
    ,
    formathistab(in_dt = .SD,
                 castvar = "QE_MAm3",
                 valuevar = "intermittent",
                 valuevarsub = as.numeric(as.character(valuevarsub)),
                 weightvar = "LENGTHKM",
                 binfunc = 'manual',
                 binarg =  binarg,
                 binlabels = binlabels,
                 datname = 'U.S. NHD'),
    by=HUC8]

  tidyperc_rivHUC  <- rivpredsubmr[
    , formathistab(in_dt = .SD,
                   castvar = 'dis_m3_pyr',
                   valuevar = 'predbasic800cat',
                   valuevarsub = valuevarsub,
                   weightvar = 'LENGTH_KM',
                   binfunc = 'manual',
                   binarg =  binarg,
                   binlabels = binlabels,
                   datname = 'Global predictions'),
    by=HUC8]

  hucmajcl <- rivpredsubmr[, .N,
                           by=.(clz_cl_cmj, HUC8)][
                             , .SD[which.max(N), clz_cl_cmj], by=HUC8] %>%
    setnames('V1', 'majcl')

  # setnames(tidyper_nhdHUC, paste0(names(tidyper_nhdHUC), 'NHD'))
  # tidyHUCmerge <- merge(tidyperc_rivHUC, tidyperc_nhdHUC,
  #                       by.x=c('HUC8','bin'), by.y=c('HUC8NHD','binNHD'),
  #                       all.x = T, all.y=T) %>%
  # ggplot(tidyHUCmerge, aes(x=perc, y=percNHD, color=factor(majcl))) +
  #   geom_point() +
  #   facet_wrap(~binformat)

  tidyperc_hucbind <- rbind(tidyperc_nhdHUC, tidyperc_rivHUC) %>%
    merge(hucmajcl, by='HUC8', all.x=T)

  clhucbox <- ggplot(tidyperc_hucbind, aes(x=factor(majcl), y=perc, fill=dat)) +
    geom_boxplot() +
    facet_wrap(~binformat)


  # --------------------- Compute national hist. comparison ----------------------------------
  #Compute bin statistics ncluding "artificial flow path"
  tidyperc_usall <- formathistab(in_dt = netmr_fsub[!is.na(intermittent),],
                                 castvar = "QE_MAm3",
                                 valuevar = "intermittent",
                                 valuevarsub = as.numeric(as.character(valuevarsub)),
                                 weightvar = "LENGTHKM",
                                 binfunc = 'manual',
                                 binarg =  binarg,
                                 binlabels = binlabels,
                                 datname = 'U.S. NHD') %>%
    .[, binsumlength := binsumlength]

  #Compute bin statistics excluding "artificial flow path"
  tidyperc_usnoartificial <- formathistab(in_dt = netmr_fsub[!is.na(intermittent) &
                                                               intermittent != -1,],
                                          castvar = "QE_MAm3",
                                          valuevar = "intermittent",
                                          valuevarsub = valuevarsub,
                                          weightvar = "LENGTHKM",
                                          binfunc = 'manual',
                                          binarg =  binarg,
                                          binlabels = binlabels,
                                          datname = 'U.S. NHD') %>%
    .[, binsumlength := binsumlength]

  #Compute bin statistics for river networks
  tidyperc_riv  <- formathistab(in_dt = rivpredsubmr,
                                castvar = 'dis_m3_pyr',
                                valuevar = 'predbasic800cat',
                                valuevarsub = valuevarsub,
                                weightvar = 'LENGTH_KM',
                                binfunc = 'manual',
                                binarg =  binarg,
                                binlabels = binlabels,
                                datname = 'Global predictions')


  #Merge statistics
  datmerge_all <- rbind(tidyperc_usall, tidyperc_riv) %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]

  datmerge_noartificial <- rbind(tidyperc_usnoartificial, tidyperc_riv) %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]

  #Create plot
  return(
    ggcompare(datmerge_all, binarg=binarg) +
      geom_bar(data=datmerge_noartificial, alpha=0.5,
               stat='identity', pattern = 'stripe', position='dodge') +
      theme(axis.title.y = element_blank(),
            legend.title = element_blank())
  )
}


#------ qc_pnw ------------
qc_pnw <- function(inp_pnwresdir, in_rivpred, interthresh=0.5) {
  in_refpts <- file.path(inp_pnwresdir, 'StreamflowPermObs_final')
  in_fulldat <- file.path(inp_pnwresdir, 'StreamflowPermObs_sub')

  valuevarsub <- "1"

  #Georeferenced/Snapped points to RiverATLAS network after removing duplicate observations at single sites
  refpts <- st_read(dsn=dirname(in_refpts),
                    layer=basename(in_refpts)
                    ) %>%
    setDT %>%
    setkey('OBJECTID')

  #All observations (excluding those that were not kept by PROSPER (aside from more recent osb, see python code for details)
  #with duplicate records for single locations (multiple observations for different dates)
  fulldat <-  st_read(dsn=dirname(in_fulldat),
                      layer=basename(in_fulldat)
  ) %>%
    setDT %>%
    setkey('OBJECTID')

  #fulldat[OBJECTID %in% refpts$OBJECTID, .N]

  #Join points and compute # of observations, min and max months of observations
  refpts_full <- merge(refpts, fulldat[, c('OBJECTID', 'dupligroup'), with=F],
                       by='OBJECTID', all.y=F) %>%
    .[fulldat, on='dupligroup', allow.cartesian=T] %>%
    .[!is.na(distatlas),]

  refpts_stats <- refpts_full[, `:=`(nobs=.N,
                                     minmonth = min(Month),
                                     maxmonth = max(Month)),
                              by = dupligroup] %>%
    .[!duplicated(dupligroup),]

  #Merge points with rivernetwork by HYRIV_ID
  refpts_join <- merge(refpts_stats,
                       in_rivpred[, .(HYRIV_ID, HYBAS_L12, predbasic800,
                                      predbasic800cat)],
                       by='HYRIV_ID', all.y=F) %>%
    .[, refinter := fifelse(Category == 'Non-perennial', 1, 0)] #Assign ephemeral and intermittent categories to 1, perennial to 0 for PNW obs


  #Consider that if one location is non-perennial, that reach is not perennial
  refpts_reach <- refpts_join[
    , refinter_reach := as.numeric(any(refinter == 1)), by=HYRIV_ID] %>%
    setorder(HYRIV_ID, -refinter) %>%
    .[!duplicated(HYRIV_ID),] %>%
    .[, IPR := predbasic800 - refinter] #Compute prediction error

  #Compute BACC
  refpts_reach[, `:=`(truth = factor(as.character(refinter_reach),levels=c('0', '1')),
                     response = factor(predbasic800cat,levels=c('0','1'))
  )]
  pnw_measures <- refpts_reach[, list(
    bacc = mlr3measures::bacc(truth, response),
    auc = mlr3measures::auc(truth, predbasic800, positive='1'),
    ce = mlr3measures::ce(truth, response),
    sen = mlr3measures::sensitivity(truth, response, positive='1'),
    spe = mlr3measures::specificity(truth, response, positive='1')
  )]

  print(paste0('Balanced accuracy of predictions based on PROSPER: ',
               round(pnw_measures$bacc, 3)))
  print(paste0('AUC based on PROSPER: ',
               round(pnw_measures$auc, 3)))
  print(paste0('Misclassification rate based on PROSPER: ',
               round(pnw_measures$ce, 3)))
  print(paste0('Sensitivity based on PROSPER: ',
               round(pnw_measures$sen, 3)))
  print(paste0('Specificity rate based on PROSPER: ',
               round(pnw_measures$spe, 3)))
  print(paste0('Median discharge of observations for PROSPER',
               refpts_reach[, round(median(HYDROSHEDSdis), 2)]))

  #Compute relative position of observation point on HydroSHEDS reach line
  refpts_reach[, RATIOLENGTH := min(c(1, (fromM*95)/LENGTH_KM)), by=OBJECTID]

  #Write data out to points
  predtowrite <- merge(
    st_read(dsn = inp_pnwresdir,
            layer = basename(in_refpts))[,c('OBJECTID', 'Shape')],
    refpts_reach[,-c('Shape','i.Shape'), with=F],
    by='OBJECTID', all.x=F, all.y=T) %>%
    setnames(old=c("HYDROSHEDSdis", "HYDROSHEDSDA",
                   "predbasic800cat", "predbasic800"),
             new=c("LINEdis", "LINEDA", "predcat", "predp")) %>%
    setnames(new=unlist(lapply(names(.), function(x) gsub('[.]','_', x))))

  st_write(obj=predtowrite[,c('OBJECTID', 'IPR',
                              "LINEdis", "LINEDA",
                              "predcat", "predp",
                               'nobs', 'refinter_reach',
                               'pop_ct_usu')],
           dsn = dirname(inp_pnwresdir),
           layer=paste0('pnwobs_IPR',
                        format(Sys.time(), '%Y%m%d%H'), '.shp'),
           driver='ESRI Shapefile')

  #Scatterpoint of x = drainage area, y = predprobability - refinter_reach
  refpts_reachmelt <- melt(
    refpts_reach[, .(OBJECTID, IPR,
                    HYDROSHEDSdis, RATIOLENGTH,
                    nobs, refinter_reach,
                    pop_ct_usu)],
    id.vars = c('OBJECTID', 'IPR', 'refinter_reach'))


  levels(refpts_reachmelt$variable) <- c(
    'HydroSHEDS Discharge (m3)',
    'Length ratio',
    '# of field obs.',
    'Watershed population density'
  )

  colorpal <- c('#1f78b4', '#ff7f00')
  rectdf <- data.table(
    xmin=rep(0, 4),
    xmax=rep(Inf, 4),
    ymin=c(-1, interthresh-1, 0, interthresh),
    ymax=c(interthresh-1, 0, interthresh, 1),
    fillpal = rep(colorpal, 2)
  )

  pnw_qcplot <- ggplot(refpts_reachmelt) +
    geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
                               ymin=ymin, ymax=ymax, fill=fillpal),
              alpha=1/4) +
    scale_fill_manual(values=colorpal,
                      name='Predicted regime',
                      labels = c('Perennial', 'Intermittent')) +
    geom_point(aes(x=value, y=IPR, color=factor(refinter_reach)),
               alpha=1/5) +
    geom_hline(yintercept=0, alpha=1/2) +
    new_scale_fill() +
    # geom_smooth(aes(x=value, y=IPR, color=factor(refinter_reach)),
    #             method='loess', span=1) + #, formula = "y ~ s(x, k=2)") +
    scale_x_sqrt() +
    # scale_x_continuous(trans='log1p',
    #                    breaks = c(1, 5, 10, 100, 1000, 10000, 100000,1000000),
    #                    labels = scientific_10) +
    geom_hline(yintercept=0) +
    scale_color_manual(values=c('#1f78b4', '#ff7f00'),
                       name='Observed regime',
                       labels = c('Perennial', 'Intermittent')) +
    coord_cartesian(expand=FALSE, clip='off') +
    labs(x='Value', y='Intermittency Prediction Residuals (IPR)') +
    theme_classic() +
    theme(legend.position = c(0.92, 0.1),
          legend.background = element_blank()) +
    facet_wrap(~variable, scales ='free_x', ncol=2)

  return(pnw_qcplot)
}


#------ qc_onde ------
qc_onde <- function(inp_ondedatdir, inp_onderesdir, inp_riveratlas,
                    in_rivpred, interthresh=0.5) {
  in_refpts <- file.path(inp_onderesdir, 'obs_finalwgs')
  in_fulldat <- file.path(inp_ondedatdir, 'onde_france_merge.csv')

  valuevarsub <- "1"

  #Georeferenced/Snapped points to RiverATLAS network after removing duplicate observations at single sites
  refpts <- st_read(dsn = inp_onderesdir,
                    layer = basename(in_refpts)) %>%
    setDT %>%
    setnames(gsub('(F_)|_', '', names(.))) %>%
    setkey('CdSiteHydro')

  #All observations with duplicate records for single locations (multiple observations for different dates)
  fulldat <- fread(in_fulldat) %>%
    setnames(gsub('<|>', '', names(.))) %>%
    setkey('CdSiteHydro')

  #Join formmatted and snapped points to full dataset of observations (keeping all
  #observation instances)
  refpts_colstokeep <- c("CdSiteHydro",
                         names(refpts)[!(names(refpts) %in% names(fulldat))])

  refpts_full <- merge(fulldat, refpts[, refpts_colstokeep, with=F], all.x=F) %>%
    .[!is.na(fromM),]

  refpts_full[, Month := format.Date(DtRealObservation, '%m')]

  #Get statistics for each site: total number of observations + N of each type of obs
  #min and max year and months of observation
  obstats <- refpts_full[, list(nobs=.N,
                                minmonth = min(as.numeric(as.character(Month))),
                                maxmonth = max(as.numeric(as.character(Month))),
                                minyear = min(Annee),
                                maxyear = max(Annee)),
                         by = .(CdSiteHydro)] %>%
    merge(dcast(refpts_full[, .N, by=.(CdSiteHydro, LbRsObservationNat)],
                CdSiteHydro~LbRsObservationNat, value.var='N'),
          by='CdSiteHydro')

  setnames(obstats,
           old=c('Assec',
                 'Ecoulement visible',
                 'Ecoulement non visible',
                 'Observation impossible'),
           new=c('Dry',
                 'Flow',
                 'NoFlow',
                 'ObsImp'))

  obstats[is.na(Dry), Dry := 0]
  obstats[is.na(Flow), Flow := 0]
  obstats[is.na(NoFlow), NoFlow := 0]
  obstats[is.na(ObsImp), ObsImp := 0]

  #Re-merge with site characteristics attributes
  uniquesites <- refpts_full[, .(CdSiteHydro, LbSiteHydro,
                                 NomEntiteHydrographique, Etat,
                                 HYRIVIDjoinedit, fromM,
                                 HYDROSHEDSDA, HYDROSHEDSdis,
                                 POINTdis, POINTDA)] %>%
    unique

  obsattri <- merge(obstats, uniquesites, by='CdSiteHydro') %>%
    .[nobs > 20,]

  obsattri[, sum(nobs)]


  #Compute (observation point/hydrosheds reach pour point) ratio of discharge and drainage area
  obsattri[, `:=`(RATIODA = POINTDA/(HYDROSHEDSDA*100),
                  RATIOdis = POINTdis/(HYDROSHEDSdis*10^5),
                  refinter_perc = (Dry + NoFlow)/nobs)]
  obsattri[, refinter := as.numeric(refinter_perc > 0)]

  #Get population density and agricultural cover
  riveratlas <- fread_cols(inp_riveratlas,
                           cols_tokeep = c('HYRIV_ID',
                                           'pop_ct_csu', 'pop_ct_usu',
                                           'ire_pc_cse', 'ire_pc_use',
                                           'glc_pc_c16', 'glc_pc_u16')
  )

  #Merge points with rivernetwork by HYRIV_ID
  refpts_join <- merge(obsattri,
                       in_rivpred[, .(HYRIV_ID, HYBAS_L12, LENGTH_KM,
                                      predbasic800, predbasic800cat, dis_m3_pyr)],
                       by.x='HYRIVIDjoinedit', by.y='HYRIV_ID', all.y=F) %>%
    merge(riveratlas, by.x='HYRIVIDjoinedit', by.y='HYRIV_ID', all.y=F) %>%
    .[, IPR := predbasic800 - refinter] #Compute prediction error

  #Consider that if one location is non-perennial, that reach is not perennial
  refpts_reach <- refpts_join[
    , refinter_reach := as.numeric(any(refinter == 1)), by=HYRIVIDjoinedit] %>%
    setorder(HYRIVIDjoinedit, -refinter) %>%
    .[!duplicated(HYRIVIDjoinedit),] %>%
    .[, IPR := predbasic800 - refinter] #Compute prediction error

  refpts_reach[, .N]
  refpts_reach[, .N, by=refinter_reach]

  #Compute BACC
  refpts_reach[, `:=`(truth = factor(as.character(refinter_reach),levels=c('0', '1')),
                     response = factor(predbasic800cat,levels=c('0','1'))
                     )]
  onde_measures <- refpts_reach[, list(
    bacc = mlr3measures::bacc(truth, response),
    auc = mlr3measures::auc(truth, predbasic800, positive='1'),
    ce = mlr3measures::ce(truth, response),
    sen = mlr3measures::sensitivity(truth, response, positive='1'),
    spe = mlr3measures::specificity(truth, response, positive='1')
  )]

  print(paste0('Balanced accuracy of predictions based on ONDE: ',
               round(onde_measures$bacc, 3)))
  print(paste0('AUC based on ONDE: ',
               round(onde_measures$auc, 3)))
  print(paste0('Misclassification rate based on ONDE: ',
               round(onde_measures$ce, 3)))
  print(paste0('Sensitivity based on ONDE: ',
               round(onde_measures$sen, 3)))
  print(paste0('Specificity rate based on ONDE: ',
               round(onde_measures$spe, 3)))
  print(paste0('Median discharge of observations',
               refpts_reach[, round(median(HYDROSHEDSdis), 2)]))

  #Compute relative position of observation point on HydroSHEDS reach line
  refpts_join[, RATIOLENGTH := (fromM/1000)/LENGTH_KM]

  #Write data out to points
  predtowrite <- merge(st_read(dsn = inp_onderesdir,
                               layer = basename(in_refpts))[, "F_CdSiteHydro_"],
                       refpts_join,
                       by.x = "F_CdSiteHydro_", by.y = 'CdSiteHydro') %>%
    setnames(old=c("HYDROSHEDSdis", "HYDROSHEDSDA",
                   "predbasic800cat", "predbasic800"),
             new=c("LINEdis", "LINEDA", "predcat", "predp"))

  write_sf(predtowrite, file.path(dirname(inp_onderesdir), 'ondeobs_IPR4.shp'))

  #Remove those whose ID doesn't fit what they were snapped to? (to inspect another time?)
  refpts_joinsub <- refpts_join[RATIOLENGTH < 1.01 & RATIODA< 1.01,]

  #Scatterpoint of x = drainage area, y = predprobability - refinter
  refpts_joinmelt <- melt(
    refpts_joinsub[, .(CdSiteHydro, IPR,
                       HYDROSHEDSdis, RATIOdis,
                       RATIODA, RATIOLENGTH,
                       nobs, refinter, refinter_perc,
                       pop_ct_usu)],
    id.vars = c('CdSiteHydro', 'IPR', 'refinter'))


  levels(refpts_joinmelt$variable) <- c(
    'HydroSHEDS Discharge (m3)',
    #'HydroSHEDS Drainage area (DA, km2)',
    'Discharge ratio',
    'Drainage area ratio',
    'Length ratio',
    '# of field obs.',
    'Percentage intermittent observations',
    'Watershed population density'
    )

  colorpal <- c('#1f78b4', '#ff7f00')
  rectdf <- data.table(
    xmin=rep(0, 4),
    xmax=rep(Inf, 4),
    ymin=c(-1, interthresh-1, 0, interthresh),
    ymax=c(interthresh-1, 0, interthresh, 1),
    fillpal = rep(colorpal, 2)
  )

  onde_qcplot <- ggplot(refpts_joinmelt) +
    geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
                               ymin=ymin, ymax=ymax, fill=fillpal),
              alpha=1/4) +
    scale_fill_manual(values=colorpal,
                      name='Predicted regime',
                      labels = c('Perennial', 'Intermittent')) +
    geom_point(aes(x=value, y=IPR, color=factor(refinter)),
               alpha=1/5) +
    geom_hline(yintercept=0, alpha=1/2) +
    new_scale_fill() +
    geom_smooth(aes(x=value, y=IPR, color=factor(refinter)),
                method='gam', formula = y ~ s(x, k=4)) +
    scale_x_sqrt() +
    geom_hline(yintercept=0) +
    scale_color_manual(values=c('#1f78b4', '#ff7f00'),
                       name='Observed regime',
                       labels = c('Perennial', 'Intermittent')) +
    coord_cartesian(expand=FALSE, clip='off') +
    labs(x='Value', y='Intermittency Prediction Residuals (IPR)') +
    theme_classic() +
    theme(legend.position = c(0.8, 0.1),
          legend.background = element_blank(),
          panel.spacing = unit(1.5, "lines")) +
    facet_wrap(~variable, scales ='free_x', ncol=2)

  return(onde_qcplot)
}



########## END #############
