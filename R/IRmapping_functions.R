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
#' zero_lomf(test1)
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
#' Read and pre-format GRDC data
#'
#' Reads text file of daily discharge data for a single GRDC station.
#' Creates columns for year, month, and date of last non-zero flow day +
#' computes yearly number of days of missing data
#'
#' @param path (character) path to the text file of daily discharge data in
#'   standard GRDC format.
#'
#' @return \link[data.table]{data.table} of daily discharge data with additional columns
#'
#' @export
readformatGRDC<- function(path) {
  #extract GRDC unique ID by formatting path
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  
  #Read GRDC text data
  gaugetab <- cbind(fread(path, header = T, skip = 40, sep=";",
                          colClasses = c('character', 'character', 'numeric',
                                         'numeric', 'integer')),
                    GRDC_NO = gaugeno)%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)
  
  #Format data
  gaugetab[, `:=`(year = as.numeric(substr(dates, 1, 4)), #Create year column
                  month = as.numeric(substr(dates, 6, 7)), #create month column
                  integervalue = fifelse(Original == round(Original), 1, 0) #Flag integer discharge values]
  )]
  
  #For each record, compute date of last non-zero flow day
  gaugetab[, prevflowdate := gaugetab[zero_lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA] #If non-zero flow, set prevflowdate to NA
  
  #Compute number of missing days per year, excluding NoData values
  gaugetab[!(Original %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Original)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']
  
  return(gaugetab)
}

#------ readformatGSIMmon -----------------
#' Read and pre-format GSIM monthly data
#'
#' Reads text file of monthly discharge indices for a single GSIM station.
#' Creates columns for year, month, and season (to join with seasonal indices)
#' Estimate the minimum possible number of monthly zero-flow days based on indices
#'
#' compute yearly number of days of missing data
#'
#' @param path (character) path to the text file of monthly discharge data in
#'   standard GSIM format (.mon file).
#'
#' @details the monthly minimum number of zero flow days n0 is estimated as follows: \cr
#'     1. If MAX == 0, then n0 = n.available
#'     2. else if MIN7 == 0, then n0 = 7
#'     3. else if MIN == 0, then n0 = 1,
#'     4. else n0 = 0
#'
#'     For more information on GSIM, see [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/).
#'
#' @return \link[data.table]{data.table} of monthly indices, with additional
#' attributes including the estimated minimum number of zero-flow days "mDur_minmo"
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
readformatGSIMmon <- function(path) {
  #extract GSIM unique ID by formatting path
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  
  #Read GSIM text data
  gaugetab <- cbind(fread(path, header=T, skip = 21, sep=",",
                          colClasses=c('Date', rep('numeric', 8),
                                       'integer', 'integer')),
                    gsim_no = gaugeno) %>%
    setnames(new=gsub('[\\\t]|["]', '', names(.))) %>% #Remove tab in field names
    setorder(date) %>%
    .[, ':='(year = as.numeric(substr(date, 1, 4)), #Create year column
             month = as.numeric(substr(date, 6, 7)))] %>% #Create month column
    .[data.table(month=c(12, 1:11),
                 season=rep(c('DJF', 'MAM', 'JJA', 'SON'), each=3)),
      on='month']  #Create season column (first create a data.table with two columns (month, season) then join it)
  
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
#' Read and pre-format GSIM seasonal data
#'
#' Reads text file of seasonal discharge indices for a single GSIM station.
#' Creates columns for year, month, and season
#' Estimate the minimum possible number of monthly zero-flow days based on indices
#'
#' @param path (character) path to the text file of seasonal discharge data in
#'   standard GSIM format (.mon file).
#'
#' @details the monthly minimum number of zero flow days n0 is estimated with
#' percentile statistics like:
#'     if P90==0, then n0 = as.integer(floor(0.9*n.available)),
#'     if else P80==0, then n0 = as.integer(floor(0.8*n.available)),
#'     if else P70 ...
#'     ...
#'     else if MIN7 == 0, then n0 = 7
#'     else if MIN == 0, then n0 = 1,
#'     else n0 = 0
#'
#' For more information on GSIM, see [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/).
#'
#' @return \link[data.table]{data.table} of monthly indices, with additional attributes including
#' the estimated minimum number of zero-flow days "mDur_minmo"
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
readformatGSIMsea <- function(path) {
  #extract GSIM unique ID by formatting path
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  
  #Read GSIM text data
  gaugetab <- cbind(fread(path, header=T, skip = 21, sep=",",
                          colClasses=c('Date', rep('numeric', 17),
                                       'integer', 'integer')),
                    gsim_no = gaugeno) %>%
    setnames(new=gsub('[\\\t]|["]', '', names(.))) %>% #Remove tab in field names
    setorder(date) %>%
    .[, ':='(year = as.numeric(substr(date, 1, 4)), #Create year column
             month = as.numeric(substr(date, 6, 7)))] %>% #Create month column
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
#' Flag GRDC outliers
#'
#' Flag potential outliers in daily discharge records for a given GRDC gauging
#' station following the criteria developed for GSIM by
#' [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/).
#'
#' @param in_gaugetab \link[data.table]{data.table} containing formatted daily
#' discharge record from GRDC gauging station (as formatted by \code{\link{readformatGRDC}}.
#'
#' @details Criteria to flag a daily discharge value (Qt) as a potential outlier include:
#' \itemize{
#'   \item Negative values (Qt < 0)
#'   \item At least ten identical consecutive discharge values (for Qt > 0)
#'   \item |log(Qt + 0.01) - mean| are larger than the mean values of log(Q + 0.01)
#'   plus or minus 6 times the standard deviation of log(Q + 0.01) computed for
#'   that calendar day for the entire length of the series. The mean and SD are
#'   computed for a 5-day window centred on the calendar day to ensure that a
#'   sufficient amount of data is considered. The log-transformation is used to
#'   account for the skewness of the distribution of daily streamflow values.
#'   \item Qt for which Original != Calculated discharge in GRDC record
#' }
#'
#' @return \link[data.table]{data.table} of daily discharge records with additional
#' columns for outlier flags
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
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
#' Plot a GRDC time series
#'
#' Creates a plot of daily discharge ({m^3}/s) for a GRDC gauging station,
#' with flags for 0-flow values and potential outliers.
#' Save plot to png if path is provided.
#'
#' @param GRDCgaugestats_record \link[data.table]{data.table} of formatted daily
#' discharge records for a single GRDC station. In this project, e.g. the output
#' from \code{\link{comp_GRDCdurfreq}}. Must contain at least five columns:
#' \code{GRDC_NO, dates, Original, flag_mathis, missingdays} \cr
#' Alternatively, a named list or vector with
#' a column called "path" towards a standard GRDC text file containing daily discharge
#' records for a single gauging station.
#' @param outpath (character) path for writing output png plot (default is no plotting).
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details the output graphs show the time series of daily streamflow values
#' for the station. For the flagging criteria, see documentation for \code{\link{flagGRDCoutliers}}.
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Individual points show daily discharge values (in {m^3}/s).
#'   \item blue lines link daily values (which may result in unusual patterns due to missing years).
#'   \item red points are zero-flow flow values.
#'   \item green points are non-zero flow daily values statistically flagged as potential outliers .
#'   \item black points are zero-flow values flagged as potential outliers.
#' }
#'
#' @return plot
#'
#' @export
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
#' Plot a GSIM time series
#'
#' Creates a plot of monthly discharge ({m^3}/s) for a GSIM gauging stations
#' Save plot to png if path is provided.
#'
#' @param GSIMgaugestats_record a data containing structure  with
#' a column called "path" towards a standard montly GSIM text file containing.
#' In this project, e.g. the output from \code{\link{comp_GSIMdurfreq}}.
#' @param outpath (character) path for writing output png plot (default is no plotting).
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details Daily streamflow records from GSIM stations are unavailable. Therefore,
#' the graph shows the following:
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Blue points: mean monthly discharge
#'   \item Light blue background shading: mean ± 2SD monthly discharge
#'   \item Black points: minimum and maximum monthly discharge
#'   \item Red points show minimum monthly discharge values equal to 0
#'   \item Purple points show months for which all daily discharge values are equal to 0.
#' }
#'
#' @return plot
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
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
#' Check zero-flow values in GRDC discharge record
#'
#' Plot and subset a given period (in days) on each side of every zero-flow
#' period in record for a given GRDC station. The output is a plot with a panel
#' showing discharge for each zero-flow period.
#'
#' @param GRDCstatsdt formatted data.table including a "path" column to access
#' GRDC standard daily discharge record text file for the gauging station of interest.
#' @param in_GRDC_NO GRDC unique identifier of the station to be investigated.
#' @param period (integer) number of days on each side of zero-flow values to subset.
#' @param yearthresh (integer) minimum year from which to analyze discharge record
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param in_scales (character) should panel scales be fixed ("fixed", the default),
#' free ("free"), or free in one dimension ("free_x", "free_y")?
#' @param labelvals (logical) whether to label individual discharge values
#'
#' @return \link[data.table]{data.table} subset of GRDCstatsdt, including only
#' records for zero-flow days +- period
#'
#' @export
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
#' Alluvium plot of gauge count over time
#'
#' Utility function: create an [alluvial plot](https://corybrunson.github.io/ggalluvial/)
#' of the number of streamflow gauging stations over time.
#'
#'
#' @param dtformat formatted data.table (generated in \code{\link{analyzemerge_gaugeir}}).
#' @param alluvvar grouping ID to identify unique records
#'
#' @return plot
#'
#' @export
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
#' Plot winter-only intermittent rivers
#'
#' Utility function: plot discharge time series for GRDC or GSIM gauging stations
#' that have been deemed to monitor winter-only intermittent rivers.
#'
#' @param dt formatted data.table (generated in \code{\link{analyzemerge_gaugeir}}).
#' @param dbname (character) 'GRDC' or 'GSIM". Will determine whether to use \code{\link{plotGRDCtimeseries}}
#' or \code{\link{plotGSIMtimeseries}}
#' @param inp_resdir (character) path to directory where to save plots. A directory will be
#' created called " 'winterir_rawplots_o(yearthresh)_(date)'
#' @param yearthresh (integer) minimum year from which to plot discharge record
#' @param plotseries (logical) whether to generate plots
#'
#' @details “Winter-only” non-perennial gauging stations were defined as those
#' whose stream record contained less than one zero-flow day per year on average
#' during months with long-term mean air temperature over 10°C (averaged across
#' the local catchment immediately draining to the river reach, according to
#' WorldClim 2). In other words, “winter-only” non-perennial gauging stations
#' were those which would not have qualified as non-perennial according to our
#' criterion if only non-winter months were taken into account.
#'
#' @return subset of input dt only including winter-only intermittent irvers
#'
#' @export
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
#' Plot coastal intermittent rivers
#'
#' Utility function: plot discharge time series for GRDC or GSIM gauging stations
#' that have been deemed to be near marine waters.
#'
#' @param in_gaugep \link[sf]{sf} object of gauges with column called class99_19_rsp9_buf3k1
#' indicating whether the gauge is located within 3 km of marine waters.
#' @param dt formatted data.table (generated in \code{\link{analyzemerge_gaugeir}}).
#' @param yearthresh (integer) minimum year from which to plot discharge record
#' @param dbname (character) 'GRDC' or 'GSIM". Will determine whether to use \code{\link{plotGRDCtimeseries}}
#' or \code{\link{plotGSIMtimeseries}}
#' @param inp_resdir (character) path to directory where to save plots. A directory will be
#' created called " 'winterir_rawplots_o(yearthresh)_(date)'
#' @param plotseries (logical) whether to generate plots
#'
#' @details so-called marine stations were defined as those within 3 km of a coastline.
#'
#' @return subset of input dt only including only intermittent rivers within 3 km of a coastline
#'
#' @export
plot_coastalir <- function(in_gaugep = in_gaugep, dt = GRDCstatsdt, yearthresh,
                           dbname = 'grdc', inp_resdir = inp_resdir,
                           plotseries = TRUE) {
  #Get column name of uniquer identifiers
  idno <- fifelse(str_to_upper(dbname) == 'GRDC', 'GRDC_NO', 'gsim_no')
  
  #Select gauges within 3 km of a coastline of the given database (GRDC or GSIM)
  #and with at least 10 valid years of data
  coastalgauges <- in_gaugep[!is.na(in_gaugep$class99_19_rsp9_buf3k1) &
                               !is.na(as.data.frame(in_gaugep)[,idno]),] %>%setDT
  
  
  coastalir <- dt[get(paste0('totalYears_kept_o', yearthresh)) >= 10 &
                    get(idno) %in% coastalgauges[,get(idno)],]
  
  #Create output directory for plots
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

#------ comp_irstats -----------
#' Compute flow intermittence statistics
#'
#' Format and compute a set of summary statistics based on a yearly streamflow 
#' gauging station time series. Used in \code{\link{comp_GRDCdurfreq}}. 
#' 
#'
#' @param tabyearly data.table 
#' @param maxgap (integer) maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics.
#' @param mdurthresh (numeric) threshold of mean annual number of zero-flow days beyond
#' which to classify gauge as intermittent.
#' @param yearthresh (integer) minimum year from which to analyze discharge record.
#' @param windowsize (integer) window size to check for zero-flow days. 
#' @param fullwindow (logical) whether years for which the window is truncated 
#' (e.g., beginning and end of time series) are taken in account in moving window analysis.
#'
#' 
#' @return link[data.table]{data.table} with 1 row and 10 columns:
#' \itemize{
#'   \item firstYear_kept: the first year in the time series with a number of missing daily values <= \code{maxgap} 
#'   (i.e., the first year in the subset of the time series kept for further analysis)
#'   \item lastYear_kept: the last year in the time series with a number of missing daily values <= \code{maxgap} 
#'   (i.e., the first year in the subset of the time series kept for further analysis)
#'   \item totalYears_kept: total number of years in the time series with a number of missing daily values <= \code{maxgap} 
#'   (i.e., the number of years in the subset of the time series kept for further analysis)
#'   \item totaldays: total number of days with data in years with a number of missing daily values <= \code{maxgap} 
#'   (i.e., total number of days of data kept for further analysis)
#'   \item integerperc: proportion of daily values which are in integer format. 
#'   Only including years with <= maxgap missing daily discharge values.
#'   \item sumDur: total number of daily zero-flow values. 
#'   Only including years with <= maxgap missing daily discharge values.
#'   \item mDur: mean annual number of daily zero-flow values. 
#'   Only including years with <= maxgap missing daily discharge values.
#'   \item mFreq: mean annual frequency of zero flow events. A zero-flow event is defined 
#'   as one or more consecutive days with a recorded daily discharge of zero. 
#'   Only including years with <= maxgap missing daily discharge values.
#'    \item intermittent: binary flow intermittence class. 1: non-perennial (if mDur >= 1, i.e., 
#'    if gauging station recorded zero-flow for at least one day per year on average); 0: perennial.
#'    \item movinginter: whether there is at least one zero-flow day in every \code{windowsize}-year (e.g., 20-year) moving window across the record.
#' } 
#' 
#'
#' @export
comp_irstats <- function(tabyearly, maxgap, mdurthresh, yearthresh,
                         windowsize, fullwindow) {
  
  checkpos <- function(x) {any(x>0)}
  
  #
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
    dor_pc_pva = dor_pc_pva/10,
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
#' Threshold data (deprecated)
#'
#' Compute Balanced accuracy (BACC) after classifying the predictions from a
#' probability ml model into binary classes for every probability threshold value
#' between 0.4 and 0.6 in 0.01 increments. \cr
#'
#' Utility function used within \code{\link{formatmisclass_bm}}.
#' Not currently used in final analysis.
#'
#' @param bmres \link[mlr3]{BenchmarkResult}
#'
#' @return
#'
#'
#' @export
threshold_dat <- function(bmres) {
  #Format the benchmark result as a data.table
  bmres_dt <- as.data.table(bmres) %>%
    .[, learner_id := unlist(lapply(learner, function(x) x$id))] %>%
    .[, task_id := unlist(lapply(task, function(x) x$id))] %>%
    .[, resampling_id := unlist(lapply(resampling, function(x) x$id))] %>%
    unique(by='uhash')
  
  #If regression task, format the predictions so that they match the format
  # of a classification task.
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
  
  outtab <- lapply(seq_len(bmres_dt[,.N]), function(rsmp_i) {
    print(rsmp_i)
    rsmp_preds <- bmres$resample_result(rsmp_i)$prediction()
    
    if (bmres$task_type == 'regr') {
      rsmp_preds <- preds
    }
    
    baccthresh <- ldply(seq(0.4,0.6,0.01), function(threshold_class) {
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
    #If in_rftuned is already a ResampleResult, simply return it
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
#' @seealso \code{\link{weighted_vimportance_nestedrf}}
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
#' See \code{\link{extract_impperf_nestedrf}} for more details on computations.
#'
#' @seealso \code{\link{extract_impperf_nestedrf}} and \code{\link{ggvimp}} and
#' \code{\link{benchmark_featsel}}
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
           imp_wsd =  weighted_sd(x=importance, w=classif.bacc)), #Compute weighted sd
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
#' @seealso \code{\link{weighted_vimportance_nestedrf}},
#'  \code{\link{ggpd_bivariate}}
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
#' \dontrun{
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "in_dt", backend = in_dt, target ='intermittent')
#' get_oversamp_ratio(task)
#' }
#'
#' @export
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
#' @seealso \code{\link{create_tasks}}
#'
#' @examples
#' \dontrun{
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "task_classif", backend = in_dt, target ='intermittent')
#' convert_clastoregrtask(task, id="task_regr")
#' }
#'
#' @export
convert_clastoregrtask <- function(in_task, in_id, oversample=FALSE) {
  #Whether the minority class needs to be oversampled to match the number of observation in the majority class
  if (oversample) {
    #Identify minority class and compute ratio between the number of observations
    #in the majority and the minority classes
    oversamp_ratio <- get_oversamp_ratio(in_task)
    
    #Subset the dataset underlying the task to only keep records for the minority class
    mino_subdat <- in_task$data()[
      eval(in_task$target_names) == oversamp_ratio$minoclass,]
    
    #Get number of obs for the minority class
    nmino <- mino_subdat[,.N]
    
    #Sample minority observation indices with replacement to make up for the difference
    #in number of observations for maj and min classes
    oversamp_index <- sample(nmino, nmino*(oversamp_ratio$ratio-1), replace=T)
    
    #Actually sample records and bind them to the original data
    newdat <- rbind(in_task$data(),
                    mino_subdat[oversamp_index,]) %>%
      .[, eval(in_task$target_names):= as.numeric(as.character(get(in_task$target_names)))] #Convert target to numeric ((required for regression))
    
  } else {
    #If no oversampling, just convert target to numeric (required for regression)
    newdat <- in_task$data()[, eval(in_task$target_names) :=
                               as.numeric(as.character(get(in_task$target_names)))]
  }
  #Create regression task
  return(mlr3::TaskRegr$new(id =in_id,
                            backend = newdat,
                            target = in_task$target_names))
}

#------ durfreq_parallel -----------------
#########Not used in the end#######################
#' Parallel wrapper for \code{\link{comp_durfreq}}
#'
#' Wrapper to run \code{\link{comp_durfreq}} (which computes the annual mean number of zero
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
#' from \code{\link{comp_durfreq}}
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
#' Resample binary brier score
#'
#' Compute the binary brier score (Brier, 1950) for binary classification
#' predictions from a single \link[mlr3]{ResampleResult} extracted from an mlr3
#'  \code{\link[mlr3]{BenchmarkResult}}. \cr
#' This is a utility function.
#'
#' @param bmres \link[mlr3]{BenchmarkResult} from either a binary classification
#' model or a regression model bounded between 0 and 1.
#' @param rsmp_i (integer) index of the resample result to process.
#'
#' @return (numeric) binary brier score
#'
#' @details See the [mlr3 reference website](https://mlr3.mlr-org.com/reference/mlr_measures_classif.bbrier.html)
#'  for a definition of the Brier score for binary classification.
#'
#' @source Brier, G. W. (1950). Verification of forecasts expressed in terms of
#' probability. Monthly Weather Review, 78(1), 1\–3.
#'
#' @export
rsmp_bbrier <- function(bmres, rsmp_i) {
  rsmp <- bmres$resample_result(rsmp_i) #Get resample result
  rsmp_pred<- rsmp$prediction() #Get predictions for resample result
  
  #Compute binary brier score
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
#' Resample AUC
#'
#' Compute the Area Under the ROC Curve or AUC (Hanley and McNiel, 1982) for binary
#' classification predictions from a single \link[mlr3]{ResampleResult} extracted
#' from an mlr3 \code{\link[mlr3]{BenchmarkResult}}. \cr
#' This is a utility function.
#'
#' @inheritParams rsmp_bbrier
#'
#' @return (numeric) auc
#'
#' @details See the [mlr3 reference website](https://mlr3.mlr-org.com/reference/mlr_measures_classif.auc.html)
#'  for a definition of the AUC.
#'
#' @source Hanley, J. A., & McNeil, B. J. (1982). The meaning and use of the
#' area under a receiver operating characteristic (ROC) curve. Radiology.
#' https://doi.org/10.1148/radiology.143.1.7063747
#'
#' @export
rsmp_auc <- function(bmres, rsmp_i) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  auc <- Metrics::auc(rsmp_pred$truth, rsmp_pred$response)
  return(auc)
}

#------ rsmp_bacc ----------------
#' Resample BACC
#'
#' Compute the Balanced class ACCuracy (Brodersen, Ong, Stephan, & Buhmann, 2010)
#' from a single \link[mlr3]{ResampleResult} extracted from an mlr3
#' \code{\link[mlr3]{BenchmarkResult}}. \cr
#' This is a utility function.
#'
#' @inheritParams rsmp_bbrier
#'
#' @return (numeric) bacc
#'
#' @details See the [mlr3 reference website](https://mlr3.mlr-org.com/reference/mlr_measures_classif.bacc.html)
#'  for a definition of the BACC.
#'
#' @source Brodersen, K. H., Ong, C. S., Stephan, K. E., & Buhmann, J. M. (2010).
#' The balanced accuracy and its posterior distribution. Proceedings -
#' International Conference on Pattern Recognition, 3121–3124.
#' https://doi.org/10.1109/ICPR.2010.764
#'
#' @export
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


#------ rsmp_sen ---------
#' Resample sensitivity
#'
#' Compute the sensitivity (proportion of correctly classified non-perennial stations)
#' from a single \link[mlr3]{ResampleResult} extracted from an mlr3
#' \code{\link[mlr3]{BenchmarkResult}}. \cr
#' This is a utility function.
#'
#' @inheritParams rsmp_bbrier
#'
#' @return (numeric) sensitivity
#'
#' @details Sensitivity is also sometimes called recall. It measures the proportion o
#' of positives that are correctly identified. [Wikipedia](https://en.wikipedia.org/wiki/Sensitivity_and_specificity)
#' actually has a nice explanation for it.
#'
#' @export
rsmp_sen <- function(bmres, rsmp_i, threshold_class) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  
  if (inherits(rsmp_pred, 'PredictionClassif')) {
    prob_resp <- rsmp_pred$prob[,'1']
  }
  
  if (inherits(rsmp_pred, 'PredictionRegr')) {
    prob_resp <- rsmp_pred$response
  }
  
  response <- as.data.table(prob_resp)[
    , c(fifelse(prob_resp>=threshold_class, 1, 0))]
  
  sen <- Metrics::recall(as.numeric(as.character(rsmp_pred$truth)),
                         response)
  
  return(sen)
}

#------ rsmp_spe ---------
#' Resample specificity
#'
#' Compute the specificity (proportion of correctly classified perennial stations)
#' from a single \link[mlr3]{ResampleResult} extracted from an mlr3
#' \code{\link[mlr3]{BenchmarkResult}}. \cr
#' This is a utility function.
#'
#' @inheritParams rsmp_bbrier
#'
#' @return (numeric) specificity
#'
#' @details Specificity measures the proportion of negatives that are correctly
#' identified. [Wikipedia](https://en.wikipedia.org/wiki/Sensitivity_and_specificity)
#' actually has a nice explanation for it.
#'
#' @export
rsmp_spe <- function(bmres, rsmp_i, threshold_class) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  
  if (inherits(rsmp_pred, 'PredictionClassif')) {
    prob_resp <- rsmp_pred$prob[,'1']
  }
  
  if (inherits(rsmp_pred, 'PredictionRegr')) {
    prob_resp <- rsmp_pred$response
  }
  
  response <- as.data.table(prob_resp)[
    , c(fifelse(prob_resp>=threshold_class, 1, 0))]
  
  spe <- Metrics::recall(1-as.numeric(as.character(rsmp_pred$truth)),
                         1-response)
  return(spe)
}

#------ rsmp_pre ---------
#' Resample precision
#'
#' Compute the precision (proportion of stations classified as non-perennial that are actually IRES)
#' from a single \link[mlr3]{ResampleResult} extracted from an mlr3
#' \code{\link[mlr3]{BenchmarkResult}}. \cr
#' This is a utility function.
#'
#' @inheritParams rsmp_bbrier
#'
#' @return (numeric) precision
#'
#' @details Precision is computed as the number of true positives (i.e. the 
#' number of items correctly labelled as belonging to the positive class) divided 
#' by the total number of elements labelled as belonging to the positive class. 
#' [Wikipedia](https://en.wikipedia.org/wiki/Precision_and_recall)
#' actually has a nice explanation for it.
#'
#' @export
rsmp_pre <- function(bmres, rsmp_i, threshold_class) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  
  if (inherits(rsmp_pred, 'PredictionClassif')) {
    prob_resp <- rsmp_pred$prob[,'1']
  }
  
  if (inherits(rsmp_pred, 'PredictionRegr')) {
    prob_resp <- rsmp_pred$response
  }
  
  response <- as.data.table(prob_resp)[
    , c(fifelse(prob_resp>=threshold_class, 1, 0))]
  
  pre <- Metrics::precision(as.numeric(as.character(rsmp_pred$truth)),
                            response)
  return(pre)
}


#------ bm_paramstime ----------------
#' Benchmark result parameters and time
#'
#' Get the parameters, their values or the range of tested values through tuning +
#' the running time for each resample result in a \code{\link[mlr3]{BenchmarkResult}}.
#' \cr
#' This is a utility function.
#'
#' @param bmres \link[mlr3]{BenchmarkResult} from either a binary classification
#' model or a regression model bounded between 0 and 1.
#'
#' @return \link[data.table]{data.table} of parameter values for each resample results
#'
#' @export
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
#' Benchmark result performance measures table
#'
#' Get hyperparameters and performance measures for every \link[mlr3]{ResampleResult}
#' in a \link[mlr3]{BenchmarkResult} of classification models.
#'
#' @param bmres \link[mlr3]{BenchmarkResult} from either a binary classification
#' model or a regression model bounded between 0 and 1.
#' @param interthresh (numeric) between 0 and 1 (inclusive), threshold to compute
#' sensitivity, specificity, and balanced accuracy. If not provided, specificity
#' and sensitivity will be provided for a 0.5 threshold and the highest BACC
#' for thresholds ranging from 0.45 to 0.55 will be provided.
#'
#' @return \link[data.table]{data.table} of model hyperparameters and
#' measure performance values for each \link[mlr3]{ResampleResult}.
#'
#' @details performance measures include Binary brier score (bbrier),
#' Area Under the ROC Curve (AUC), sensitivity, specificity, and Balanced
#' classification ACCuracy (BACC).
#'
#' @seealso \code{\link{bm_paramstime}}, \code{\link{rsmp_bbrier}},
#' \code{\link{rsmp_auc}}, \code{\link{rsmp_sen}}, \code{\link{rsmp_spe}},
#' \code{\link{rsmp_bacc}}
#'
#' @export
bm_msrtab <- function(bmres, interthresh=NULL) {
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
  
  print('Getting sensitivity')
  sen_vec <- lapply(1:bmres$n_resample_results, function(i) {
    print(i)
    rsmp_sen(bmres=bmres, rsmp_i=i, threshold_class=interthresh)
  }) %>%
    unlist
  
  print('Getting specificity')
  spe_vec <- lapply(1:bmres$n_resample_results, function(i) {
    print(i)
    rsmp_spe(bmres=bmres, rsmp_i=i, threshold_class=interthresh)
  }) %>%
    unlist
  
  print('Getting precision')
  pre_vec <- lapply(1:bmres$n_resample_results, function(i) {
    print(i)
    rsmp_pre(bmres=bmres, rsmp_i=i, threshold_class=interthresh)
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
    if (is.null(interthresh)) {
      ldply(seq(0.43,0.55,0.01), function(threshold_class) {
        print(threshold_class)
        rsmp_bacc(bmres, rsmp_i=i, threshold_class=threshold_class)}) %>%
        setorder(bacc) %>%
        setDT %>%
        .[, .SD[.N, .(bacc, threshold_class)]]
    } else {
      rsmp_bacc(bmres, rsmp_i=i, threshold_class=interthresh)
    }
  }) %>%
    rbindlist
  
  print('Getting inner sampling params, general params, and time &  add to rest')
  moddt <- bm_paramstime(bmres) %>%
    cbind(., outer_smp) %>%
    cbind(., bacc_vec) %>%
    .[,`:=`(bbrier = bbrier_vec,
            auc = auc_vec,
            sen = sen_vec,
            spe = spe_vec,
            pre = pre_vec)]
  
  return(moddt)
}

#------ reset_tuning ----------------------
#' Reset AutoTuner
#'
#' Create a new \link[mlr3tuning]{AutoTuner} based on the \link[mlr3]{Learner}
#' and other hyperparameters for an input AutoTuner, but replace with a new
#' \link[mlr3]{Task}.
#'
#' @param in_autotuner \link[mlr3tuning]{AutoTuner} (single or list)
#' @param in_task new \link[mlr3]{Task} to feed into the \link[mlr3tuning]{AutoTuner}.
#' @param in_lrnid ID for the \link[mlr3]{Learner} to select if a list of \link[mlr3tuning]{AutoTuner} was provided.
#'
#' @return \link[data.table]{data.table} of model hyperparameters and
#' measure performance values for each \link[mlr3]{ResampleResult}.
#'
#' @details performance measures include Binary brier score (bbrier),
#' Area Under the ROC Curve (AUC), sensitivity, specificity, and Balanced
#' classification ACCuracy (BACC).
#'
#' @seealso \code{\link{bm_paramstime}}, \code{\link{rsmp_bbrier}},
#' \code{\link{rsmp_auc}}, \code{\link{rsmp_sen}}, \code{\link{rsmp_spe}},
#' \code{\link{rsmp_bacc}}
#'
#' @export
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
#' Format model computational metadata
#'
#' Format model Learner names for each learner in a \link[mlr3]{BenchmarkResult}. \cr
#' Utility function.
#'
#' @param bmres \link[mlr3]{BenchmarkResult}
#' @param typecomp (character) one of 'classif1', 'regr1', 'classif2'. Type of model.
#'
#' @return \link[data.table]{data.table} of BenchmarkResult with formatted model names
#'
#' @export
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
  
  return(bmres)
}

#------ compute_confustats ----------------
#' Compute confusion statistics
#'
#' Compute standard set of statistics from a confusion matrix based on classification
#' results (including misclassification rate, sensitivity, specificity, precision,
#' number of perennial and non-perennial stations correctly and incorrectly classified,
#' total number of stations, and the predicted prevalence of flow intermittence compared
#' to the observed prevalence). 
#' 
#'
#' @param in_gselpreds data.table with a row for each record and two columns, 'truth' and 'response'.
#' @param ndigits number of decimals to include in misclassification, sensitivity, 
#' specificity, and precision metrics.

#'  
#' @details Formatted statistics for Extended Data Table 2 in Messager et al. (2021).
#' 
#' @return vector of statistics
#' 
#' 
#' @export
compute_confustats <- function(in_gselpreds, ndigits=2) {
  confu <- in_gselpreds[, .N, by=.(truth, response)]
  
  outvec <- data.table(
    inter_confu =  paste0(confu[truth==1 & response==1, max(0L, N)], '|',
                          confu[truth==1 & response==0, max(0L, N)]),
    pere_confu = paste0(confu[truth==0 & response==1, max(0L, N)], '|',
                        confu[truth==0 & response==0, max(0L, N)]),
    misclas = round(confu[truth != response, sum(N)] / confu[, sum(N)], ndigits),
    sens = round(
      confu[truth == '1' & response == '1', max(0L, N)]/confu[truth=='1', sum(N)],
      ndigits),
    spec  = round(
      confu[truth=='0' & response==0, max(0L, N)]/confu[truth=='0', sum(N)],
      ndigits),
    prec = round(
      confu[truth == '1' & response == '1', max(0L, N)]/confu[response == '1', sum(N)],
      ndigits),
    N = sum(confu$N),
    predtrue_inter = paste0(round(confu[response==1, 100*sum(N)/sum(confu$N)]),
                            '|',
                            round(confu[truth==1, 100*sum(N)/sum(confu$N)]))
  )
  return(outvec)
}

#------ ggmisclass_single -----------------
#' Single ggplot of misclassification 
#'
#' Plot cross-validation sensitivity, specificity, and misclassification rate 
#' based on probability threshold for classifying a reach as non-perennial. 
#'
#' @param in_predictions Either:
#' 1. a \link[mlr3]{PredictionClassif} or
#' 2. a data.table of predictions for a set of CV repetitions as formatted by
#' \code{\link{analyze_benchmark}}.
#' @param in_rftuned Output from \link{selecttrain_rf}; 
#' list containing inner and outer resampling results + task.
#' @param spatial_rsp (boolean) whether to use spatial or non-spatial 
#' cross-validation results
#'  
#' @details in_rftuned is only needed if in_predictions is not provided.
#' 
#' @return list containing a ggplot and the threshold for which sensitivity == specificity (numeric)
#' 
#' 
#' @export
ggmisclass_single <- function(in_predictions=NULL, in_rftuned=NULL, spatial_rsp=FALSE) {
  #Get predicted probabilities of intermittency for each gauge
  # in_gaugestats[!is.na(cly_pc_cav), intermittent_predprob :=
  #                 as.data.table(in_predictions)[order(row_id), mean(prob.1), by=row_id]$V1]
  #Get misclassification error, sensitivity, and specificity for different classification thresholds
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class
  
  #If provided resampling results rather than prediction table, extract 
  if (!is.null(in_rftuned)) {
    rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
    in_predictions <- rsmp_res$prediction()
  }
  
  #Get confusion matrices for range of thresholds (i.e., probability of flow intermittence
  #above which a watercourse is classified as non-perennial)
  threshold_confu_dt <- ldply(seq(0,1,0.01), threshold_misclass, in_predictions) %>%
    setDT
  
  #Get classification threshold at which sensitivity and specificity are the most similar
  balanced_thresh <- threshold_confu_dt[which.min(abs(spec-sens)),]
  print(paste('Sensitivity =', round(balanced_thresh$sens,2),
              'and Specificity =', round(balanced_thresh$spec,2),
              'at a classification threshold of', balanced_thresh$i))
  
  #Plot trends in confusion matrix metrics with increasing threshold
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

#------ sfformat_wintri ----------------------
#' Format to sf and project to Winkel Tripel
#'
#' Convert an object (e.g. data.frame or data.table) to simple feature (spatial object)
#' and project it to [Winkel Tripel](https://en.wikipedia.org/wiki/Winkel_tripel_projection).
#'
#' @param in_sp object to be converted into an object class \code{sf}
#'
#' @return \code{[sf](https://r-spatial.github.io/sf/articles/sf1.html)} object
#'
#' @export
sfformat_wintri <- function(in_sp) {
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"
  
  return(st_as_sf(in_sp) %>%
           st_transform(crs_wintri, use_gdal = FALSE)
  )
}

#------ bin_dt -----------------
#' Bin data.table
#'
#' Bins a data.table over a numeric column.
#'
#' @param in_dt \link[data.table]{data.table} to bin.
#' @param binvar (character) column that will be used to define bins.
#' @param binfunc (character) binning approach. One of 'manual', 'equal_length', 'equal_freq'.
#' @param binarg (numeric) binning argument, depends on binning approach (\code{binfunc}).
#' @param bintrans (character or numeric) transformation of \code{binvar}, default is NULL.
#' @param ndigits (integer) number of decimals to keep for displaying formatted bin limits
#' @param na.rm (logical) whether to include NAs.
#' @param valuevar (character) if na.rm = FALSE, column to use to detect NAs and remove records.
#'
#' @return input \link[data.table]{data.table} with four new columns:
#' \itemize{
#'   \item bin - bin number (1 being the lowest value bin)
#'   \item bin_lmin - bin lower limit
#'   \item bin_lmax - bin higher limit
#'   \item bin_lformat - formatted character of bin limits
#'   (format: \code{round(bin_lmin, ndigits) - round(bin_lmax, ndigits))
#' }
#'
#' @details inspired from [rbin package](https://github.com/rsquaredacademy/rbin).
#' Differences include that it concentrates all binning approaches within a single
#' function and works on a data.table.
#'
#' binfunc: \cr
#' \itemize{
#'   \item 'manual' - bin continuous data manually. \code{binarg} sets the inner bin limits,
#'  such that the final table will have \code{length(binarg) + 1} bins. The lower end of the
#'  first bin is automatically set to be the minimum value in \code{binvar} and the upper end of
#'  the last bin is set to be the maximum value in \code{binvar}
#'
#'   \item 'equal_length' - Bin continuous data such that each bin has the same \code{binvar} interval length.
#'   If \code{bintrans} is not \code{NULL}, then interval length is computed on transformed scale.
#'   \code{binarg} (automatically rounded to the nearest integer) sets the number of bins.
#'
#'   \item 'equal_freq' - Bin continuous data such that each bin has the same number of records.
#'   \code{binarg} (automatically rounded to the nearest integer) sets the number of bins.
#' }
#'
#' bintrans: can either be 'log' (for natural log) or a numeric exponent to transform
#' according to x^bintrans.
#'
#' @seealso for examples, see applications in \code{\link{bin_misclass}},
#' \code{\link{eval_watergap}}, \code{\link{tabulate_globalsummary}},
#' \code{\link{formathistab}}, \code{\link{compare_us}}
#'
#' @export
bin_dt <- function(in_dt, binvar, binfunc, binarg,
                   bintrans=NULL, ndigits=2,
                   na.rm=FALSE, valuevar=NULL) {
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
      .[bin == i, bin_lformat := paste(round(bin_lmin, ndigits),
                                       round(bin_lmax, ndigits),
                                       sep='-')]
    
    if (i == bins) {
      in_dt[get(eval(binvar)) == u_freq[i],  bin := i]
    }
  }
  
  if (binfunc == 'equal_freq') {in_dt[, binid := NULL]}
  
  return(in_dt)
}

#------ label_manualbins ------------
#' Label manual bins
#'
#' Utility function: label bin limits for \code{\link{formathistab}}.
#'
#' @param binarg (character vector) Arguments for bin_dt manual
#' @param minval (numeric) value to set for lower limit of first bin.
#'
#' @return vector of labels
#'
#' @export
label_manualbins <- function(binarg, minval) {
  minlabel <- paste(minval, binarg[1], sep=" - ")
  otherlabels <- mapply(function(x, y) {paste(x, y-1, sep=" - ")},
                        binarg[1:(length(binarg)-1)], binarg[2:length(binarg)])
  return(c(minlabel, otherlabels))
}

#------ formathistab -----------------
#' Format histogram table
#'
#' Creates a binary frequency histogram table by computing the proportion of records
#' for which a selected column has a given value
#' (in this study, the percentage length of rivers that are deemed intermittent)
#' after binning the table by a continuous variable (e.g. river discharge).
#'
#' @inheritParams bin_dt
#' @param castvar (character) olumn that will be used to define bins —
#' the equivalent of \code{binvar} in \code{bin_dt}.
#' @param valuevar (character) variable to summarize (e.g. intermittency class).
#' @param valuevarsub (character) value of \code{valuevar} to summarize (e.g. '1' for intermittent rivers)
#' @param weightvar (character) variable to weigh proportion by. (e.g. river reach length).
#' @param binlabels (character vector) formatted labels for bins.
#' @param datname (character) optional, adds a column to output table describing what the data describes (e.g. France)
#'
#' @return \link[data.table]{data.table} with five columns:
#' \itemize{
#'   \item bin - bin number
#'   \item perc - percentage of records (or of \code{weightvar}) that meet \code{valuevar == valuevarsub}.
#'   \item binsumlength - number of records (if default \code{weightvar}) or sum of \code{weightvar} by bin (e.g. total river length).
#'   \item dat - \code{datname} argument
#'   \item binformat - label for each bin provided through \code{binlabels} argument
#' }
#'
#' @seealso used in \code{\link{compare_fr}}, \code{\link{compare_us}},
#' and \code{\link{compare_au}}
#'
#' @export
formathistab <- function(in_dt, castvar, valuevar, valuevarsub,
                         weightvar=1, binfunc, binarg, binlabels,
                         datname=NULL) {
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
#' Format number to nearest log10 bin in scientific format
#'
#' Format a number to nearest order of magnitude (log10 integer) and gives out a
#' scientific format expression.
#'
#' @param x (numeric)
#'
#' @return (expression)
#'
#' @example
#' scientific_10(45435) #expression(10^+04)
#' scientific_10(0.0000431) #expression(10^-05)
#'
#' @export
scientific_10 <- function(x) {
  parse(text=gsub(".*e", "10^", scales::scientific_format()(x)))
}
#------ ggcompare -------------------
#' Comparison ggplot histogram
#'
#' utility function: creates a ggplot histogram showing side-by-side bars
#' to compare two sets of data for the same bins. Based on \link[data.table]{data.table}
#' containing formatted histogram data.
#'
#' @param datmerge \link[data.table]{data.table}
#' @param binarg labels for inner bin limits, see \link{bin_dt}
#' @param insetx (numeric) proportion of main plot width to occupy with inset (between 0 and 1)
#' @param insety (numeric) proportion of main plot height to occupy with inset (between 0 and 1)
#'
#' @return histogram \link[ggplot2]{ggplot} object
#'
#' @details \code{datmerge} \link[data.table]{data.table} must contain all the columns
#' returned by \link{formathistab}.
#'
#' @seealso used in \code{\link{compare_fr}}, \code{\link{compare_us}},
#' and \code{\link{compare_au}}
#'
#' @export
ggcompare <- function(datmerge, binarg, insetx = 0.4, insety = 0.8) {
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
          text = element_text(size=10),
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
  
  
  xmin_inset <- datmerge[, insetx*(max(bin)-min(bin))]
  xmax_inset <- datmerge[, 1*max(bin)]
  ymin_inset <- datmerge[,round((insety * max(perc, na.rm=T) - min(perc, na.rm=T)) +
                                  min(perc, na.rm=T))]
  plot_join <- plot_inter +
    annotation_custom(grob = ggplotGrob(plot_size),
                      xmin = xmin_inset,
                      xmax = xmax_inset,
                      ymin = ymin_inset,
                      ymax = Inf)
  return(plot_join)
}

#------ test_joincount -------------------
#' Join count test
#'
#' Test auto-correlation in observed and predicted flow intermittence class among
#' gauging stations using a "join count test" for each HydroBASINS level 3.
#'
#' @param in_gauges data.table of gauging stations data, including reference and 
#' predicted flow intermittence, hydrobasin membership, and WGS84 coordinates. 
#' Formatted internally in \link{map_basinBACC}.
#'
#' @return data.table with the p-value and standard deviate of the join count test
#' for reference and predicted flow intermittence at gauging stations in each
#' HydroBASINS level 3.
#'
#' @details For each river basin that included both IRES and perennial stations and contained at least 20
#' gauging stations, we tested whether spatial predictions of intermittence differed further from a random
#' spatial distribution than the observed patterns. We did so in the following steps:
#' * We measured the degree of clustering separately for the observed and predicted 
#' flow intermittence class of gauging stations — by computing the join-count statistics
#'  (Cliff & Ord, 1981) based on four nearest neighbors (see Salima & de Bellefon, 2018 for an example implementation).
#' * We assessed whether the predicted spatial distribution of intermittence differed more from what 
#' would be expected by chance (i.e., a random distribution) than the observed distribution. This 
#' assessment was based on the standard score between the estimated join-count statistics and the joincount 
#' statistics that would be obtained based on a random spatial distribution of flow intermittence 
#' classes among the stations, using 1000 permutations.
#' 
#' The join-count statistics and permutations were computed with the spatial-cross validation predictions,
#' using the joincount.mc function from the spdep package (Bivand et al., 2009).
#'
#' @export
test_joincount <- function(in_gauges) {
  predsp_gaugep <- st_as_sf(in_gauges)
  #Convert to SpatialPointDataFrame to use in gstats and project to Goode Homolosine
  crs_aeqd <- "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  predsp_gaugep_df <- st_transform(predsp_gaugep, crs=crs_aeqd)  %>%
    as_Spatial()
  
  baslist = names(which(table(predsp_gaugep_df$PFAF_ID03) > 20))
  jclist <- lapply(baslist, function(bas) {
    print(bas)
    subdat <- subset(predsp_gaugep_df, PFAF_ID03 == bas)
    g.nb <- knn2nb(knearneigh(subdat,k=4))
    
    if ((length(unique(subdat$intermittent_o1800)) > 1) &
        all(table(subdat$intermittent_o1800)>1)) {
      jc_ref <- joincount.mc(subdat$intermittent_o1800,
                             listw2U(nb2listw(g.nb)),
                             alternative = "greater",
                             nsim = 1000)
      
      stdeviate_ref <- with(jc_ref[[2]],
                            (statistic-estimate[1])/sqrt(estimate[2]))
      pref <- jc_ref[[2]]$p.value
      
    } else {
      jc_ref <- NA
      stdeviate_ref <- NA
      pref <- NA
    }
    
    if ((length(unique(subdat$IRpredcat_CVsp)) > 1) &
        all(table(subdat$IRpredcat_CVsp)>1)) {
      jc_pred <- joincount.mc(as.factor(subdat$IRpredcat_CVsp),
                              listw2U(nb2listw(g.nb)),
                              alternative = "greater",
                              nsim = 1000)
      
      stdeviate_pred <- with(jc_pred[[2]],
                             (statistic-estimate[1])/sqrt(estimate[2]))
      ppred <- jc_pred[[2]]$p.value
    } else {
      jc_pred <- NA
      stdeviate_pred <- NA
      ppred <- NA
    }
    
    outdat <- list(
      PFAF_ID03 = bas,
      #jc_ref = jc_ref,
      #jc_pred = jc_pred,
      p_value_ref = pref,
      p_value_pred = ppred,
      stdeviate_ref = stdeviate_ref,
      stdeviate_pred = stdeviate_pred
    )
    
    return(outdat)
  }
  ) %>%
    rbindlist
  
  return(jclist)
}
##### -------------------- Workflow functions ---------------------------------
#------ def_filestructure -----------------
#' Define file structure (deprecated)
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
#' Import streamgauging station spatial data and attributes for GRDC and GSIM.
#' Add new and updated environmental predictors from RiverATLAS v1.0.9 and, optionally,
#' estimates of monthly discharge from WaterGAP v2.2.
#'
#' @param inp_GRDCgaugep path to point data for GRDC gauging stations.
#' @param inp_GSIMgaugep path to point data for GSIM gauging stations.
#' @param inp_riveratlas2 path to attribute table of RiverATLAS v1.0.9
#' @param in_monthlydischarge data frame of naturalized monthly discharge for
#' every river reach in HydroSHEDS
#'
#' @return object of class \link[sf]{sf}
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
    gaugep_outformat <- merge(gaugep_attriall, in_monthlydischarge,
                              by.x = 'HYRIV_ID', by.y = 'REACH_ID',
                              all.x=TRUE, all.y=FALSE)
  } else {
    gaugep_outformat <- gaugep_attriall
  }
  
  return(gaugep_outformat)
}

#------ read_GRDCgauged_paths -----------------
#' Read file paths to streamflow data from GRDC gauging stations
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @param inp_GRDCgaugedir path to directory containing streamflow data GRDC standard files.
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables, assuming
#' that files are called "GRDC_NO.txt", GRDC_NO being replaced with a 7-digit integer.
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
#' Read file paths to streamflow data from GSIM gauging stations
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @param inp_GSIMindicesdir path to directory containing directories for different
#' GSIm indices.
#' @param in_gaugep table containing column named \code{gsim_no} with the
#' gauge IDs that will be used to generate file path.
#' @param timestep which indices to get, 'monthly' or 'seasonal'
#'
#' @details
#' For an explanation of GSIM data structure, see [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/) and
#' the [GSIM repository](https://doi.pangaea.de/10.1594/PANGAEA.887470).
#'
#' \code{inp_GSIMindicesdir} may for example be "D://GSIM/GSIM_indices". Containing: \cr
#' D://GSIM/GSIM_indices/ \cr
#' ----  HOMOGENEITY \cr
#' ----  TIMESERIES \cr
#' --------  monthly \cr
#' --------  seasonal \cr
#' --------  yearly \cr
#' ----  README.txt
#'
#' @return vector of paths to GSIM-formatted streamflow time series indices tables
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) - Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787-804. https://doi.org/10.5194/essd-10-787-2018
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
#' @param path (character) file path to a GRDC-formatted streamflow time series table
#' @param in_gaugep (table; data.frame or data.table) of gauges' hydro-environmental attributes (including mean monthly temperature data)
#' @param maxgap (integer) maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param mdurthresh (numeric) threshold of mean annual number of zero-flow days beyond
#' which to classify gauge as intermittent.
#' @param windowsize (integer) window size to check for zero-flow days. 
#' @param fullwindow (logical) whether years for which the window is truncated 
#' (e.g., beginning and end of time series) are taken in account in moving window analysis.
#' @param monthsel (integer vector) selected months to compute the statistics over
#' @param verbose whether to print input path
#'
#' @return One row data.table with 110 columns: \cr
#' \describe{
#' \itemize{
#'   \item{GRDC_NO} - (char) unique identifier for the gauge
#'   \item{firstYear} - (num) first year on full record
#'   \item{lastYear} - (num) last year on full record
#'   \item{totalYears} - (int) total number of years on full record \cr
#'   For three subsets of the time series post-1800, post-1961, and post-1971 (e.g., suffix "mDur_o1800"):
#'   \itemize{
#'     \item{firstYear_kept} - (num) first year on record with < maxgap missing days
#'     \item{lastYear_kept} - (num) first year on record with < maxgap missing days
#'     \item{totalYears_kept} - (int) total number of years with < maxgap missing days
#'     \item{totaldays} - (num) total number of days with discharge data
#'     \item{integerperc} - (num) proportion of daily values which are in integer format. 
#'     Only including years with <= maxgap missing days.
#'     \item{sumDur} - (int) total number of days with discharge = 0
#'     \item{mDur} - (num) mean number of days/year with discharge = 0
#'     \item{mFreq} - (num) mean number of periods with discharge = 0
#'     \item{intermittent} - (factor): binary flow intermittence class. 1: non-perennial (if mDur >= 1, i.e., 
#'    if gauging station recorded zero-flow for at least one day per year on average); 0: perennial.
#'    `\item{movinginter} - (logical): whether there is at least one zero-flow day 
#'    in every \code{windowsize}-year (e.g., 20-year) moving window across the record. (with at least one day of flow between periods)
#'    `\item{monthly statistics} - (numeric) long-term mean frequency of zero-flow events 
#'    and number of zero-flow days for each month (column name example: 'Jan_mdur_o1800')
#'     \item{winteronlyir} - Indicates whether gauge is only non-perennial during winter months. 
#'     if intermittent == 1 AND the average annual number of zero-flow days during warm months < 1. 
#'     Warm months are those with mean monthly catchment air temperature >= 10 (WorldClim v2; Fick and Hijmans 2017).
#'     0: either perennial, or non-perennial outside of winter months.
#'     Only including years with <= maxgap missing daily discharge values.
#' }
#' }
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
#' Compute intermittency statistics for GSIM gauging stations
#'
#' Determine general characteristics of the whole time series and of the subset
#' of years that have less than a given threshold of missing data as well as
#' intermittency statistics. The intermittency statistics can be computed for a
#' subset of months of the year (e.g. only winter months)
#'
#' @param path_mo (character) file path to a GSIM-formatted monthly streamflow time series table (e.g., ""C:/globalIRmap/data/GSIM/GSIM_indices/TIMESERIES/monthly/MA_0000010.mon")
#' @param path_sea (character) file path to a GSIM-formatted seasonal streamflow time series table (e.g., ""C:/globalIRmap/data/GSIM/GSIM_indices/TIMESERIES/seasonal/MA_0000010.seas")
#' @param in_gaugep (table; data.frame or data.table) of gauges' hydro-environmental attributes (including mean monthly temperature data)
#' @param maxgap (integer) maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param mdurthresh (numeric) threshold of mean annual number of zero-flow days beyond
#' which to classify gauge as intermittent.
#' @param windowsize (integer) window size to check for zero-flow days. 
#' @param fullwindow (logical) whether years for which the window is truncated 
#' (e.g., beginning and end of time series) are taken in account in moving window analysis.
#' @param monthsel (integer vector) selected months to compute the statistics over
#' @param verbose whether to print input path
#'
#'
#'
#' @return One row data.table with 110 columns: \cr
#' \describe{
#' \itemize{
#'   \item{gsim_no} - (char) unique identifier for the gauge
#'   \item{firstYear} - (num) first year on full record
#'   \item{lastYear} - (num) last year on full record
#'   \item{totalYears} - (int) total number of years on full record \cr
#'   For three subsets of the time series post-1800, post-1961, and post-1971 (e.g., suffix "mDur_o1800"):
#'   \itemize{
#'     \item{firstYear_kept} - (num) first year on record with < maxgap missing days
#'     \item{lastYear_kept} - (num) first year on record with < maxgap missing days
#'     \item{totalYears_kept} - (int) total number of years with < maxgap missing days
#'     \item{totaldays} - (num) total number of days with discharge data
#'     \item{sumDur} - (int) total number of days with discharge = 0
#'     \item{mDur} - (num) mean number of days/year with discharge = 0
#'     \item{intermittent} - (factor): binary flow intermittence class. 1: non-perennial (if mDur >= 1, i.e., 
#'    if gauging station recorded zero-flow for at least one day per year on average); 0: perennial.
#'    `\item{movinginter} - (logical): whether there is at least one zero-flow day 
#'    in every \code{windowsize}-year (e.g., 20-year) moving window across the record. (with at least one day of flow between periods)
#'     \item{monthly statistics} - (numeric) long-term mean monthly number of zero-flow days (column name example: 'Jan_mdur_o1800')
#'     \item{winteronlyir} - Indicates whether gauge is only non-perennial during winter months. 
#'     if intermittent == 1 AND the average annual number of zero-flow days during warm months < 1. 
#'     Warm months are those with mean monthly catchment air temperature >= 10 (WorldClim v2; Fick and Hijmans 2017).
#'     0: either perennial, or non-perennial outside of winter months.
#'     Only including years with <= maxgap missing daily discharge values.
#' }
#' }
#' }
#'
#' @export
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
#' Plot all GRDC time series with data quality flags
#'
#' Creates pngs of streamflow time series plots of daily discharge ({m^3}/s; with flags 
#' for 0-flow values and potential outlier) for all GRDC gauges with at least 
#' 10 years of data, excluding years with more than 20 days of missing data. 
#'
#' @param in_GRDCgaugestats data.table (or list of data.tables) with time series 
#' and intermittency statistics for all GRDC gauging stations. 
#' @param yearthresh  (integer) minimum year from which to plot discharge record.
#' @param inp_resdir (character) path to the results directory in which to create folder and write output plots
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details the output graphs are written into two separate newly created directories in \code{inp_resdir} called
#' GRDCir_rawplots_\code{yearthresh}_\code{YYYYMMDD} and GRDCper_rawplots_\code{yearthresh}_\code{YYYYMMDD}
#' (e.g., GRDCir_rawplots_1800_20200512 and GRDCper_rawplots_1800_20200512). The first contains plots for non-perennial
#' gauging stations and the second contains plots for perennial gauging stations. \cr
#' \cr
#' Each plot is generated by \code{\link{plotGRDCtimeseries}} and shows
#' the time series of daily streamflow values for a station. 
#' For the flagging criteria, see documentation for \code{\link{flagGRDCoutliers}}.
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Individual points show daily discharge values (in {m^3}/s).
#'   \item blue lines link daily values (which may result in unusual patterns due to missing years).
#'   \item red points are zero-flow flow values.
#'   \item green points are non-zero flow daily values statistically flagged as potential outliers .
#'   \item black points are zero-flow values flagged as potential outliers.
#' }
#'
#' @return nothing (empty data.table)
#'
#' @export
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
#' Plot all GSIM time series with data quality flags
#'
#' Creates pngs of streamflow time series plots of monthly discharge ({m^3}/s) 
#' for all GSIM gauges with at least 10 years of data, excluding years with more 
#' than 20 days of missing data. 
#'
#' @param in_GSIMgaugestats data.table (or list of data.tables) with time series 
#' and intermittency statistics for all GSIM gauging stations. 
#' @param yearthresh  (integer) minimum year from which to plot discharge record.
#' @param inp_resdir (character) path to the results directory in which to create folder and write output plots
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details the output graphs are written into two separate newly created directories in \code{inp_resdir} called
#' GSIMir_rawplots_\code{yearthresh}_\code{YYYYMMDD} and GSIMper_rawplots_\code{yearthresh}_\code{YYYYMMDD}
#' (e.g., GSIMir_rawplots_1800_20200512 and GSIMper_rawplots_1800_20200512). The first contains plots for non-perennial
#' gauging stations and the second contains plots for perennial gauging stations. \cr
#' \cr
#' 
#' Each plot is generated by \code{\link{plotGSIMtimeseries}}.
#' Daily streamflow records from GSIM stations are unavailable. Therefore,
#' the graph shows the following:
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Blue points: mean monthly discharge
#'   \item Light blue background shading: mean ± 2SD monthly discharge
#'   \item Black points: minimum and maximum monthly discharge
#'   \item Red points show minimum monthly discharge values equal to 0
#'   \item Purple points show months for which all daily discharge values are equal to 0.
#' }
#'
#' @return nothing (empty data.table)
#'
#' @export
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


#------ get_USdata ---------------------------------
#' Compute intermittency statistics for US gauging stations
#'
#' Download original daily data from the USGS National Water Information System
#' for U.S. stations and compute intermittency statistics. 
#'
#' @param in_ids (character vector) unique USGS identifiers for gauging stations (e.g., "05540130")
#' @param maxyear (integer) maximum year to include in computation of intermittency statistics
#' @param verbose whether to print the \code{in_ids}
#'
#' @details This function is implemented in \code{\link{analyzemerge_gaugeir}} to compare
#' intermittency statistics obtained directly from USGS data compared to those obtained from
#' GSIM. This was implemented because it appears that GSIM rounded discharge data to two decimals
#' when converting cubic feet to cubic meters per second. Therefore, any discharge alue < 0.005 was considered zero.
#'
#' @return data.table with time series record size and intermittency statistics computed by 
#' \code{\link{comp_irstats}}
#'
#' @export
get_USdata <- function(in_ids, maxyear, verbose = F) {
  if (verbose) {
    print(paste0('Getting USGS data for', in_ids))
  }
  
  gaugetab <- dataRetrieval::readNWISdv(
    siteNumbers = in_ids,
    startDate = "1750-01-01",
    endDate = paste0(maxyear, "-01-01"),
    parameterCd = "00060",
    statCd = "00003"
  ) %>%
    as.data.table
  
  if ("X_..2.._00060_00003" %in% names(gaugetab)) {
    gaugetab[, X_00060_00003 := `X_..2.._00060_00003`]
  }
  
  if (nrow(gaugetab) == 0) {
    gaugetab[, `:=`(year=1800, month=1, Date="1800-01-01", X_00060_00003=0)]
  }
  
  gaugetab[,`:=`(year = as.numeric(substr(Date, 1, 4)), #Create year column
                 month = as.numeric(substr(Date, 6, 7)),
                 integervalue = fifelse(X_00060_00003 == round(X_00060_00003), 1, 0)
  )]
  
  #For each record, compute date of last non-zero flow day
  gaugetab[, prevflowdate := gaugetab[zero_lomf(X_00060_00003),'Date', with=F]] %>% #Get previous date with non-zero flow
    .[X_00060_00003 != 0, prevflowdate:=NA] #If non-zero flow, set prevflowdate to NA
  
  #Compute number of missing days per year, excluding NoData values
  gaugetab[!(X_00060_00003 %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Date)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']
  
  #Compute number of days of zero flow for years with number of gap days under threshold
  gaugetab_yearly <- merge(gaugetab[, .(missingdays=max(missingdays, na.rm=T),
                                        datadays=max(datadays, na.rm=T),
                                        integerperc=sum(integervalue,na.rm=T)
  ), by='year'],
  gaugetab[X_00060_00003 == 0, .(dur=.N,
                                 freq=length(unique(prevflowdate))), by='year'],
  by = 'year', all.x = T) %>%
    .[!is.finite(missingdays), missingdays := diny(year)] %>%
    .[!is.finite(datadays), datadays := diny(year)] %>%
    .[dur==diny(year), freq:=1] %>%
    .[is.na(dur), `:=`(dur=0, freq=0)]
  
  irstats <- comp_irstats(tabyearly = gaugetab_yearly, 
                          maxgap=20,
                          mdurthresh = 1,
                          yearthresh = 1800,
                          windowsize = 20,
                          fullwindow = F) %>%
    .[, totaldays_o1800 := as.integer(totaldays_o1800)]
  
  return(irstats)
  
}

#------ analyze_gaugeir ----------------------------
#' Quality assurance/quality checking of streamflow time series
#'
#' Record of quality-checking procedure which aimed to ensure the validity of 
#' zero-flow readings and the flow intermittence class assigned to each gauge 
#' (i.e., perennial or non-perennial).
#' 
#'
#' @param in_GRDCgaugestats list of data.tables of intermittency statistics, output from 
#' \code{\link{comp_GRDCdurfreq}} applied to all GRDC gauging stations.
#' @param in_GSIMgaugestats list of data.tables of intermittency statistics, output from 
#' \code{\link{comp_GSIMdurfreq}} applied to all GSIM gauging stations.
#' @param yearthresh (integer) minimum year from which to analyze discharge record.
#' @param in_gaugep \link[sf]{sf} object of gauging stations.
#' @param inp_resdir (character) path to the results directory in which to create folder and write output plots
#' @param plotseries (logical) whether to create plots of streamflow time series for gauging stations 
#' that are non-perennial only in winter or within 3 km from the coast.
#'
#' @details For details on the QA/QC procedure, see Supplementary Information at 
#' \link{https://www.nature.com/articles/s41586-021-03565-5}. The reason for 
#' station exclusion is also presented in an interactive online map at  
#' \link{https://messamat.github.io/globalIRmap/} in the Data and methods/Reference streamflow gauging stations tab.
#'
#' @return list of two data.tables. 
#' \itemize{
#'   \item A data.table listed as 'data' which contains the collated outputs from \code{\link{comp_GRDCdurfreq}}
#'   and\code{\link{comp_GSIMdurfreq}} for the stations that were not excluded through this QA/QC process.
#'   \item A data.table listed as flags, listing all stations excluded from further analysis and the reason their exclusion.
#' }
#'
#' @export
analyzemerge_gaugeir <- function(in_GRDCgaugestats, in_GSIMgaugestats, yearthresh,
                                 in_gaugep, inp_resdir, plotseries = FALSE) {
  ### Analyze GSIM data ########################################################
  GSIMstatsdt <- rbindlist(in_GSIMgaugestats)
  
  #------ Remove stations with unstable intermittent flow regime
  #Remove those which have at least one day per year of zero-flow day but instances
  #of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
  GSIMtoremove_unstableIR <- data.table(
    gsim_no = GSIMstatsdt[(mDur_o1800 >= 1) & (!movinginter_o1800), gsim_no],
    flag = 'removed',
    comment = 'automatic filtering: at least one no-flow day/year on average but no zero-flow event during >= 20 years of data'
  )
  
  #------ Remove stations based on examination of plots and data series
  #Outliers from examining plots of ir time series (those that were commented out were initially considered)
  GSIMtoremove_irartifacts <- list(
    c('AR_0000014', 'removed', "large gaps in data, changed flow permanence"),
    c('AT_0000021', 'removed', "single flow intermittency event, probably gap in data"),
    c('AT_0000026', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
    c('AT_0000038', 'removed', "large gaps in data, single flow intermittency event at the end"),
    c('AT_0000059', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
    c('AT_0000080', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
    c('AT_0000124', 'removed', "0 flow only at beginning of record"),
    c('BE_0000032', 'removed', 'only one zero flow event'),
    c('BR_0000286', 'removed', "Tapajos river. Impossible that it dries out."),
    c('BR_0000557', 'inspected', "Confirmed dry channel visually on satellite imagery"),
    c('BR_0000581', 'removed', "Single flow intermittency event at the end"),
    c('BR_0000593', 'removed', "short record but probably changed from perennial to IRES"),
    c('BR_0000620', 'removed', "short record but probably changed from perennial to IRES"),
    c('BR_0000645', 'to inspect', "looks flow regulated"),
    c('BR_0000651', 'to inspect', "looks flow regulated"),
    c('BR_0000662', 'removed', "flow regulated. unsure when it started, large data gap"),
    c('BR_0000664', 'inspected', "change of regime in last few years due to regulation, but originally IRES"),
    c('BR_0000706', 'removed', "too many data gaps to determine long-term flow permanence"),
    c('BR_0000717', 'removed', "flow regulated, changed from perennial to IRES"),
    c('BR_0000726', 'removed', "changed from perennial to IRES"),
    c('BR_0000778', 'removed', "too many data gaps to determine long-term flow permanence"),
    c('BR_0000786', 'removed', "only one flow intermittency event"),
    c('BR_0000862', 'removed', "only one flow intermittency event, too many data gaps"),
    c('BR_0001011', 'removed', "only one flow intermittency event, too many data gaps"),
    c('BR_0001013', 'removed', "only one flow intermittency event, too many data gaps"),
    c('BR_0001094', 'to inspect', "looks regulated"),
    c('BR_0001104', 'removed', "only one flow intermittency event at the end"),
    c('BR_0001115', 'removed', "seems to have changed flow permanence"),
    c('BR_0001116', 'removed', "seems to have changed flow permanence, many gata gaps"),
    c('BR_0001132', 'removed', "many data gaps, probably changed flow permanence"),
    c('BR_0001133', 'removed', "changed flow permanence"),
    c('BR_0001193', 'removed', "changed flow permanence"),
    c('BR_0001199', 'removed', "changed flow permanence"),
    c('BR_0001206', 'removed', "changed flow permanence"),
    c('BR_0001208', 'removed', "changed flow permanence"),
    c('BR_0001603', 'removed', "unreliable record"),
    c('BR_0002600', 'to inspect', "maybe regulated"),
    c('BR_0002666', 'to inspect', "odd patterns"),
    c('BR_0002683', 'to inspect', "odd patterns"),
    c('BR_0002687', 'removed', "only one real flow intermittency event"),
    c('BR_0003120', 'removed', "only one flow intermittency event"),
    c('CA_0000341', 'to inspect', "looks regulated"),
    c('CA_0000390', 'to inspect', "looks regulated"),
    c('CA_0000406', 'to inspect', "looks regulated"),
    c('CA_0000414', 'to inspect', "looks regulated"),
    c('CA_0000464', 'to inspect', "looks regulated"),
    c('CA_0000735', 'removed', "only one flow intermittency event since 1974"),
    c('CA_0000891', 'removed', "changed flow permanence"),
    c('CA_0000977', 'to inspect', "looks regulated"),
    c('CA_0000980', 'removed', "rounded values, only one flow intermittency event"),
    c('CA_0001025', 'to inspect', "maybe regulated"),
    c('CA_0001057', 'removed', "only one flow intermittency event in 28 years"),
    c('CA_0001130', 'removed', "either regulated or unreliable record"),
    c('CA_0001254', 'removed', "abrupt decrease to 0 flow, large gaps"),
    c('CA_0001289', 'removed', "changed flow permanence, probably regulated"),
    c('CA_0001293', 'to inspect', "looks regulated"),
    c('CA_0001319', 'removed', "abrupt decreases to 0 flow"),
    c('CA_0001386', 'removed', "abrupt decreases to 0 flow, large gap at the end"),
    c('CA_0001399', 'removed', "abrupt decreases to 0 flow, looks regulated"),
    c('CA_0001413', 'removed', "zero-flow values looks rounded"),
    c('CA_0001423', 'removed', "zero-flow only at the end"),
    c('CA_0001497', 'removed', "zero-flow only at the end"),
    c('CA_0001526', 'removed', "zero-flow are gaps, unreliable record"),
    c('CA_0001556', 'removed', "zero-flow are probably gaps, unreliable record"),
    c('CA_0001556', 'removed', "changed flow permanence, probably due to regulation"),
    c('CA_0001592', 'removed', "abrupt decreases to 0 flow, probably regulated"),
    c('CA_0001604', 'removed', "too many data gaps, difficult to assess"),
    c('CA_0001606', 'removed', "changed flow permanence, zero-flow only at the beginning of record"),
    c('CA_0001669', 'removed', "looks regulated, or rounded values"),
    c('CA_0001690', 'removed', "looks regulated"),
    c('CA_0001691', 'removed', "regulated"),
    c('CA_0001757', 'removed', "regulated"),
    c('CA_0001758', 'removed', "regulated"),
    c('CA_0001759', 'removed', "regulated"),
    c('CA_0001888', 'removed', "only one flow intermittency event"),
    c('CA_0001923', 'removed', "only one flow intermittency event at the end"),
    c('CA_0001975', 'removed', "regulated"),
    c('CA_0001989', 'removed', "regulated"),
    c('CA_0002002', 'to inspect', "maybe regulated"),
    c('CA_0002016', 'removed', "regulated"),
    c('CA_0002019', 'to inspect', "maybe regulated"),
    c('CA_0002808', 'removed', "insufficient data"),
    c('CA_0002867', 'removed', "difficult to determine flow intermittency class"),
    c('CA_0003279', 'to inspect', "changed flow permanence, looks regulated"),
    c('CA_0003290', 'to inspect', "maybe regulated"),
    c('CA_0003315', 'to inspect', "looks regulated"),
    c('CA_0003454', 'removed', "insufficient data to determine flow intermittency class"),
    c('CA_0003473', 'removed', "only one flow intermittency event"),
    c('CA_0003488', 'removed', "only one flow intermittency event"),
    c('CA_0003526', 'removed', "only one flow intermittency event over the winter"),
    c('CA_0003544', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CA_0004132', 'removed', "urban, probably changed flow permanence"),
    c('CA_0004177', 'removed', "data look rounded"),
    c('CA_0004191', 'removed', "data look rounded"),
    c('CA_0004799', 'removed', "looks regulated"),
    c('CA_0005031', 'removed', "only one flow intermittency event"),
    c('CA_0005743', 'removed', "changed flow permanence"),
    c('CN_0000002', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000004', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000009', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000010', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000012', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000013', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000020', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000021', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000022', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000026', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000029', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000032', 'removed', "no flow events do not look reliable"),
    c('CN_0000038', 'removed', "only one flow intermittency event, temporary dewatering by dam"),
    c('CN_0000043', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000047', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000062', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CN_0000063', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('CY_0000005', 'removed', "only one flow intermittency event"),
    c('CY_0000014', 'removed', "regulated"),
    c('DE_0000107', 'removed', "changed flow permanence, no flow intermittence in last 30 years"),
    c('DE_0000175', 'removed', "abrupt decreases to 0 flow"),
    c('DE_0000334', 'removed', "changed flow permanence"),
    c('DE_0000347', 'removed', "only one flow intermittence event"),
    c('DE_0000506', 'removed', "regulated"),
    c('DE_0000621', 'removed', "abrupt decrease to 0 flow, only one flow intermittence event"),
    c('ES_0000078', 'removed', "changed flow permanence, first 20 years with 0 flow"),
    c('ES_0000087', 'removed', "0 flow seems driven by exceptional drought, not throughout record"),
    c('ES_0000221', 'removed', "regulated"),
    c('ES_0000237', 'removed', "only one flow intermittence event"),
    c('ES_0000238', 'removed', "downstream of reservoir"),
    c('ES_0000249', 'removed', "only one flow intermittence event"),
    c('ES_0000354', 'removed', "only one flow intermittence event, regulated"),
    c('ES_0000362', 'removed', "regulated"),
    c('ES_0000363', 'removed', "changed flow permanence, regulated"),
    c('ES_0000364', 'removed', "changed flow permanence, regulated"),
    c('ES_0000381', 'removed', "only one flow intermittence event, regulated"),
    c('ES_0000387', 'removed', "regulated"),
    c('ES_0000388', 'removed', "too many data gaps to determine long-term flow permanence; and looks regulated"),
    c('ES_0000393', 'removed', "regulated"),
    c('ES_0000394', 'removed', "regulated"),
    c('ES_0000403', 'removed', "regulated"),
    c('ES_0000404', 'removed', "regulated"),
    c('ES_0000405', 'removed', "regulated"),
    c('ES_0000407', 'removed', "regulated"),
    c('ES_0000408', 'removed', "regulated"),
    c('ES_0000410', 'removed', "regulated"),
    c('ES_0000411', 'removed', "regulated"),
    c('ES_0000412', 'removed', "regulated"),
    c('ES_0000418', 'removed', "regulated, changed flow permanence"),
    c('ES_0000421', 'removed', "regulated"),
    c('ES_0000423', 'removed', "regulated"),
    c('ES_0000424', 'removed', "regulated"),
    c('ES_0000426', 'removed', "regulated"),
    c('ES_0000429', 'inspected', "not regulated, confirmed dry river bed"),
    c('ES_0000444', 'removed', "only one flow intermittency event"),
    c('ES_0000452', 'removed', "regulated"),
    c('ES_0000455', 'removed', "looks regulated"),
    c('ES_0000460', 'removed', "regulated, molino de chincha"),
    c('ES_0000461', 'removed', "regulated"),
    c('ES_0000525', 'removed', "very discontinued record. Seems to have changed flow permanence"),
    c('ES_0000528', 'removed', "regulated"),
    c('ES_0000529', 'removed', "regulated"),
    c('ES_0000530', 'removed', "regulated"),
    c('ES_0000542', 'removed', "regulated"),
    c('ES_0000544', 'removed', "regulated"),
    c('ES_0000581', 'removed', "appears unreliable"),
    c('ES_0000660', 'removed', "regulated, changed flow permanence"),
    c('ES_0000661', 'removed', "regulated, changed flow permanence"),
    c('ES_0000663', 'removed', "regulated, changed flow permanence"),
    c('ES_0000665', 'removed', "regulated, changed flow permanence"),
    c('ES_0000676', 'inspected', "flow regulated, but originally IRES"),
    c('ES_0000677', 'removed', "regulated"),
    c('ES_0000679', 'removed', "regulated"),
    c('ES_0000713', 'removed', "regulated, changed flow permanence"),
    c('ES_0000729', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('ES_0000733', 'removed', "regulated, changed flow permanence"),
    c('ES_0000766', 'removed', "looks regulated"),
    c('ES_0000770', 'removed', "looks regulated, impossible to tell flow permanence"),
    c('ES_0000784', 'removed', "same as ES_0000785"),
    c('ES_0000785', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('ES_0000786', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
    c('ES_0000787', 'removed', "regulated"),
    c('ES_0000794', 'removed', "same as ES_0000832"),
    c('ES_0000795', 'removed', "regulated, changed flow permanence"),
    c('ES_0000796', 'removed', "regulated"),
    c('ES_0000802', 'removed', "regulated, changed flow permanence"),
    c('ES_0000816', 'removed', "abrupt decreases to 0 flow"),
    c('ES_0000818', 'removed', "same as ES_0000785 and ES_0000784, which are not intermittent"),
    c('ES_0000830', 'removed', "regulated, changed flow permanence"),
    c('ES_0000832', 'removed', "regulated"),
    c('ES_0000841', 'removed', "regulated, changed flow permanence"),
    c('ES_0000844', 'inspected', "appears upstream of reservoir, dry bed confirmed"),
    c('ES_0000856', 'removed', "only one flow intermittency event"),
    c('ES_0000870', 'removed', "changed flow permanence, regulated"),
    c('ES_0000892', 'removed', "appears perennial with only one period of flow intermittency, data gaps"),
    c('ES_0000903', 'inspected', "not regulated, 0 flows before"),
    c('ES_0000906', 'removed', "abrupt decreases to 0 flow"),
    c('ES_0000910', 'removed', "only one flow intermittency event"),
    c('ES_0000921', 'removed', "regulated"),
    c('ES_0000931', 'inspected', "not super reliable but multiple instances of zero-flow, not regulated"),
    c('ES_0000933', 'removed', "only two flow intermittency events"),
    c('ES_0000951', 'removed', "only one flow intermittency event"),
    c('ES_0000956', 'removed', "regulated, changed flow permanence"),
    c('ES_0000958', 'removed', "only one flow intermittency event, abrupt decrease"),
    c('ES_0000959', 'removed', "only one flow intermittency event, abrupt decrease"),
    c('ES_0000962', 'removed', "only few flow intermittency event, abrupt decrease"),
    c('ES_0000973', 'removed', "only flow intermittency event at the end"),
    c('ES_0000974', 'removed', "only one flow intermittency event"),
    c('ES_0000976', 'removed', "only one flow intermittency event"),
    c('ES_0000977', 'removed', "only one flow intermittency event"),
    c('ES_0000986', 'removed', "only one flow intermittency event, abrupt decrease"),
    c('ES_0000996', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('ES_0001002', 'removed', "only one flow intermittency event"),
    c('ES_0001004', 'removed', "only one flow intermittency event"),
    c('ES_0001005', 'removed', "only one flow intermittency event"),
    c('ES_0001020', 'removed', "too many data gaps to determine long-term flow permanence"),
    c('ES_0001052', 'removed', "regulated, changed flow permanence"),
    c('ES_0001056', 'removed', "regulated, changed flow permanence"),
    c('ES_0001079', 'removed', "looks regulated"),
    c('ES_0001082', 'removed', "regulated, changed flow permanence"),
    c('ES_0001085', 'removed', "only one flow intermittency event at the end"),
    c('ES_0001105', 'removed', "regulated"),
    c('ES_0001116', 'removed', "only one year of flow intermittency events at the end"),
    c('ES_0001162', 'removed', "abrupt decreases to 0 flow"),
    c('FI_0000015', 'removed', "abrupt decrease to 0 flow"),
    c('FI_0000102', 'removed', "regulated"),
    c('FI_0000104', 'removed', "regulated"),
    c('FI_0000107', 'removed', "lake inlet, maybe just become standing water when high water level"),
    c('FI_0000119', 'removed', "regulated"),
    c('FI_0000156', 'removed', "no 0 flow for first 60% of record"),
    c('FR_0000052', 'removed', "no 0 flow for first 20 years of record"),
    c('FR_0000389', 'removed', "only one flow intermittency event"),
    c('FR_0000455', 'removed', "abrupt decrease to 0"),
    c('FR_0000486', 'removed', "abrupt decrease to 0, probably data gap"),
    c('FR_0000489', 'removed', "abrupt decrease to 0, probably data gap, only one event"),
    c('FR_0000622', 'removed', "only one flow intermittency event"),
    c('FR_0000730', 'removed', "only one flow intermittency period in 42 years"),
    c('FR_0000622', 'removed', "only one flow intermittency event, too short of a period to tell"),
    c('FR_0000820', 'removed', "only one flow intermittency event"),
    c('FR_0000821', 'removed', "only one flow intermittency event"),
    c('FR_0001240', 'removed', "perennial for first 22 years"),
    c('GB_0000196', 'removed', "only one flow intermittency event"),
    c('HU_0000017', 'removed', "only one 0 flow occurence in 25 years"),
    c('IE_0000014', 'removed', "only one flow intermittency event"),
    c('IN_0000014', 'removed', "no 0 flow occurrence in first 25 years"),
    c('IN_0000022', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000023', 'removed', "no 0 flow occurrence in first 20 years"),
    c('IN_0000024', 'removed', "only two zero-flow event in 30 years"),
    c('IN_0000032', 'removed', "only two zero-flow event in 30 years"),
    c('IN_0000045', 'removed', "no 0 flow occurrence in first 20 years"),
    c('IN_0000046', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000050', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000062', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000063', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000064', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000074', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000075', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000085', 'removed', "changed flow permanence, all tributaries are regulated"),
    c('IN_0000094', 'removed', "only one zerof-low event in last 18 years"),
    c('IN_0000105', 'removed', "regulated -- no records pre-1968 time of dam building"),
    c('IN_0000113', 'removed', "only one 0 flow occurence"),
    c('IN_0000117', 'removed', "looks regulated and changed flow permanence"),
    c('IN_0000121', 'removed', "only one 0 flow occurence"),
    c('IN_0000122', 'removed', "changed flow permanence"),
    c('IN_0000124', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000125', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000127', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000134', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000136', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000142', 'removed', "only one 0 flow occurence in first 15 years"),
    c('IN_0000159', 'removed', "regulated -- no records pre-1968 time of dam building"),
    c('IN_0000168', 'removed', "only one 0 flow occurence in first 20 years"),
    c('IN_0000170', 'removed', "only one 0 flow occurence"),
    c('IN_0000172', 'removed', "only one 0 flow occurence"),
    c('IN_0000174', 'removed', "changed flow permanence"),
    c('IN_0000190', 'removed', "abrupt decreases to 0"),
    c('IN_0000198', 'removed', "only one 0 flow occurence"),
    c('IN_0000202', 'removed', "only one 0 flow occurence"),
    c('IN_0000215', 'inspected', "doesn't seem voerly regulated"),
    c('IN_0000247', 'removed', "only one 0 flow occurrence in first 20 years"),
    c('IN_0000249', 'removed', "changed flow permanence"),
    c('IN_0000250', 'removed', "changed flow permanence"),
    c('IN_0000255', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000257', 'removed', "only one 0 flow occurence"),
    c('IN_0000283', 'removed', "changed flow permanence"),
    c('IN_0000105', 'removed', "regulated -- no records pre-1968 time of dam building"),
    c('IN_0000280', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000309', 'inspected', "before building of reservoir in 1988"),
    c('IN_0000312', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000313', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000315', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('IN_0000317', 'removed', "only one 0 flow occurence"),
    c('IT_0000161', 'removed', "only two 0 flow occurence"),
    c('JM_0000004', 'removed', "all integers"),
    c('MA_0000002', 'inspected', 'confirmed dry bed with imagery'),
    c('MX_0000032', 'removed', 'abrupt decreases to 0'),
    c('MZ_0000010', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('NA_0000050', 'removed', "unreliable record, data gaps and interpolated values"),
    c('NO_0000004', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('NO_0000018', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('NO_0000020', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('NO_0000024', 'removed', "only one 0 flow occurence"),
    c('NO_0000028', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('NO_0000029', 'removed', "errorneous 0 flow occurrence"),
    c('NO_0000030', 'removed', "only one 0 flow occurence"),
    c('NO_0000044', 'removed', "0 flow values only at the beginning, probably data gaps"),
    c('NO_0000060', 'removed', "only one 0 flow occurence"),
    c('NO_0000090', 'removed', "0 flow values at the beginning probably data gaps, otherwise only one 0 flow event"),
    c('NO_0000107', 'removed', "only one 0 flow occurence"),
    c('NO_0000134', 'removed', "only one 0 flow occurence"),
    c('OM_0000001', 'inspected', "downstream of reservoir but confirmed dry channels upstrea,"),
    c('RU_0000024', 'removed', "0 flow probably due to integer records"),
    c('RU_0000026', 'removed', "0 flow probably due to integer records"),
    c('RU_0000062', 'removed', "0 flow probably due to integer records"),
    c('RU_0000089', 'removed', "0 flow probably due to integer records"),
    c('RU_0000093', 'removed', "0 flow probably due to integer records"),
    c('RU_0000136', 'removed', "no zero flow"),
    c('RU_0000189', 'removed', "only one 0 flow occurence"),
    c('RU_0000250', 'removed', "only two 0 flow occurence"),
    c('RU_0000265', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('RU_0000269', 'removed', "abrupt decreases to zero flow"),
    c('RU_0000278', 'removed', "insufficient record to tell flow intermittency class"),
    c('RU_0000350', 'removed', "only two 0 flow events in 50 years"),
    c('RU_0000358', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('RU_0000361', 'removed', "only two 0 flow events in 50 years"),
    c('RU_0000363', 'removed', "insufficient record to tell flow intermittency class"),
    c('RU_0000374', 'removed', "only one 0 flow occurence"),
    c('RU_0000441', 'removed', "only one 0 flow occurence"),
    c('SE_0000053', 'removed', "only one 0 flow occurrence in last 30 years"),
    c('SE_0000058', 'removed', "downstream of dam, changed flow permanence. and experiences ice"),
    c('SG_0000001', 'removed', "0 flow probably due to interger records"),
    c('TH_0000034', 'removed', "0 values seem abrupt"),
    c('TZ_0000018', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('UA_0000024', 'removed', "only one 0 flow occurence"),
    c('UA_0000044', 'removed', "only one 0 flow occurence"),
    c('UA_0000064', 'removed', "only one 0 flow occurence"),
    c('US_0000056', 'removed', "regulated"),
    c('US_0000071', 'removed', "regulated"),
    c('US_0000492', 'removed', "regulated"),
    c('US_0000526', 'removed', "regulated"),
    c('US_0000761', 'removed', "abrupt decreases to 0, regulated"),
    c('US_0001665', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('US_0001856', 'removed', "only one 0 flow occurence"),
    c('US_0002247', 'removed', "regulated"),
    c('US_0002248', 'removed', "regulated"),
    c('US_0002720', 'removed', "onle two 0 flow events in 42 years"),
    c('US_0002791', 'removed', "looks regulated"),
    c('US_0003775', 'removed', "record is not long enough to determine whether flow intermittency is outlier"),
    c('US_0003835', 'removed', "regulated"),
    c('US_0003836', 'removed', "regulated"),
    c('US_0003837', 'removed', "regulated"),
    c('US_0004668', 'removed', "regulated"),
    c('US_0004711', 'removed', "only one 0 flow occurrence"),
    c('US_0004738', 'removed', "only one 0 flow occurrence"),
    c('US_0004783', 'removed', "changed from non-perennial to perennial"),
    c('US_0004882', 'removed', "changed from perennial to non-perennial"),
    c('US_0005054', 'removed', "changed from perennial to non-perennial"),
    c('US_0005064', 'removed', "probably regulated, changed from perennial to non-perennial"),
    c('US_0005069', 'removed', "probably regulated, looks heavily altered"),
    c('US_0005090', 'removed', "probably regulated, changed from perennial to non-perennial"),
    c('US_0005103', 'removed', "regulated"),
    c('US_0005104', 'removed', "probably regulated, changed from perennial to non-perennial"),
    c('US_0005105', 'removed', "regulated"),
    c('US_0005106', 'removed', "probably changed flow permanence from perennial to non-perennial"),
    c('US_0005123', 'removed', "probably regulated, changed from perennial to non-perennial"),
    c('US_0005559', 'removed', "regulated"),
    c('US_0005583', 'removed', "regulated"),
    c('US_0005709', 'inspeced', "confirmed dry"),
    c('US_0005874', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('US_0005876', 'removed', "probably changed flow permanence from perennial to non-perennial"),
    #'US_0001855', #maybe rounded values --- checked
    #'US_0001861', #maybe rounded values --- checked
    #'US_0001868', #maybe rounded avlues --- checked
    c('US_0002247', 'removed', "regulated, GRDC 4149415 upstream not intermittent"),
    c('US_0002248', 'removed', "regulated, just downstream of US_0002247"),
    c('US_0002791', 'removed', "on usgs website: close proximity to the Ohio River. During periods of high water on the Ohio River, the computed discharge at this site may be incorrectly displayed due to the backwater effect created."),
    #'US_0003591', #maybe rounded values --- checked
    # 'US_0003774', #maybe rounded values --- checked
    # 'US_0003836', #maybe rounded values --- checked
    # 'US_0004023', #maybe rounded values --- checked
    # 'US_0004216', #maybe rounded values --- checked
    # 'US_0004232', #maybe rounded values --- checked
    # 'US_0004658', #maybe rounded values --- checked
    c('US_0004773', 'removed', "regulated, no records prior to reservoir buiding"),
    c('US_005099', 'inspected', "confirmed dry bed on satellite imagery"),
    # 'US_0005161', #maybe rounded values --- checked
    # 'US_0005177', #maybe rounded values --- checked
    c('US_0005303', 'inspected', "looks fine on usgs website and imagery. just small"),
    c('US_0005596', 'inspected', "looks fine on usgs website and imagery. just small"),
    c('US_0005597', 'inspected', "looks fine on usgs website and imagery. just small"),
    # 'US_0005622', #maybe rounded values --- checked
    # 'US_0005623', #maybe rounded values --- checked
    # 'US_0005684', #maybe rounded values --- checked
    # 'US_0005687', #maybe rounded values --- checked
    c('US_0005732', 'removed', "regulated by Lake Arcadia reservoir, changed flow permanence"),
    #'US_0005859', #maybe rounded values --- checked
    #'US_0005879', #maybe rounded values --- checked
    # 'US_0006073', #maybe rounded values --- checked
    c('US_0006103', 'removed', "regulated, changed flow permanence"),
    # 'US_0006109', #maybe rounded values --- check
    # 'US_0006154', #maybe rounded values --- check
    c('US_0006105', 'removed', "regulated"),
    c('US_0006155', 'removed', "regulated, did not change flow permanence but same as 0006156"),
    c('US_0006156', 'removed', "regulated"),
    c('US_0006206', 'inspected', "looks fine on usgs website and imagery"),
    #'US_0006301', #maybe rounded values --- checked
    #'US_0006327', #maybe rounded values --- checked
    c('US_0006396', 'removed', "probably changed from perennial to non-perennial"),
    #'US_0006387', #maybe rounded values --- checked
    c('US_0006396', 'removed', "erroneous values pre-1960"),
    c('US_0006440', 'removed', "probably changed flow permanence from perennial to non-perennial"),
    c('US_0006537', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('US_0006575', 'removed', "only one 0 flow event"),
    c('US_0007055', 'removed', "regulated"),
    #'US_0006975', #maybe rounded values --- checked
    #'US_0006984', #maybe rounded values --- checked
    #'US_0006985', #maybe rounded values --- checked
    #'US_0006986', #maybe rounded values --- checked
    c('US_0008607', 'removed', "regulated since 1978, schanged flow permanence from perennial to non-perennial"),
    c('US_0008642', 'removed', "changed flow permanence from perennial to non-perennial"),
    c('US_0008687', 'removed', "regulated since 1935, schanged flow permanence from perennial to non-perennial"),
    # 'US_0008726', #maybe rounded values --- checked
    # 'US_0008779', #maybe rounded values --- checked
    c('ZA_0000007', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
    c('ZA_0000008', 'removed', "0 flow record are gaps in data"),
    c('ZA_0000047', 'removed', "too many gaps in data"),
    c('ZA_0000067', 'removed', "only two 0 flow events"),
    c('ZA_0000073', 'removed', "0 flow records are gaps in data"),
    c('ZA_0000074', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
    c('ZA_0000146', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
    c('ZA_0000159', 'removed', "abrupt decrease to 0 flow, regulated"),
    c('ZA_0000164', 'removed', "abrupt decrease to 0 flow"),
    c('ZA_0000268', 'removed', "probably changed flow permanence"),
    c('ZA_0000294', 'removed', "abrupt decrease to 0 flow"),
    c('ZA_0000084', 'removed', "probably changed flow permanence"),
    c('ZA_0000270', 'removed', "changed flow permanence"),
    c('ZW_0000009', 'removed', "regulated"),
    c('ZW_0000052', 'removed', "regulated"),
    c('ZW_0000075', 'removed', "regulated")
  ) %>%
    do.call(rbind, .) %>%
    as.data.table %>%
    setnames(c('gsim_no', 'flag', 'comment'))
  
  
  GSIMtoremove_perartifacts <- list(
    c('AR_0000023', 'removed', 'regulated'),
    c('BR_0000809', 'removed', 'regulated'),
    c('BR_0002561', 'removed', 'regulated'),
    c('BR_0003275', 'removed', 'regulated'),
    c('CA_0000215', 'inspected', 'natural lake outlet'),
    c('CA_0000665', 'removed', 'regulated, urban'),
    c('CA_0000765', 'removed', 'regulated, urban'),
    c('CA_0000838', 'removed', 'regulated, urban'),
    c('CA_0000844', 'removed', 'regulated, urban'),
    c('CA_0004040', 'removed', 'regulated'),
    c('CA_0004808', 'removed', 'regulated'),
    c('ES_0000899', 'removed', 'regulated'),
    c('IN_0000087', 'removed', 'appears downstream of dam, unsure'),
    c('IN_0000276', 'removed', 'both tributaries are regulated'),
    c('NO_0000061', 'removed', 'regulated'),
    c('RU_0000263', 'removed', 'regulated'),
    c('RU_0000264', 'removed', 'regulated'),
    c('RU_0000297', 'removed', 'regulated'),
    c('RU_0000360', 'removed', 'reservoirs on two of the main tributaries'),
    c('US_0000083', 'removed', 'regulated'),
    c('US_0000085', 'removed', 'regulated'),
    c('US_0000150', 'removed', 'regulated, urban'),
    c('US_0005594', 'removed', 'regulated, urban'),
    c('US_0005866', 'removed', 'reservoirs on headwaters of every tributary totalling 1/3 of watershed'),
    c('US_0006172', 'removed', 'regulated, urban'),
    c('US_0006173', 'removed', 'regulated, urban'),
    c('US_0006174', 'removed', 'regulated, urban'),
    c('US_0006175', 'removed', 'regulated, urban'),
    c('US_0006177', 'removed', 'regulated, urban'),
    c('US_0006180', 'removed', 'regulated, urban'),
    c('US_0006181', 'removed', 'regulated, urban'),
    c('US_0006182', 'removed', 'regulated, urban'),
    c('US_0006186', 'removed', 'regulated, urban'),
    c('US_0006187', 'removed', 'regulated, urban'),
    c('US_0006328', 'removed', 'regulated, urban'),
    c('US_0006454', 'removed', 'changed flow permanence, reservoir upstream')
  ) %>%
    do.call(rbind, .) %>%
    as.data.table %>%
    setnames(c('gsim_no', 'flag', 'comment'))
  
  
  
  #-----  Check that US stations actually have 0 flow events -------
  USstats <- as.data.table(in_gaugep)[
    substr(gsim_no, 1, 2) == 'US' & 
      !(gsim_no %in% c(GSIMtoremove_irartifacts$gsim_no, 
                       GSIMtoremove_unstableIR$gsim_no)),
    cbind(gsim_no,
          get_USdata(
            in_ids = fifelse(nchar(reference_no)==8, 
                             reference_no,
                             paste0(0, reference_no)
            ),
            maxyear = max(GSIMstatsdt$lastYear_kept_o1800) + 1
          )
    ), by=reference_no]
  
  GSIMstatscheck_US <- merge(
    GSIMstatsdt[!(gsim_no %in% c(GSIMtoremove_irartifacts$gsim_no, 
                                 GSIMtoremove_unstableIR$gsim_no)
    ),],
    USstats, by='gsim_no', all.x=F, all.y=F)
  
  GSIMtoremove_USroundingfuckup <- GSIMstatscheck_US[
    mDur_o1800.x >= 1 & mDur_o1800.y < 1,
    list(gsim_no = gsim_no,
         flag = 'removed',
         comment = 'automatic filtering: incorrect rounding of flow values during conversion by GSIM'
    )
  ]

  #-----  Check flags in winter IR for GSIM
  wintergaugesall_GSIM <- plot_winterir(
    dt = GSIMstatsdt, dbname = 'gsim', inp_resdir = inp_resdir,
    yearthresh = 1800, plotseries = plotseries)
  #
  # #Check suspicious canadian ones
  # GSIMwintermeta <- in_gaugep[in_gaugep$gsim_no %in% wintergaugesall_GSIM$gsim_no,]
  # canadians_toinspect <- in_gaugep[in_gaugep$gsim_no %in%
  #                                    paste0('CA_000', c(3469, 3473, 3526, 3544, 6082, 6122)),]$reference_no
  #
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
  #
  #Check others 'CN_0000047', 'NO_0000018', 'RU_0000089',
  #'RU_0000391', 'RU_0000393', 'RU_00000395', 'RU_0000436',
  #'RU_0000470', 'US_0008687')
  #check <- readformatGSIMmon(GSIMstatsdt[gsim_no == 'US_0008687',path])
  
  #Remove
  GSIMtoremove_winterIR <- list(
    c('CA_0003473', 'removed', "sudden peak — unsure about estimated discharge under ice conditions"),
    c('CN_0000047', 'removed', "Anomalous change from near 0 discharge to 150 m3/s, no explanation"),
    c('RU_0000391', 'removed', "Stopped recording during the winter the last ~10 years. Maybe questionable winter data"),
    c('RU_0000393', 'removed', "Didn't record during the winter for the first 20 years. Maybe questionable winter data"),
    c('RU_0000436', 'removed', "Didn't record during the winter for the first 20 years. Maybe questionable winter data"),
    c('RU_0000470', 'removed', "Didn't record during the winter for the first 20 years. Maybe questionable winter data")
  ) %>%
    do.call(rbind, .) %>%
    as.data.table %>%
    setnames(c('gsim_no', 'flag', 'comment'))
  
  #-----  Check flags in coastal IR for GSIM
  GSImcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GSIMstatsdt,
                                     dbname = 'gsim', inp_resdir = inp_resdir,
                                     yearthresh = 1800, plotseries = plotseries)
  #Already removed suspect ones
  
  ### Analyze GRDC data ########################################################################
  GRDCstatsdt <- rbindlist(in_GRDCgaugestats)
  
  #Remove all gauges with 0 values that have at least 99% of integer values as not reliable (see GRDC_NO 6140700 as example)
  GRDCtoremove_allinteger <- data.table(
    GRDC_NO = GRDCstatsdt[integerperc_o1800 >= 0.95 &
                            intermittent_o1800 == 1, GRDC_NO],
    flag = 'removed',
    comment = 'All integer discharge values'
  )
  
  #Remove those which have at least one day per year of zero-flow day but instances
  #of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
  GRDCtoremove_unstableIR <- data.table(
    GRDC_NO = GRDCstatsdt[(mDur_o1800 >= 1) & (!movinginter_o1800) &
                            !(GRDC_NO %in% c(1160115, 1160245, 4146400)), GRDC_NO],
    flag = 'removed',
    comment = 'automatic filtering: at least one no-flow day/year on average but no zero-flow event during >= 20 years'
  )
  
  #Outliers from examining plots of ir time series (those that were commented out were initially considered)
  GRDCtoremove_irartifacts <- list(
    c(1104800, 'inspected', 'low-quality gauging but confirmed seasonally intermittent flow'),
    c(1134300, 'removed', 'changed flow permanence from perennial to non-perennial, large data gaps'),
    c(1134500, 'removed', 'only 1 occurrence of 0 flow values'),
    c(1159110, 'inspected', 'regulated but drying was confirmed upstream of reservoir'),
    c(1159120, 'inspected', 'at weird or dam; but drying was confirmed upstream of reservoir'),
    c(1159132, 'inspected', 'regulated but drying was confirmed upstream of dam'),
    c(1159302, 'removed', 'abrupt decreases to 0 flow values'),
    c(1159303, 'removed', 'unreliable record, isolated 0s, sudden jumps and capped at 77'),
    c(1159320, 'inspected', "0 values for the first 14 years but still apparently originally IRES"),
    c(1159325, 'removed', "0 values for most record. probably episodic and due to series of agricultural ponds"),
    c(1159510, 'removed', "0 values for most record, on same segment as 1159511 but seems unreliable"),
    c(1159520, 'inspected', "values seem capped after 1968, otherwise seem fine. Could just be rating curve"),
    c(1159830, 'removed', 'only one occurence of 0 flow values'),
    c(1160101, 'removed', 'abrupt decreases to 0 flow values'),
    c(1160210, 'inspected', 'lower plateaus in the beginning of record are likely 0-flow values'),
    c(1160245, 'removed', 'regulated, at reservoir'),
    c(1160301, 'inspected', 'lower plateaus in the beginning of record are likely 0-flow values'),
    c(1160340, 'removed', 'abrupt decreases to 0 flow values'),
    c(1160378, 'inspected', 'appears IRES before regulation, reservoir fully dry on satellite imagery, so keep as intermittent even when not regulated'),
    c(1160420, 'removed', 'decrease to 0 appears a bit abrupt'),
    c(1160435, 'removed', 'unreliable record, abrupt decreases to 0, capped'),
    c(1160470, 'removed', 'unreliable record, probably change of rating curve in 1947, mostly missing data until 1980 but truly intermittent based on imagery'),
    c(1160540, 'inspected', 'only 0 - nodata for first 15 years. seemingly good data post 1979 and IRES'),
    c(1160635, 'removed', 'valid 0 values only during early 80s'),
    c(1160670, 'removed', 'regulated'),
    c(1160675, 'removed', 'some outliers but otherwise most 0 values seem believable'),
    c(1160780, 'removed', 'unreliable record, abrupt decreases to 0 flow values, large data gaps'),
    c(1160775, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1160785, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1160793, 'removed', 'changed flow permanence from perennial to non-perennial, unreliable record'),
    c(1160795, 'removed', 'abrupt decreases to 0 flow values'),
    c(1160800, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1160840, 'removed', 'only 2 zero flow values are believable, others are outliers'),
    c(1160850, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1160880, 'removed', 'unreliable record. Tugela river, perennial'),
    c(1160881, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1160900, 'removed', 'most 0 values look like outliers, abrupt decreases'),
    c(1160911, 'removed', 'most 0 values look like outliers, abrupt decreases'),
    c(1160971, 'removed', 'most 0 values look like outliers, abrupt decreases'),
    c(1160975, 'removed', 'most 0 values look like outliers, abrupt decreases'),
    c(1196102, 'removed', 'unreliable record, large data gaps, hard to tell original flow permanence'),
    c(1196141, 'removed', "doesn't look reliable, hard to assess long term flow permanence"),
    c(1196160, 'removed', 'changed flow permanence, some outlying 0 flow values but most are good'),
    c(1197500, 'removed', 'only one flow intermittency event, abrupt decrease to 0'),
    c(1197540, 'removed', 'abrupt decreases to 0'),
    c(1197591, 'removed', 'abrupt decreases to 0'),
    c(1197700, 'removed', 'abrupt decreases to 0'),
    c(1197740, 'removed', 'some outlying 0 flow values but most are good'),
    c(1199100, 'removed', 'most 0 values look like outliers'),
    c(1199200, 'removed', 'abrupt decreases to 0'),
    c(1199410, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1234130, 'inspected', 'low quality record but confirmed intermittent by https://doi.org/10.3390/w11010156'),
    c(1259500, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1259800, 'removed', '0 values come from integer-based part of the record'),
    c(1286690, 'removed', 'changed flow permanence, record too short to determine original flow permanence'),
    c(1289230, 'removed', 'unreliable record, gaps, shifts'),
    c(1259800, 'removed', 'changed flow permanence, only one 0 flow value post 1963'),
    c(1428400, 'removed', '0 values come from integer-based part of the record'),
    c(1428500, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1434200, 'removed', 'almost all integers'),
    c(1434300, 'removed', 'almost all integers'),
    c(1434810, 'removed', '0 values come from integer-based part of the record'),
    c(1491790, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1491815, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1491870, 'removed', 'some outlying 0 flow values but most are good'),
    c(1494100, 'removed', 'abrupt decreases to 0'),
    c(1494100, 'inspected', 'regulated but naturally intermittent'),
    c(1495360, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1495700, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(1495720, 'removed', 'too many gaps, probably changed flow permanence from perennial to non-perennial'),
    c(1591110, 'removed', "doesn't look reliable, changed flow permanence from perennial to non-perennial"),
    c(1591730, 'removed', 'abrupt decreases to 0'),
    c(1733600, 'removed', '0 values come from integer-based part of the record and outliers'),
    c(1837410, 'removed', 'abrupt decreases to 0'),
    c(1837430, 'inspected', 'nearly same as 1837410. Naturally intermittent before dam'),
    c(1897550, 'removed', 'abrupt decreases to 0'),
    c(1898501, 'removed', 'abrupt decreases to 0'),
    c(1992840, 'removed', 'changed flow permanence'),
    c(1992400, 'removed', 'most 0 values look like outliers'),
    c(2181960, 'removed', 'series of small reservoirs upstream on both tributaries'),
    c(2588500, 'removed', 'abrupt decreases to 0'),
    c(2588551, 'removed', 'abrupt decreases to 0'),
    c(2588630, 'removed', 'abrupt decreases to 0'),
    c(2588640, 'removed', 'abrupt decreases to 0'),
    c(2588708, 'removed', 'abrupt decreases to 0'),
    c(2588820, 'removed', 'abrupt decreases to 0'),
    c(2589230, 'removed', 'abrupt decreases to 0'),
    c(2589370, 'removed', 'abrupt decreases to 0'),
    c(2591801, 'removed', 'abrupt decreases to 0'),
    c(2694450, 'removed', 'abrupt decreases to 0'),
    c(2969081, 'removed', 'abrupt decreases to 0'),
    c(2999920, 'removed', '0 values come from integer-based part of the record'),
    c(3650380, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(3650460, 'removed', 'large portion of integer values, probably changed flow permanence'),
    c(3650470, 'removed', '0 values come from integer-based part of the record'),
    c(3650475, 'removed', 'large portion of integer values, probably changed flow permanence'),
    c(3650610, 'removed', 'integers pre-1960s but still intermittent after'),
    c(3650640, 'removed', 'abrupt decreases to 0'),
    c(3650649, 'inspected', 'change of flow regime due to dam building but intermittent before'),
    c(3650690, 'inspected', 'abrupt decreases to 0'),
    c(3650860, 'removed', 'only one 0 flow event in first 28 years'),
    c(3650928, 'removed', 'most 0 values look like outliers'),
    c(3652050, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(3652135, 'removed', 'only one valid 0-flow event'),
    c(3652200, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(3844460, 'removed', 'abrupt decreases to 0'),
    c(3844460, 'removed', 'abrupt decreases to 0'),
    c(4101451, 'inspected', 'station downstream also has 0s'),
    c(4103700, 'removed', '0 values come from integer-based part of the record'),
    c(4150605, 'inspected', 'just downstream of lwesville dam in Dallas. previously intermittent as well but will be removed anyways as >50% dor'),
    c(4151513, 'inspected', 'looks regulated but will be removed as > 50% regulated'),
    c(4185051, 'removed', 'regulated'),
    c(4203820, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4208043, 'removed', 'changed from perennial to non-perennial'),
    c(4208195, 'removed', 'unreliable record, 0 flow values stem from interpolation'),
    c(4208372, 'removed', 'abrupt decreases to 0 flow, probably data gaps'),
    c(4208585, 'removed', 'changed from perennial to non-perennial'),
    c(4208655, 'removed', 'insufficient data to tell flow permanence'),
    c(4208855, 'removed', 'insufficient data to tell flow permanence'),
    c(4208857, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4213090, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4213091, 'removed', 'only 2 occurrence of 0 flow values'),
    c(4213566, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4213905, 'removed', "regulated, changed flow permanence"),
    c(4214075, 'removed', '0 flow values are data gaps'),
    c(4214200, 'removed', 'changed from perennial to non-perennial'),
    c(4214297, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4214298, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4234300, 'removed', "regulated, changed flow permanence"),
    c(4243610, 'removed', "regulated, abrupt decrease to 0 probably due to reservoir filling/construction"),
    c(4351710, 'removed', '0 values come from integer-based part of the record and outliers'),
    c(4355500, 'removed', "regulated, outlier 0 flow values"),
    c(4357510, 'removed', "single flow intermittency event, probably gap in data"),
    c(4769200, 'removed', 'only 1 occurrence of 0 flow values'),
    c(4773050, 'removed', 'abrupt decreases to 0'),
    c(5101020, 'removed', "single flow intermittency event, probably gap in data"),
    c(5101101, 'removed', "single flow intermittency event, probably gap in data"),
    c(5101130, 'removed', 'abrupt decreases to 0'),
    c(5101201, 'removed', 'unreliable record, large data gap as 0 flow values'),
    c(5101290, 'removed', '0 flow values before 2000 are outliers, changed flow permanence'),
    c(5101305, 'removed', 'most 0 values look like outliers'),
    c(5101380, 'removed', 'abrupt decreases to 0'),
    c(5109200, 'removed', 'unreliable record, interpolation, large data gap as 0 flow values'),
    c(5109230, 'removed', 'abrupt decreases to 0'),
    c(5202140, 'removed', 'abrupt decreases to 0'),
    c(5202145, 'removed', 'most 0 values look like outliers'),
    c(5202185, 'only one maybe erroneous 0 flow values in first 34 years'),
    c(5204140, 'only one 0 flow values in first 30 years'),
    c(5202228, 'removed', 'maybe regulated, unreliable record post 1983 accounts for 0 flow values'),
    c(5204170, 'removed', 'changed flow permanence'),
    c(5302251, 'removed', 'large data gap as 0 flow values, otherwise only one flow intermittency event'),
    c(5302261, 'removed', 'large data gap as 0 flow values'),
    c(5405095, 'removed', 'changed flow permanence from perennial to non-perennial'),
    c(5405105, 'only one 0 flow value'),
    c(5608100, 'removed', 'large data gap as 0 flow values'),
    c(5708200, 'removed', 'changed flow permanence'),
    c(5803160, 'removed', 'large data gap as 0 flow values'),
    c(5864500, 'removed', 'rounded to 10L/s'),
    c(5870100, 'removed', 'rounded to 100L/s'),
    c(6119100, 'removed', 'rounded to 10L/s'),
    c(6125680, 'removed', 'changed flow permanence'),
    c(6139140, 'removed', 'only one 0 flow occurrence in first 30 years'),
    c(6233410, 'removed', '0 flow values from tidal reversals'),
    c(6442300, 'removed', '0 values come from integer-based part of the record and outliers'),## perfect example of what an integer-based record involves
    c(6444250, 'removed', '0 values come from integer-based part of the record and outliers'),
    c(6444350, 'removed', '0 values come from integer-based part of the record and outliers'),
    c(6444400, 'removed', 'abrupt decreases to 0'),
    c(6935570, 'removed', 'rounded to 10L/s')
  ) %>%
    do.call(rbind, .) %>%
    as.data.table %>%
    setnames(c('GRDC_NO', 'flag', 'comment'))
  
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
  GRDCtoremove_pereartifacts <- list(
    c(1159800, 'removed', 'regulated'),
    c(1160324, 'removed', 'regulated'),
    c(1160331, 'removed', 'low-flow plateaus are likely overestimated 0 values'),
    c(1160520, 'removed', 'regulated by diversions'),
    c(1160602, 'removed', 'regulated'),
    c(1160709, 'removed', 'regulated'),
    c(1160788, 'removed', 'low-flow plateaus may be overestimated 0 values'),
    c(1197310, 'removed', 'regulated'),
    c(1255100, 'inspected', 'now regulated, but naturally perennial, Cunene River'),
    c(1593100, 'inspected', 'bad quality but clearly not IRES'),
    c(1593751, 'inspected', 'missing values may contain intermittency, strange regular patterns, maybe interpolated'),
    c(2357750, 'inspected', 'not regulated, in Sri Lanka'),
    c(2588707, 'removed', 'regulated'),
    c(3628200, 'removed', 'appears to change flow permanence'),
    c(3650634, 'removed', 'regulated, probably changed flow permanence'),
    c(3652030, 'removed', '0s in missing years and low flows in other years may also be 0s'),
    c(4101200, 'removed', 'low-flow plateaus are likely overestimated 0 values'),
    c(4115225, 'removed', 'regulated, may have been IRES otherwise'),
    c(4118850, 'removed', 'low-flow plateaus may be overestimated 0 values'),
    c(4125903, 'removed', 'regulated, may have been IRES otherwise'),
    c(4126351, 'removed', 'regulated, may have been IRES otherwise'),
    c(4148850, 'removed', 'many 0 flow values in missing years'),
    c(4151801, 'removed', 'regulated, may have been IRES otherwise, Rio Grande'),
    c(4152651, 'removed', 'regulated by blue mesa reservoir, may have been IRES otherwise, Rio Grande'),
    c(4208610, 'removed', 'too much missing data but if not would be IRES'),
    c(4213055, 'removed', 'too much missing data but if not would be IRES'),
    c(4213802, 'removed', 'identical to 4213801'),
    c(4214320, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
    c(4362100, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
    c(5606090, 'removed', 'low-flow plateaus may be overestimated 0 values'),
    c(5606414, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
    c(6123630, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
    c(6335020, 'removed', 'identical to 6335060'),
    c(6335050, 'removed', 'identical to 6335060'),
    c(6337503, 'removed', 'regulated, cannot tell whether may have been intermittent before'),
    c(6442100, 'removed', 'identical to 6442600'),
    c(6935146, 'removed', 'identical to 6935145'),
    c(6335050, 'removed', 'identical to 6335060'),
    c(6935600, 'removed', 'identical to 6935145')
  ) %>%
    do.call(rbind, .) %>%
    as.data.table %>%
    setnames(c('GRDC_NO', 'flag', 'comment'))
  
  #---------- Check flags in winter IR
  plot_winterir(dt = GRDCstatsdt, dbname = 'grdc', inp_resdir = inp_resdir,
                yearthresh = 1800, plotseries = plotseries)
  #Checked for seemingly anomalous 0s. Sudden decreases.
  #Check for flags, check satellite imagery, station name, check for construction of reservoir
  
  #------ Check time series of stations within 3 km of seawater
  GRDCcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GRDCstatsdt,
                                     dbname = 'grdc', inp_resdir = inp_resdir,
                                     yearthresh = 1800, plotseries = plotseries)
  #GRDCcoastalirall[, unique(readformatGRDC(path)$Flag), by=GRDC_NO]
  #Nothing obviously suspect beyond those that ad already been flagged
  
  #Inspect statistics for 4208857, 4213531 as no flow days occurred only one year
  # ID = '6976300'
  # GRDCstatsdt[GRDC_NO == ID,]
  # check <- readformatGRDC(GRDCstatsdt[GRDC_NO == ID,path])
  # unique(check$Flag)
  #
  # plotGRDCtimeseries(GRDCstatsdt[GRDC_NO == ID,], outpath=NULL)
  
  
  ### Summarize removal ########################################################################
  
  #Before cleaning
  GRDCflags <- rbindlist(list(GRDCtoremove_allinteger,
                              GRDCtoremove_unstableIR,
                              GRDCtoremove_irartifacts,
                              GRDCtoremove_pereartifacts
  ))
  
  GRDCtoremove_all <- GRDCflags[flag=='removed', GRDC_NO]
  
  GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10, .N]
  GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10 &
                !(GRDC_NO %in% GRDCtoremove_all), .N]
  
  ### Check changes in GSIM discharge data availability and flow regime over time ####
  GSIMflags <- rbindlist(list(GSIMtoremove_irartifacts,
                              GSIMtoremove_winterIR,
                              GSIMtoremove_unstableIR,
                              GSIMtoremove_USroundingfuckup,
                              GSIMtoremove_perartifacts
  ))
  
  #Add US_0004962 as it used to be intermittent https://pubs.usgs.gov/pp/1277/report.pdf
  
  GSIMstatsdt_clean <- GSIMstatsdt[!(gsim_no %in%  GSIMflags[flag=='removed', gsim_no]),]
  
  
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
  
  # plots <- grid.arrange(
  #     ggalluvium_gaugecount(dtformat = alluv_formatGRDC, alluvvar = 'GRDC_NO'),
  #     ggalluvium_gaugecount(dtformat = alluv_formatGSIM, alluvvar = 'gsim_no'),
  #     ggirsensi
  # )
  
  ### Bind GRDC and GSIM records ####################################
  databound <- rbind(GRDCstatsdt_clean,
                     GSIMstatsdt_clean,
                     use.names=TRUE, fill=T)
  
  return(list(#plots=plots,
    data=databound,
    flags=rbind(GRDCflags, GSIMflags,
                use.names=TRUE, fill=T)))
}

#------ format_gaugestats --------------------------------------------------------
#' Format gauge statistics
#'
#' Final selection of gauging stations and formatting of gauge statistics + 
#' hydro-environmental attributes.
#'
#' @param in_gaugestats data.table of summary time series statistics and intermittency statistics for 
#' all gauging stations (GRDC +  GSIM) not excluded in QA/QC process.
#' @param in_gaugep table of hydro-environmental attributes for all gauging stations.
#' @param yearthresh (integer) minimum year from which to analyze/consider discharge record.
#' 
#' @details An additional set of gauges are excluded in this step:
#' - All gauges with less than 10 years of data are excluded (considering only years with no more than 20 days of missing data)
#' - All gauges with a Degree of Regulation >= 50% are excluded (Lehner et al. 2011)
#' 
#' @return data.table of gauging stations included in all subsequent analysis, with
#' summary time series statistics and hydro-environmental attributes.
#'
#' @source Lehner, B., Liermann, C.R., Revenga, C., Vörösmarty, C., Fekete, B., 
#' Crouzet, P., Döll, P., Endejan, M., Frenken, K., Magome, J., Nilsson, C., 
#' Robertson, J.C., Rödel, R., Sindorf, N. and Wisser, D. (2011), 
#' High-resolution mapping of the world's reservoirs and dams for sustainable 
#' river-flow management. Frontiers in Ecology and the Environment, 9: 494-502. 
#' \link{https://doi.org/10.1890/100125}
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
    .[, DApercdiff := (area_correct-UPLAND_SKM)/UPLAND_SKM]
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
  
  print(paste0('Removing ', gaugestats_joinsel[dor_pc_pva >= 500, .N],
               ' stations with >=50% flow regulation'))
  gaugestats_derivedvar <- gaugestats_joinsel[dor_pc_pva < 500, ] %>% #Only keep stations that have less than 50% of their discharge regulated by reservoir
    comp_derivedvar #Compute derived variables, rescale some variables, remove -9999
  
  return(gaugestats_derivedvar)
}

#------ selectformat_predvars -----------------
#' Select + format predictor variables
#'
#' Select candidate predictor variables to use in subsequent modelling and create
#' a table of formatted names.
#'
#' @param inp_riveratlas_meta path to metadata table for river atlas variables
#' @param in_gaugestats data.table of formatted gauging station summary statistics and hydro-environmental attributes.
#' 
#' 
#' @return data.table of selected predictor variable codes, names, category, 
#' attribute, sources, references, etc. Most sources are from RiverATLAS 
#' technical documentation available at https://www.hydrosheds.org/page/hydroatlas.
#' 
#' @export
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
#' Create tasks
#'
#' Create machine learning \link[mlr3]{Task} list for random forest models.
#'
#' @param in_gaugestats data.table of formatted gauging station summary statistics and hydro-environmental attributes. Output from \link{format_gaugestats}.
#' @param in_predvars data.table of candidate predictor variable codes and names. Output from \link{selectformat_predvars}.
#' @param id_suffix (character) suffix to add after task names
#' @param include_discharge (logical) whether to include predictor variables of long-term discharge modeled from WaterGAP v2.2. 
#' 
#' @details Tasks are objects that contain the (usually tabular) data and 
#' additional meta-data to define a machine learning problem. The meta-data is,
#'  for example, the name of the target variable for supervised machine learning 
#'  problems, or the type of the dataset (e.g. a spatial or survival). 
#'  This information is used for specific operations that can be performed on a task. 
#' For more info on tasks, see \link{https://mlr3book.mlr-org.com/tasks.html}.
#' 
#' @return list of three \link[mlr3]{Task}:
#' \itemize{
#'   \item(a spatiotemporal classification task) |link[mlr3spatiotempcv]{TaskClassifST}; 
#'   for classification random forest with spatiotemporal CV (standard probability and probability CIF)
#'   \item(a standard regression tasks) |link[mlr3]{TaskRegr}; 
#'   for regression random forest (MAXSTAT)
#'   \item(a regression task where the minority class has already been oversampled) 
#'   |link[mlr3]{TaskRegr}; for regression random forest (MAXSTAT). 
#'   Implemented because class oversampling pipelines are not available for regression tasks.
#' } 
#' 
#' @export
create_tasks <- function(in_gaugestats, in_predvars,
                         id_suffix=NULL, include_discharge = TRUE) {
  #Create subset of gauge data for analysis (in this case, remove records with missing soil data)
  datsel <- in_gaugestats[, c('intermittent_o1800',in_predvars$varcode, 'X', 'Y'),
                          with=F] %>%
    na.omit
  
  #Remove WaterGAP variables if include_discharge == FALSE
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
#' Create learners
#'
#' Create machine learning \link[mlr3]{Learner} list for random forest model.
#'
#' @param in_task \link[mlr3]{Task} of type |link[mlr3spatiotempcv]{TaskClassif}, 
#' |link[mlr3spatiotempcv]{TaskClassifST}, or |link[mlr3]{TaskRegr}.
#' 
#' @details Learners are methods to train and predict a model for a Task and provide 
#' meta-information about the learners, such as the hyperparameters you can set.
#' For more info on mlr3 Learners, see \link{https://mlr3book.mlr-org.com/learners.html}.
#' 
#' 
#' @return This function will return different types of learners depending on the type of task input. \cr
#' If a |link[mlr3spatiotempcv]{TaskClassif} or |link[mlr3spatiotempcv]{TaskClassifST} is provided, 
#' then a list is returned with six learners, a classif.ranger (standard probability random forest learner) and 
#' classif.cforest (probability conditional inference forest learner), each in three formats for class
#' imbalance correction: 1. without any class imbalance correction (e.g., lrn_ranger), 
#' 2. with random oversampling of minority class (e.g., 'lrn_ranger_overp), and 3. with class weighting (e.g., 'lrn_ranger_weight'). \cr
#' 
#' If a |link[mlr3]{TaskRegr} is provided, then a single learner is returned of type 'regr.ranger'
#' with maximally selected rank statistics (maxstat) as splitting rule. \cr
#' 
#' For all learners, 800 trees are grown.
#' 
#' @export
create_baselearners <- function(in_task) {
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
#' Set tuning strategy
#'
#' Create a \link[mlr3tuning]{AutoTuner} to set the random forest hyperparameter
#' tuning strategy.
#'
#' @param in_learner \link[mlr3]{Learner} whose hyperparameters to tune.
#' @param in_measure performance measure to use for hyperparameter tuning. 
#' The hyperparameter configuration that optimizes this measure will be selected and used for final model training.
#' @param nfeatures total number of predictor variables in the model.
#' @param insamp_nfolds number of cross-validation folders used for tuning
#' @param insamp_neva number of times the cross-validation must be conducted (e.g. if 2, then twice-repeated CV)
#' @param insamp_nbatch number of hyperparameter configurations to evaluate at the same time. This will dictate how many
#' processing cores will be used in hyperparameter tuning. The greater this number, the faster the tuning, but also the 
#' more computing intensive. This shouldn't be set higher than the number of cores on the computer used.
#' 
#' @details For more information on hyperparameter tuning and a table of the range of 
#' hyperparameters and tuning strategies used, see sections IVa and IVb in the Supplementary
#' Information of Messager et al. 2021 at \link{https://www.nature.com/articles/s41586-021-03565-5}.
#' 
#' 
#' @return a \link[mlr3tuning]{AutoTuner}.
#' 
#' @export
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
  
  #learnertune$store_tuning_instance = FALSE
  learnertune$id <- in_learner$id
  
  return(learnertune)
}

#------ set_cvresampling ---------------
#' Set cross-validation resampling
#'
#' Create a \link[mlr3]{Resampling} for model cross-validation.
#'
#' @param rsmp_id type of resampling, 'cv', 'repeated_cv', or 'repeated-spcv-coords'
#' @param in_task \link[mlr3]{Task} to resample
#' @param outsamp_nrep number of cross-validation repetitions
#' @param outsamp_nfolds number of cross-validation folds.
#' 
#' @details For more information on the cross-validation strategy used in this study, see 
#' section IV in the Supplementary Information of Messager et al. 2021 at \link{https://www.nature.com/articles/s41586-021-03565-5}.
#' For general information on resampling in mlr3, see \link{https://mlr3book.mlr-org.com/resampling.html}.
#' 
#' @return a \link[mlr3]{mlr_resamplings_repeated_cv} or \link[mlr3spatiotempcv]{ResamplingRepeatedSpCVCoords}.
#' 
#' @export
set_cvresampling <- function(rsmp_id, in_task, outsamp_nrep, outsamp_nfolds) {
  #repeated_cv or repeated-spcv-coords
  outer_resampling = rsmp(rsmp_id,
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)
  outer_resampling$instantiate(in_task)
  
  return(outer_resampling)
}


#------ dynamic_resample ------------------
#' Dynamic resampling
#'
#' Runs a train-predict-test routine a.k.a. a resampling (possibly in parallel): 
#' Repeatedly apply learner on a training set of task to train a model, 
#' then use the trained model to predict observations 
#' of a test set. Training and test sets are defined by the resampling.
#'
#' @param in_task \link[mlr3]{Task} to resample
#' @param in_learner \link[mlr3]{Learner} to use in resampling (e.g., conditionel inference forest with minority class oversampling). 
#' Can be a simple learner or an \link[mlr3tuning]{Autotuner} already including hyperparameter tuning strategy.
#' @param in_resampling \link[mlr3]{Resampling} cross-validation strategy
#' @param type (character) Type of learner 'classif' or 'regr'. Other values are not accepted.
#' @param store_models  whether to keep the fitted model after the test set has been predicted. Set to TRUE if you want to further analyse the models or want to extract information like variable importance.
#' 
#' @details The dynamic aspect of this model is that it runs 'reset_tuning' on the fly
#' to make sure that the hyperparameter search space matches the task (e.g., if the number of candidate 
#' predictor variables has been reduced, it adjusts mtry)
#' 
#' @return \link[mlr3]{ResampleResult}
#' 
#' @export
dynamic_resample <- function(in_task, in_learner, in_resampling, type,
                             store_models = FALSE) {
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
#' Dynamic resampling using learners from benchmarking results
#'
#' Run a train-predict-test routine using a learner from a \link[mlr3]{BenchmarkResult}: 
#' Repeatedly apply learner on a training set of task to train a model, 
#' then use the trained model to predict observations 
#' of a test set. Training and test sets are defined by the resampling.
#'
#' @param in_task \link[mlr3]{Task} to resample
#' @param in_bm \link[mlr3]{BenchmarkResult} from which to extract learner or name of 
#' serialized file (.qs) on disk (accessed through \link[qs]{qread})
#' @param in__lrnid id of learner to extract from \link[mlr3]{BenchmarkResult} (e.g., "oversample.classif.ranger")
#' @param in_resampling \link[mlr3]{Resampling} cross-validation strategy
#' @param type (character) Type of learner 'classif' or 'regr'. Other values are not accepted.
#' @param inp_resdir (character) path to where qs file is located (excluding the name of the qs file)
#' @param store_models  whether to keep the fitted model after the test set has been predicted. Set to TRUE if you want to further analyse the models or want to extract information like variable importance.
#' 
#' @details The dynamic aspect of this model is that it runs 'reset_tuning' on the fly
#' to make sure that the hyperparameter search space matches the task (e.g., if the number of candidate 
#' predictor variables has been reduced, it adjusts mtry)
#' 
#' @return \link[mlr3]{ResampleResult}
#' 
#' @export
#Run resample on in_task, selected learner (in_lrnid) from in_bm, in_resampling
dynamic_resamplebm <- function(in_task, in_bm, in_lrnid, in_resampling, type,
                               inp_resdir = NULL, store_models = FALSE) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
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
#' Combine benchmarking results
#'
#' Combine a list of \link[mlr3]{ResampleResult} into a single \link[mlr3]{BenchmarkResult}
#'  object and save the object in serialized form on disk.
#
#' @param in_resampleresults one or a list of \link[mlr3]{ResampleResult} objects to convert to \link[mlr3]{BenchmarkResult}
#' @param write_qs whether to write out serialized (qs) file
#' @param inp_resdir path to directory where to write out qs file
#' 
#' @details qs file is writen to \code{inp_resdir/combine_bmYYMMDDHHSS.qs}
#' 
#' 
#' @return file name of qs file
#' 
#' @export
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
  } 
  else {
    warning('You provided only one resample result to combine_bm,
            simply returning output from as_benchmark_result...')
    bmrbase = as_benchmark_result(in_resampleresults[[1]])
  }
  print('Done combining, now writing to qs...')
  if (write_qs) {
    out_filen <- paste0('combine_bm', format(Sys.time(), '%Y%m%d%H%M%s'), '.qs')
    out_qs <- file.path(inp_resdir, out_filen)
  }
  
  qs::qsave(bmrbase, out_qs)
  
  return(out_filen)
}

#------ select_features ------------------
#' Select features
#'
#' Select a subset of model predictor variables (aka features) based on variable importance p-value.
#
#' @param in_bm \link[mlr3]{BenchmarkResult} containing the \link[mlr3]{ResampleResult} 
#' for the learner of interest, based on which to select predictor variables. Can also be name of 
#' serialized file (.qs) on disk (accessed through \link[qs]{qread}) containing BenchmarkResult object.
#' @param in_lrnid id of learner to extract from \link[mlr3]{BenchmarkResult} (e.g., "oversample.classif.ranger")
#' @param in_task \link[mlr3]{Task} containing predictor variables to subset
#' @param pcutoff variable p-value cut-off under which to select variable, see \link{extract_impperf_nestedrf} and
#' \link[ranger]{importance_pvalues} for more information on p-value calculation.
#' @param inp_resdir (character) path to where qs file is located (excluding the name of the qs file)
#' 
#' @return list containing the original \code{in_task} and the updated task only included selected predictor variables.
#' 
#' @export
select_features <- function(in_bm, in_lrnid, in_task, pcutoff, inp_resdir = NULL) {
  
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  #get desired resampled_results/learner
  if (inherits(in_bm, "BenchmarkResult")) {
    in_rf <- in_bm$filter(learner_ids = in_lrnid)
  } else {
    in_rf <- as_benchmark_result(in_bm)
  }
  
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
#' Select and train final random forest
#'
#' Select final random forest learner and train it, optionally adjusting inner sampling strategy
#' for hyperparameter tuning (i.e., cross-validation strategy).
#
#' @param in_rf \link[mlr3]{ResampleResult} of learner to use or \link[mlr3]{BenchmarkResult} from which to extract learner.
#' @param in_lrnid id of learner to extract from \link[mlr3]{BenchmarkResult} (e.g., "oversample.classif.ranger")
#' @param in_task \link[mlr3]{Task} containing predictor variables to subset.
#' @param insamp_nfolds (optional) number of cross-validation folds to adjust in inner (hyperparameter tuning) cross-validation
#' @param insamp_nevals (optional) number of cross-validation repetitions in inner (hyperparameter tuning) cross-validation
#' 
#' @return list containing the outer resampling (i.e. performance cross-validation) results named ('rf_outer'), 
#' the trained learner ('rf_inner') and the task on which it was trained ('task').
#' 
#' @export
selecttrain_rf <- function(in_rf, in_learnerid=NULL, in_task = NULL,
                           insamp_nfolds =  NULL, insamp_nevals = NULL) {
  
  outlist <- list()
  
  ######### Prepare autotuner for full training ####################
  # If a ResampleResult was provided
  if (inherits(in_rf, 'ResampleResult')) {
    in_bmsel <- in_rf$clone()
    lrn_autotuner <- in_bmsel$learner
    in_task <- in_bmsel$task
    outlist[['rf_outer']] <- in_rf
    
    # If a BenchmarkResult was provided
  } else if (inherits(in_rf, 'BenchmarkResult')) {
    in_bmsel <- in_rf$clone()$filter(learner_ids = in_learnerid,
                                     task_id = in_task)
    
    lrn_autotuner <- in_bmsel$clone()$learners$learner[[1]]
    in_task <-in_bmsel$tasks$task[[1]]
    
    #Return outer sampling object for selected model (or list of outer sampling objects)
    uhashes <- unique(as.data.table(in_bmsel)$uhash)
    if (length(uhashes) == 1) {
      outlist[['rf_outer']] <- in_bmsel$resample_result(uhash=uhashes)
    } else {
      outlist[['rf_outer']] <- lapply(uhashes, function(x) {
        in_bmsel$resample_result(uhash=x)
      })
    }
  } else if (inherits(in_rf, 'AutoTuner')) {
    lrn_autotuner <- in_rf
  }
  
  if (!is.null(insamp_nfolds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- insamp_nfolds
  }
  
  if (!is.null(insamp_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- insamp_nevals
  }
  
  ######### Train it ####################
  lrn_autotuner$param_set$values = mlr3misc::insert_named(
    lrn_autotuner$param_set$values,
    list(classif.ranger.importance = 'permutation')
  )
  lrn_autotuner$train(in_task)
  
  outlist[['task']] <- in_task
  outlist[['rf_inner']] <- lrn_autotuner
  
  return(outlist)
}

#------ change_taskbackend -----------------------------------
#' Change task back end
#'
#' Change target variable for task. This is implemented in this model to adjust
#' the criterion used to classify rivers and streams as non-perennial. Changing 
#' the mean annual number of zero-flow day threshold (in the case of this study from 1 to 30 days).
#
#' @param in_task \link[mlr3]{Task} whose target variable to update.
#' @param in_gaugestats data.table of task back-end content (i.e., data to update task with), 
#' ontaining the column 'mDur_o1800' and all hydro-environmental predictor variables.
#' @param newmdurthresh (numeric) criterion in terms of mean annual number of zero-flow days over which to classify a river or stream as non-perennial (e.g., 30)
#' 
#' @return \link[mlr3]{Task} with updated target variable.
#' 
#' @export

change_tasktarget <- function(in_task, in_gaugestats, newmdurthresh) {
  in_gaugestats_new <- copy(in_gaugestats)
  in_gaugestats_new[, intermittent_o1800 := factor(
    as.numeric(mDur_o1800 >= newmdurthresh),
    levels=c('0', '1'))]
  in_task_new <- in_task$clone()
  
  in_task_new$backend <-  as_data_backend(
    in_gaugestats_new[,
                      c(names(in_task_new$data()),in_task$coordinate_names),
                      with=F])
  return(in_task_new)
}

#------ change_oversamp ----------------------
#' Change oversampling ratio
#'
#' Change the class oversampling ratio of an existing random forest \mlr3[mlr3]{Learner} object.
#
#' @param in_oversamprf \link[mlr3]{Learner} whose oversample.ratio parameter to update.
#' @param in_task \link[mlr3]{Task} based on which to determine minority class oversampling ratio.
#' 
#' @details In this study, this function was created to adjust the oversampling ratio after updating
#' the criterion for classifying rivers and streams as non-perennial (see \link{change_taskbackend}).
#' The ratio of perennial to non-perennial gauges increases a lot when only considering watercourses
#' that cease to flow at least 30 days per year as non-perennial. This function enables one to keep
#' all settings and hyperparameters intact for a given Learner while updating that ratio.
#' 
#' @return input \link[mlr3]{Learner} (\code{in_oversamprf}) with updated oversampling ratio.
#' 
#' @export

change_oversamp <- function(in_oversamprf, in_task) {
  oversamprf <- in_oversamprf$clone()
  oversamprf$param_set$values$oversample.ratio <- get_oversamp_ratio(in_task)$ratio
  return(oversamprf)
}


##### --------------------- Predictions ---------------------
#------ rformat_network ------------------
#' Format river network
#'
#' Format river network, joining data, selecting variables, computing derived variables, etc.
#
#' @param in_predvars data.table of predictor variable codes, names and attributes. See \link{selectformat_predvars}.
#' @param in_monthlydischarge optional data.table of long-term mean monthly 
#' discharge values (from WaterGAP v2.2) for global river reaches.
#' @param inp_riveratlasmeta path to RiverATLAS attributes metadata table.
#' @param inp_riveratlas path to RiverATLAS (v1.0) attribute table (for the entire global river network).
#' @param inp_riveratlas2 path to new attributes calculated for global river network following RiverATLAS methodology
#' 
#' 
#' @return data.table of identifiers and hydro-environmental attribute for all reaches in the
#' RiverATLAS global river
#' 
#' @export
rformat_network <- function(in_predvars, in_monthlydischarge=NULL,
                            inp_riveratlasmeta, inp_riveratlas, inp_riveratlas2) {
  cols_toread <-  unique(
    c("HYRIV_ID", "HYBAS_L12", "HYBAS_ID03", "LENGTH_KM",'INLAKEPERC',
      'PFAF_ID05',
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
#' Make gauge predictions
#'
#' Make predictions for all gauging stations used in model training based on final trained model.
#' Join outputs of prediction and cross-validations to original gauging stations statistics table.
#
#' @param in_rftuned trained random forest model, output from \link{ selecttrain_rf}
#' @param in_gaugestats data.table of formatted gauging station summary statistics and hydro-environmental attributes. Output from \link{format_gaugestats}.
#' @param in_res_spcv \link[mlr3]{ResampleResult} from spatial cross-validation, output from \link{dynamic_resamplebm}.
#' @param in_predvars data.table of predictor variable codes, names and attributes. See \link{selectformat_predvars}.
#' @param interthresh (numeric) between 0 and 1 (inclusive), probability threshold 
#' above which to classify gauging stations as non-perennial (e.g., 0.50).
#' @param simple (logical) if TRUE, return only final model predictions, if FALSE return cross-validation predictions
#' 
#' @return \code{in_gaugestats} data.table with one new column (IRpredprob_full, if simple == T)
#' or five new columns (if simple == F). \cr
#' \itemize{
#'   \item IRpredprob_full: predicted probability that station is non-perennial based on final model training. 
#'   \item IRpredprob_CVsp: predicted probability that station is non-perennial based on spatial cross-validation
#'    (i.e., prediction when fold whose station is a member was excluded from training)
#'   \item IRpredprob_CVnosp: predicted probability that station is non-perennial based on spatial cross-validation
#'   \item inter_class_u10_featsel_spfold orinter_class_o1_featsel_spfold: cross-validation fold membership (id)
#'\}
#' 
#' @export
make_gaugepreds <- function(in_rftuned, in_gaugestats, in_res_spcv = NULL,
                            in_predvars, interthresh, simple = FALSE) {
  
  #Get fully trained predictions
  gpreds_full <- in_rftuned$rf_inner$predict(in_rftuned$task)
  
  #Format gauge data.table prior to merging
  datsel <- na.omit(in_gaugestats, c('intermittent_o1800',
                                     in_predvars$varcode, 'X', 'Y'))
  
  datsel[, IRpredprob_full := gpreds_full$prob[,2]]
  
  #If not simple (include CV predictions)
  if (!simple) {
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
    rsmp_res_sp <- get_outerrsmp(in_rftuned, spatial_rsp=TRUE)
    gpreds_CV_sp <-  in_res_spcv$prediction()$set_threshold(1-interthresh) %>%
      as.data.table %>%
      .[, list(truth=first(truth), prob.1=mean(prob.1)), by=row_id] %>%
      .[, response := fifelse(prob.1 >= interthresh, 1, 0)] %>%
      setorder(row_id)
    
    #Merge all predictions
    datsel[, `:=`(
      IRpredprob_CVsp = gpreds_CV_sp$prob.1,
      IRpredprob_CVnosp = gpreds_CV_nosp$prob.1,
      spfold_unique =  spcv_clustersformat$spfold_unique
    )] %>%
      setnames(old='spfold_unique', new=paste0(in_rftuned$task$id, '_spfold'))
    
  }
  return(datsel)
}

#------ bind_gaugepreds -----------------
#' Bind gauge predictions
#'
#' Compile predictions from multiple models. This function was used in the case 
#' when multiple models were developed for different (possibly overlapping) subsets
#' of gauging stations. 
#
#' @param in_gpredsdt either a single data.table of predictions for gauges (output from \link{make_gaugepreds}) 
#' or a list of data.tables of model predictions to bind.
#' @param interthresh either a single value or a data.table of values of the
#'  probability threshold to assign predicted flow intermittence classes to 
#'  gauging stations (based on the predicted probability of flow intermittence).
#'  
#' @details when several predictions are provided for the same gauging station (as was the
#' case for this study for gauges with a long-term mean annual flow between 1 and 10 m3/s), the
#' average predicted probability of flow intermittence was computed.
#' 
#' @return data.table compiling all predictions with new columns of the predicted 
#' categorical flow intermittence class ("IRpredcat_') and flow intermittence 
#' prediction residuals (IPR; named preduncert_).
#' The column suffixes reflect whether the predictions are based on the
#'  final model training or on cross-validation results (see \link{make_gaugepreds}).
#' 
#' @export
bind_gaugepreds <- function(in_gpredsdt, interthresh) {
  #Get binary classification threshold
  if (inherits(interthresh, 'data.table')) {
    interthresh <- interthresh[learner == in_learnerid, thresh]
  }
  
  if (is.list(in_gpredsdt) & !inherits(in_gpredsdt, 'data.table')) {
    in_gpredsdt <- rbindlist(in_gpredsdt,
                             use.names = TRUE, fill = TRUE, idcol = "modelgroup")
  }
  
  
  #Compute mean predicted probability if overlapping models
  out_gpredsdt <- in_gpredsdt[, IRpredprob_full := mean(IRpredprob_full, na.rm=T),
                              by=.(GAUGE_NO)]
  
  if ('IRpredprob_CVnosp' %in% names(in_gpredsdt)) {
    out_gpredsdt[, IRpredprob_CVnosp := mean(IRpredprob_CVnosp, na.rm=T),
                 by=.(GAUGE_NO)]
  }
  
  if ('IRpredprob_CVsp' %in% names(in_gpredsdt)) {
    out_gpredsdt[, IRpredprob_CVsp := mean(IRpredprob_CVsp, na.rm=T),
                 by=.(GAUGE_NO)]
  }
  
  out_gpredsdt <- out_gpredsdt[!duplicated(GAUGE_NO),]     #Remove duplicates (if overlapping model)
  
  out_gpredsdt[,`:=`(IRpredcat_full = fifelse(IRpredprob_full >= interthresh,
                                              '1', '0'),
                     preduncert_full = IRpredprob_full -
                       as.numeric(as.character(intermittent_o1800))
  )]
  
  if ('IRpredprob_CVnosp' %in% names(out_gpredsdt)) {
    out_gpredsdt[,`:=`(
      IRpredcat_CVnosp =fifelse(IRpredprob_CVnosp >= interthresh,
                                '1', '0'),
      preduncert_CVnosp = IRpredprob_CVnosp -
        as.numeric(as.character(intermittent_o1800))
    )]
  }
  
  if ('IRpredprob_CVsp' %in% names(out_gpredsdt)) {
    out_gpredsdt[,`:=`(
      IRpredcat_CVsp =fifelse(IRpredprob_CVsp >= interthresh,
                              '1', '0'),
      preduncert_CVsp = IRpredprob_CVsp -
        as.numeric(as.character(intermittent_o1800))
    )]
  }
  
  return(out_gpredsdt)
}



#------ write_gaugepreds -----------------
#' Write gauge predictions
#'
#' Join model predictions for gauging stations to point feature class of gauging stations,
#' then write the result of the merging operation to disk as a GeoPackage.
#
#' @param in_gaugep \link[sf]{sf} points and hydro-environmental attributes for all gauging stations.
#' @param in_gpredsdt final formatted model predictions for all gauging stations. Output from \link{bind_gaugepreds}.
#' @param outp_gaugep full path to which the geopackage will be written, including the file name and extension.
#'  
#' @details for more details on GeoPackage (extension .gpkg), see \link{https://www.geopackage.org/}.
#' 
#' @return path to output geopackage i.e. \code{outp_gaugep}
#' 
#' @export
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
#' Write river network predictions
#'
#' Make model predictions for global river network and write out csv table.
#
#' @param in_network data.table of formatted global river network with all predictor
#'  hydro-environmental variables. Output from \link{rformat_network}.
#' @param in_rftuned trained random forest model. Output from \link{selecttrain_rf}; 
#' list containing inner resampling results.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param predcol (character) name of the column to which the predicted probability 
#' of flow intermittence for each reach be written.
#' @param discharge_interval (vector of two numeric values) discharge criterion 
#' to subset global river network. Range of long-term mean annual discharge values 
#' for which to produce model predictions.
#' @param interthresh (numeric) between 0 and 1 (inclusive), probability threshold 
#' above which to classify river and stream reaches as non-perennial (e.g., 0.50).
#' @param outp_riveratlaspred (character) ull path to which the output csv file will be written,
#'  including the file name and extension.
#'  
#' @details the date is added to the output file name in the format 
#' \code{outp_riveratlaspred_YYYYMMDD}.
#' 
#' @return path to output csv file.
#' 
#' @export
write_netpreds <- function(in_network, in_rftuned, in_predvars, predcol,
                           discharge_interval = c(-Inf, Inf),
                           interthresh = 0.5, outp_riveratlaspred) {
  
  if (is.character(in_network)) {
    in_network <- fread(in_network)
  }
  
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
                 pred := in_rftuned$rf_inner$predict_newdata(
                   as.data.frame(.SD))$prob[,2],]
      toc()
    }
    
    #Label each reach categorically based on threshold
    in_network[, predcat := fifelse(pred>=interthresh, 1, 0)]
    
    return(in_network[, c('HYRIV_ID', 'HYBAS_L12',
                          'pred', 'predcat'), with=F])
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
    networkpreds <- networkpreds[, pred := mean(pred, na.rm=T),
                                 by=.(HYRIV_ID)] %>% #Compute mean predicted probability if overlapping models
      .[!duplicated(HYRIV_ID),] %>%     #Remove duplicates (if overlapping model)
      .[, predcat := fifelse(pred >= interthresh, '1', '0')]
  }
  
  setnames(networkpreds, old=c('pred', 'predcat'),
           new=c(predcol, paste0(predcol, 'cat'))
  )
  
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
#' Format network predictions
#'
#' Read model predictions for network and add river reach and basin identifiers, 
#' geographic information, and hydrological attributes.
#
#' @param outp_riveratlaspred full path to model predictions for river network. 
#' @param in_rivernetwork  data.table of formatted global river network with all 
#' predictor hydro-environmental variables. Output from \link{rformat_network}.
#'  
#' @details the date is added to the output file name in the format \code{outp_riveratlaspred_YYYYMMDD}.
#' 
#' @return data.table with model predictions and other river network attributes, including
#' river reach length (LENGTH_KM), long-term mean annual discharge (dis_m3_pyr), 
#' the percentage of the river reach length that intersects with a lake (INLAKEPERC),
#' the surface area of the river reach's full upstream drainage area (UPLAND_SKM), and
#' the climate zone where the river reach is located (clz_cl_cmj).
#' 
#' @export
netpredformat <- function(outp_riveratlaspred, in_rivernetwork) {
  fread(file_in(outp_riveratlaspred))%>%
    .[in_rivernetwork[, c('HYRIV_ID', 'HYBAS_L12', 'PFAF_ID05',
                          'LENGTH_KM', 'dis_m3_pyr', 'INLAKEPERC',
                          'UPLAND_SKM', 'clz_cl_cmj'),
                      with=F], on='HYRIV_ID']
}
#------ extrapolate_networklength --------
#' Extrapolate river network length
#'
#' Train models and use them to extrapolate the global length of rivers.
#
#' @param inp_riveratlas full path to RiverATLAS (v1.0) attribute table.
#' @param min_cutoff (numeric) minimum discharge to include in model training.
#' @param dispred (numeric vector) discharge values for which to produce model predictions.
#' @param interactive (logical) whether to print results and make plots for user
#'  to interactively evaluate results and troubleshoot.
#' @param grouping_var variable with which to subset the dataset into groups. A separate
#' model is trained on each group.
#'  
#' @details The prevalence of IRES  was independently extrapolated for a total 
#' of 465 spatial sub-units representing all occurring intersections of 
#' 62 river basin regions (BasinATLAS level 2 subdivisions) and 18 climate zones
#'  (Global Environmental Stratification). For each basin–climate sub-unit, we 
#'  first extrapolated the empirical cumulative distribution of total stream length 
#'  (of all reaches with MAF ≥ 0.1 m3 s−1) down to 0.01 m3 s−1 MAF using a generalized 
#'  additive model (GAM). We excluded reaches larger than the 95th percentile of 
#'  MAF (that is, the largest rivers) within the sub-unit from model fitting to avoid common
#'  discontinuities at the high end of the empirical distribution
#'  that can affect the low end of the power-law-like trendline.
#' 
#' @return list containing a data.table and two plots.
#' The data.table contains an estimate stream length for each discharge size class
#' and climate-basin sub-unit. cumL_pred is the cumulative length of all rivers 
#' and streams with MAF > dispred. cumL_cutoffref is the cumulative river length
#' at MAF == min_cutoff. cumL_predextra = cumL_pred - cumL_cutoffref.
#' 
#' 
#' @export
extrapolate_networklength <- function(inp_riveratlas,
                                      min_cutoff = 0.1,
                                      dispred = seq(0.01, 0.09, 0.01),
                                      interactive = T,
                                      grouping_var = 'PFAF_IDclz') {
  
  #Glossary - CCDF: Complementary Cumulated Distribution Function
  
  #Columns to import from full network
  incols <- c('HYRIV_ID', 'INLAKEPERC', 'PFAF_ID05',
              'dis_m3_pyr', 'LENGTH_KM', 'clz_cl_cmj')
  #Import global river network and join to predictions
  riveratlas <- fread_cols(file_name=inp_riveratlas,
                           cols_tokeep = incols) %>%
    setorder(dis_m3_pyr) %>%
    .[, LENGTH_KM_NOLAKE := LENGTH_KM*(1-INLAKEPERC)] %>%
    .[LENGTH_KM_NOLAKE > 0,]
  
  if (grouping_var == 'PFAF_IDclz') {
    riveratlas[, PFAF_IDclz := paste0(floor(PFAF_ID05/1000), '_', clz_cl_cmj)]
  }
  
  #Use level 3 basins because otherwise, too many level 5 basins don't have streams with >= 0.1 m3/s discharge and can't be modeled
  
  #----------- Function to plot empirical CCDF ------------------------------------------
  plot_ccdf <- function(dt_format, minthresh, maxthresh) {
    ggplot(dt_format[(dis_m3_pyr >= minthresh) & (dis_m3_pyr < maxthresh)],
           aes(x = dis_m3_pyr, y = cumP)) +
      geom_step(size=1.5) +
      geom_step(data=dt_format[(dis_m3_pyr < minthresh),],
                linetype='dashed') +
      geom_step(data=dt_format[(dis_m3_pyr >= maxthresh),],
                linetype='dashed') +
      geom_smooth(method='gam', color='red', fullrange = T, alpha=0.7) +
      scale_x_log10(name = bquote('Naturalized long-term mean annual discharge'~(m^3~s^-1)),
                    breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      annotation_logticks() +
      theme_bw()
  }
  
  #----------- Visualize CCDF plot for the world ------------------------------
  if (interactive) {
    ccdf_datall  <- riveratlas[dis_m3_pyr >= 0.001, .N, by=dis_m3_pyr] %>%
      setorder(-dis_m3_pyr) %>%
      .[, cumN := cumsum(N)] %>%
      .[, cumP := cumN/riveratlas[,.N]]
    
    plot_ccdf(ccdf_datall, minthresh = min_cutoff, maxthresh = 10000)
  }
  
  
  #----------- Compute empirical cumulative distribution by basin --------------
  #Compute cumulative length by discharge
  cumL_bas03 <- riveratlas[dis_m3_pyr > 0, list(suml = sum(LENGTH_KM_NOLAKE)),
                           by=c('dis_m3_pyr', grouping_var)] %>%
    setorder(-dis_m3_pyr) %>%
    .[, cumL := cumsum(suml), by = grouping_var]
  
  #Compute cumulative distribution by number and merge to length
  print('Compute empirical cumulative distribution')
  ccdf_datbas03  <- riveratlas[dis_m3_pyr > 0, .N,
                               by=c('dis_m3_pyr', grouping_var)] %>%
    setorder(-dis_m3_pyr) %>%
    .[, cumN := cumsum(N), by = grouping_var] %>%
    merge(riveratlas[,list(totalN = .N), by=grouping_var], by=grouping_var) %>%
    .[, cumP := cumN/totalN] %>%
    merge(cumL_bas03, by = c(grouping_var, 'dis_m3_pyr')) %>%
    merge(.[dis_m3_pyr >= min_cutoff, list(mindis_o01 = min(dis_m3_pyr)),
            by=grouping_var]) %>%
    merge(.[(dis_m3_pyr >= min_cutoff & dis_m3_pyr < quantile(dis_m3_pyr, 0.95)),
            list(nuniquedis_o01 = length(unique(dis_m3_pyr))),
            by=grouping_var], by=grouping_var) %>% #Number of unique discharge values >= 0.1 m3/s
    .[!is.na(get(grouping_var)),]
  
  #Check distribution of cut-off values for discharge when computing ECDF GAM
  if (interactive) {
    qplot(ccdf_datbas03[nuniquedis_o01 >= 20, quantile(dis_m3_pyr, 0.95),
                        by=grouping_var]$V1) + scale_x_log10()
  }
  
  #----------- Fit GAM to cumL ~ dist_m3_pyr -----------------------------------
  #Check how many basins/climate zones drop off from the cut off criterion
  if (interactive) {
    ccdf_datbas03[!duplicated(get(grouping_var)),.N]
    ccdf_datbas03[!duplicated(get(grouping_var)) & nuniquedis_o01 >= 20, .N]
    ccdf_datbas03[nuniquedis_o01 < 20, sum(suml)]/ccdf_datbas03[nuniquedis_o01 >= 20, sum(suml)]
    #~3000 km out of 33 million kilometers
  }
  
  #Simple function to fit a gam with basic params
  gamfit_util <- function(in_dt, in_k=-1) {
    fit <- mgcv::gam(log10(cumL) ~ s(log10(dis_m3_pyr), k=in_k, bs = "cs"), #cubic splines are good when lots of data points
                     method = "REML",
                     data=in_dt[dis_m3_pyr < quantile(dis_m3_pyr, 0.95),])
  }
  
  #Function to iterate over k values and check whether dimension is too low
  gamfit_kiter <- function(in_dt, verbose = F, kstep=1, kmax=20) {
    gam_out <- gamfit_util(in_dt)
    k_check <- k.check(gam_out)
    k_set <- k_check[1]
    if (k_check[4] < 0.05) {
      while ((k_check[4] < 0.05) & (k_set < kmax)) {
        if (verbose) print(k_set)
        gam_out <- gamfit_util(in_dt, in_k=k_set)
        k_check <- k.check(gam_out)
        k_set <- k_set + kstep
      }
    }
    
    return(gam_out)
  }
  
  #-----Across all basins:
  #Only train GAM models for basins that have at least 20 unique discharge values >= 0.1
  gambas <- ccdf_datbas03[
    (dis_m3_pyr >= min_cutoff) & nuniquedis_o01 >= 20,
    list(mod = list(gamfit_kiter(.SD, verbose = F, kstep = 2, kmax=20))),
    by = grouping_var
  ]
  
  #------Illustrate approach on a single basin
  fit_rivlengam <- function(dt_format, pfaf_id) {
    gamsubdat <- dt_format[(get(grouping_var) == pfaf_id),]
    gamsubdatsub <- gamsubdat[(dis_m3_pyr >= min_cutoff) &
                                (dis_m3_pyr < quantile(dis_m3_pyr, 0.95)),]
    
    gamfit <- gamfit_kiter(in_dt=gamsubdatsub, verbose=T, kstep = 2, kmax=9)
    
    print(summary(gamfit))
    preds <- data.frame(
      dis_m3_pyr=seq(0.01,
                     quantile(gamsubdat$dis_m3_pyr, 0.95), 0.01))
    preds[, c('cumL', 'cumL_se')] <- predict(gamfit, preds, se.fit=TRUE)
    preds$cumL_95CIlow <- 10^(preds$cumL - 1.96*preds$cumL_se)
    preds$cumL_95CIhigh <- 10^(preds$cumL + 1.96*preds$cumL_se)
    preds$cumL <- 10^preds$cumL
    
    ggplot(gamsubdat, aes(x=dis_m3_pyr, y=cumL)) +
      geom_step(data=gamsubdat, size=1.5, color='black', alpha=1/3) +
      geom_step(data=gamsubdatsub, size=1.5, color='black') +
      geom_line(data=preds, color='red') +
      geom_ribbon(data=preds[preds$dis_m3_pyr<=0.1,],
                  aes(ymin=cumL_95CIlow , ymax=cumL_95CIhigh ),
                  fill='blue', alpha=1/4) +
      geom_line(data=preds[preds$dis_m3_pyr<=0.1,], color='blue') +
      scale_x_log10(name=bquote('Naturalized long-term \nmean annual discharge'~(m^3~s^-1)),
                    limits=c(0.01,
                             2*max(gamsubdat$dis_m3_pyr))) +
      scale_y_log10(name='Cumulative river length (km)') +
      annotation_logticks() +
      theme_classic() +
      theme(text=element_text(size=12),
            axis.title.x = element_text(vjust=-2))
  }
  plot_62_18 <- fit_rivlengam(dt_format = ccdf_datbas03, pfaf_id="62_18")
  plot_14_17 <- fit_rivlengam(dt_format = ccdf_datbas03, pfaf_id="14_17")
  
  if (interactive) {
    print(plot_62_18)
    print(plot_14_17)
  }
  
  #Check the distribution of adjusted R-squares across all basins
  rsqdt <- gambas[, list(rsq_adj = summary(mod[[1]])$r.sq), by=grouping_var]
  ggplot(rsqdt, aes(x=rsq_adj)) + geom_histogram()
  print(paste0(
    "The mean the GAM models' adjusted-R squares across all global basins is ",
    round(rsqdt[, mean(rsq_adj)], 4)
  ))
  
  #Fitting function
  predict_baslen <- function(in_gambas, discheck) {
    predcheck <- in_gambas[, 10^predict(mod[[1]], data.frame(dis_m3_pyr=discheck)),
                           by=grouping_var] %>%
      setnames(old='V1', new='cumL_pred')
  }
  
  #Check fit
  if (interactive) {
    predcheck <- predict_baslen(in_gambas = gambas,
                                discheck = 0.2) %>%
      merge(ccdf_datbas03[dis_m3_pyr==0.2,], by=grouping_var)
    
    summary(lm(log10(cumL_pred)~log10(cumL), data=predcheck))$r.sq
    
    ggplot(predcheck, aes(x=cumL_pred, y=cumL)) +
      geom_point() +
      scale_x_log10(name = 'Predicted cumulative length') +
      scale_y_log10(name = 'Observed cumulative length (RiverATLAS') +
      geom_abline()
  }
  
  #----------- Compute predictions and return them -----------------------------
  print('Extrapolating total river length for large basins...')
  #Predict for basins with at least 20 unique discharge size classes
  pred_largebas <-lapply(dispred, function(x) {
    predict_baslen(in_gambas = gambas, discheck = x) %>%
      .[, dispred := x]
  }) %>%
    rbindlist
  
  pred_all <- pred_largebas
  
  #For really tiny basins, don't even bother
  ccdf_datbas03_tinybas <- ccdf_datbas03[dis_m3_pyr >= min_cutoff &
                                           nuniquedis_o01 < 20 &
                                           mindis_o01 > 0.5,]
  if (ccdf_datbas03_tinybas[, .N] > 0) {
    ccdf_datbas03_tinybas[, sum(suml)]
    extend_tinybas <- lapply(dispred, function(x) {
      sub <- ccdf_datbas03_tinybas[, list(
        cumL_pred = .SD[dis_m3_pyr==min(dis_m3_pyr), cumL]),
        by=grouping_var]
      sub[, dispred := x]
      return(sub)
    }) %>%
      rbindlist
    
    pred_all <- rbind(pred_all, extend_tinybas)
  }
  
  #For basins with less than 20 unique discharge values, use average global ratio of
  #predicted cumulative length at 0.01 m3/s to that at minimum discharge over 0.1 m3/s
  #available in that basin
  ccdf_datbas03_smallbas <- ccdf_datbas03[
    dis_m3_pyr >= min_cutoff & nuniquedis_o01 < 20 & mindis_o01 <= 0.5,]
  
  if (ccdf_datbas03_smallbas[, .N] > 0) {
    print('Extrapolating total river length for small basins...')
    #Predict for basins with at least 20 unique discharge values
    pred_smallbas <- lapply(dispred, function(x) {
      print(x)
      ccdf_datbas03_smallbas[
        , list(cumL_pred = mean(
          predict_baslen(in_gambas = gambas, discheck = x)$cumL_pred/
            predict_baslen(in_gambas = gambas, discheck = min(dis_m3_pyr))$cumL_pred
        ) * .SD[dis_m3_pyr == min(dis_m3_pyr), cumL]
        ), by=grouping_var] %>%
        .[, dispred := x]
    }) %>%
      rbindlist
    
    pred_all <- rbind(pred_all, pred_smallbas)
  }
  
  #Merge predictions for small and large basins
  pred_all_format <- pred_all %>%
    merge(riveratlas[dis_m3_pyr >= min_cutoff,
                     list(cumL_cutoffref = sum(LENGTH_KM_NOLAKE)), by=grouping_var],
          by=grouping_var) %>%
    .[, `:=`(cumL_predextra = round(cumL_pred - cumL_cutoffref, 2),
             cumL_pred = round(cumL_pred, 2))]
  
  #Check total
  pred_all_format[dispred == 0.01, sum(cumL_pred)]
  
  return(list(preds=pred_all_format,
              plot_62_18 = plot_62_18,
              plot_14_17 = plot_14_17)
  )
  ################### EXTRA STUFF NOT USED#############################
  # #
  # #Look at distribution for a given basin
  # if (plot) {
  #   #Of length
  #   ggplot(ccdf_datbas03[(dis_m3_pyr >= 0.2) & (dis_m3_pyr < 10000) &
  #                          (clz_cl_cmj == '62_18'),], aes(x=dis_m3_pyr, y=cumL)) +
  #     geom_line() +
  #     scale_x_log10() +
  #     scale_y_log10()
  #
  #   #Plot CCDF
  #   plot_ccdf(ccdf_datbas03[clz_cl_cmj == '62_18',], minthresh = 0.1, maxthresh = 10000)
  # }
  #
  #
  # #----------- Fit distributions
  # #Test poweRlaw package on example basin
  # testdat <- riveratlas[(dis_m3_pyr >= min_cutoff) &
  #                         (dis_m3_pyr < 10000) &
  #                         (clz_cl_cmj == '62_18'),]
  #
  # testcpl <- testdat[, conpl$new(dis_m3_pyr)]
  # #Set xmin
  # xminvec = seq(0.1, 2, 0.05)
  # testcpl$setXmin(max(min_cutoff,
  #                     testdat[, min(dis_m3_pyr)],
  #                     estimate_xmin(testcpl, xmins = xminvec)$xmin))
  # #Estimate continuous power law distribution parameters (fit max likelihood)
  # test_est = estimate_pars(testcpl)
  # testcpl$setPars(test_est)
  #
  # #Model with other distributions
  # testlnorm <- testdat[, conlnorm$new(dis_m3_pyr)]
  # testexp <- testdat[, conexp$new(dis_m3_pyr)]
  # testweibull <- testdat[, conweibull$new(dis_m3_pyr)]
  #
  # testlnorm$setXmin(testcpl$xmin)
  # testexp$setXmin(testcpl$xmin)
  # testweibull$setXmin(testcpl$xmin)
  #
  # test_lnormest = estimate_pars(testlnorm)
  # test_expest = estimate_pars(testexp)
  # test_weibullest = estimate_pars(testweibull)
  #
  # testlnorm$setPars(test_lnormest)
  # testexp$setPars(test_expest)
  # testweibull$setPars(test_weibullest)
  #
  # #Select distribution
  # sel_mod_util <- function(in_compcomb) {
  #   xwin <- as.numeric(in_compcomb[1, 'pval']<0.05)
  #
  #   win <- in_compcomb[1, 2-xwin][[1]]
  #   lose <- in_compcomb[1, xwin+1][[1]]
  #
  #   if (nrow(in_compcomb) > 1) {
  #     sel_mod_util(in_compcomb = in_compcomb[!((in_compcomb$x %in% lose) |
  #                                                (in_compcomb$y %in% lose)),])
  #   } else{
  #     return(win)
  #   }
  # }
  #
  # sel_mod <- function(modlist) {
  #   modcomb <- combn(modlist, 2)
  #   compcomb <- mapply(function(x, y) {
  #     pval <- compare_distributions(x, y)$p_one_sided
  #     data.frame(x=class(x)[1], y=class(y)[1], pval)
  #   }, modcomb[1,], modcomb[2,]) %>%
  #     t %>%
  #     as.data.frame
  #
  #   selmod_char <- sel_mod_util(in_compcomb=compcomb)
  #   selmod_num <- which(unlist(lapply(modlist, function(x) class(x)==selmod_char)))
  #
  #   return(modlist[[selmod_num]])
  # }
  #
  # mod <- sel_mod(modlist = c(testcpl, testlnorm, testexp, testweibull))
  #
  # check <- dist_rand(mod, n=100000)
  #
  #
  # qplot(check) + scale_x_log10() + scale_y_log10()
  # ggplot(check2) + scale_x_log10() + scale_y_log10() +
  #   geom_histogram(data=jitter(testcpl$dat, 0.05))
  # qplot() + scale_x_log10() + scale_y_log10()
  #
  # ks.test(check, "plnorm",
  #         mod$pars[1], mod$pars[2], alternative="two.sided")
  #
  # ## Plot the data (from xmin)
  # plot(testcpl)
  # lines(testcpl, col = 'red')
  # lines(testlnorm, col = 'blue')
  # lines(testweibull, col = 'green')
  #
  # plot.new()
  # ggplot(plotdat, aes(x, y)) +
  #   geom_line() +
  #   scale_x_log10() +
  #   scale_y_log10() +
  #   stat_smooth(method='gam', formula=y ~ s(x, k=6), color='orange') +
  #   geom_line(data = lines(testcpl, col = 'blue'), color='green') +
  #   geom_line(data = lines(testlnorm, col = 'blue'), color='red') +
  #   geom_line(data = lines(testweibull, col = 'blue'), color='blue')
  #
  # #Extremely slow, unfeasible for some reason
  # tic()
  # check <- bootstrap_p(mod, xmins=c(mod$xmin-0.01, mod$xmin),
  #                      no_of_sims = 10, threads = 5)
  # toc()
  #
  # #Check uncertainty in alpha with different xmin
  # parvar_xmin <- lapply(xminvec, function(xmin) {
  #   testcpl$setXmin(xmin)
  #   estimate_pars(testcpl)$pars
  # }) %>%
  #   unlist() %>%
  #   data.table(xmin=xminvec, alpha=.)
}

#------ extrapolate_IRES -----------
#' Extrapolate prevalence of flow intermittence
#'
#' Train models and use them to extrapolate the global prevalence of intermittent
#' rivers and ephemeral streams (as a percentage of network length).
#'
#' @param in_rivpred output from \link{netpredformat}.
#' @param in_extranet output from \link{extrapolate_networklength}. 
#' @param min_cutoff (numeric) minimum discharge to include in model training.
#' @param interactive (logical) whether to print results and make plots for user
#'  to interactively evaluate results and troubleshoot.
#' @param valuevar (character) name of the column containing the beta-distributed 
#' response variable to extrapolate (based on mean annual flow as the indepedent variable).
#' @param grouping_var variable with which to subset the dataset into groups. A separate
#' model is trained on each group.
#'  
#' @details The prevalence of IRES  was independently extrapolated for a total 
#' of 465 spatial sub-units representing all occurring intersections of 
#' 62 river basin regions (BasinATLAS level 2 subdivisions) and 18 climate zones
#'  (Global Environmental Stratification). For each basin–climate sub-unit, after 
#'  extrapolating the empirical cumulative distribution of total stream length 
#'  (of all reaches with MAF ≥ 0.1 m3 s−1) down to 0.01 m3/s MAF, we extrapolated
#'  the prevalence of flow intermittence (in percentage of stream length) down to 
#'  0.01 m3/s MAF. wefitted a GAM for beta-distributed data—that is, 
#'  with a (0, 1) range—to the prevalence of intermittence in each logarithmic 
#'  MAF size bin of the sub-unit. 
#' 
#' @return list containing a data.table and two plots.
#' The data.table contains an estimate 
#' 
#' 
#' @export
extrapolate_IRES <- function(in_rivpred,
                             in_extranet,
                             min_cutoff = 0.1,
                             interactive = F,
                             valuevar = 'predbasic800cat',
                             grouping_var = 'PFAF_IDclz') {
  
  in_extranet <- in_extranet$preds
  #Format rivp
  binlog <- function(x) {
    #log10 scale
    10^(floor(log10(x)*10)/10)
  }
  
  netsub <- in_rivpred %>%
    .[dis_m3_pyr >= min_cutoff & INLAKEPERC < 1,] %>% #Exclude reaches below discharge cutoff and within lakes
    .[, `:=`(PFAF_IDclz = paste0(floor(PFAF_ID05/1000), '_', clz_cl_cmj), #Get PFAF_ID for HydroBASINS level 2
             dislogbin = binlog(dis_m3_pyr), #Get discharge log bins e.g. [0.1,0.2) ; [0.2,0.3) ... [1,2) ; [2,3) ... [10) ...
             LENGTH_KM_NOLAKE = LENGTH_KM*(1-INLAKEPERC)
    )]
  
  #-------- Compute prevalence of IRES (by length) by discharge bin
  netsub_binir <- netsub[
    , list(percinter_bin = sum(get(valuevar)*LENGTH_KM_NOLAKE)/sum(LENGTH_KM_NOLAKE),
           Nreach = .N),
    by = .(dislogbin)] %>%
    setorder(dislogbin)
  
  #-------- Compute prevalence of IRES by discharge bin AND BASIN
  netsub_basbinir <- netsub[
    , list(percinter_bin = sum(get(valuevar)*LENGTH_KM_NOLAKE)/sum(LENGTH_KM_NOLAKE),
           Nreach = .N),
    by = c('dislogbin', grouping_var)]
  
  #----------- Fit GAM to IRES prevalence ~ binned discharge  -----------------------------------
  fit_basIRESgam <- function(dt, verbose=T, onlyextra=F) {
    #Truncate data to fall between 0.001 and 0.999
    dt_gamformat <- dt[, percinter_bin := fcase(percinter_bin == 0, 0.001,
                                                percinter_bin == 1, 0.999,
                                                default = percinter_bin),
                       by = dislogbin]
    
    preds_input <- data.frame(dislogbin=c(seq(0.01, 0.09, 0.01),
                                          unique(dt$dislogbin)))
    preds <-  try(
      {
        if (dt_gamformat[, length(unique(percinter_bin)) > 1  &
                         diff(range(percinter_bin)) > 0.01]
        ) {
          gamfit <- mgcv::gam(percinter_bin ~ s(log10(dislogbin), bs = "cs"),
                              method = "REML", family= betar(),
                              data=dt_gamformat)
          
          if (verbose){
            print(summary(gamfit))
          }
          
          preds_input$percinter_bin <- as.vector(predict(gamfit, preds_input, type='response'))
          preds_input$adjr2 <- summary(gamfit)$r.sq
        } else {
          preds_input$percinter_bin <- dt_gamformat[, max(percinter_bin)]
        }
        
        preds_input
      }
    )
    if (inherits(preds, 'try-error')) {
      preds <- preds_input
      preds$percinter_bin <- dt_gamformat[, max(percinter_bin)]
      
    }
    
    setDT(preds)
    
    preds[dislogbin < min(dt_gamformat$dislogbin),
          percinter_bin_adjust := max(c(percinter_bin,
                                        dt_gamformat[dislogbin == min(dislogbin),
                                                     percinter_bin])),
          by=dislogbin]
    
    if (onlyextra) {
      preds <- preds[!is.na(percinter_bin_adjust),]
    }
    
    return(preds)
  }
  
  #----------- Plot GAM model  -----------------------------------
  plot_basIRESgam <- function(dt_format, pfaf_id=NULL) {
    
    if (!is.null(pfaf_id)) {
      gamsubdat <- dt_format[(get(grouping_var) == pfaf_id),]
    } else {
      gamsubdat <- dt_format
    }
    
    
    tryCatch(
      {
        preds <- fit_basIRESgam(dt=gamsubdat)
        
        ggplot(gamsubdat, aes(x=dislogbin, y=100*percinter_bin, group=1)) +
          geom_step(size=1, stat='identity') +
          geom_line(data=preds, color='red') +
          geom_line(data=preds[dislogbin <= min(gamsubdat$dislogbin),],
                    aes(y=100*percinter_bin_adjust, group=1),color='blue') +
          scale_x_log10(name=bquote(
            'Binned naturalized long-term \nmean annual discharge'~(m^3~s^-1)),
            limits=c(0.01, 2*max(gamsubdat$dis_m3_pyr))) +
          scale_y_continuous(name='Prevalence of intermittence \nwithin bin (% of length)') +
          annotation_logticks(sides='b') +
          theme_classic() +
          theme(text=element_text(size=12),
                axis.title.x = element_text(vjust=-2))
        
      },  error=function(cond) {
        message(cond)
        ggplot(gamsubdat, aes(x=dislogbin, y=percinter_bin)) +
          geom_bar(size=1.5, stat='identity', alpha=1/3) +
          scale_x_log10(name=bquote(
            'Binned naturalized long-term mean annual discharge'~(m^3~s^-1)),
            limits=c(0.001, 10000)) +
          scale_y_continuous(name='Prevalence of intermittence within bin (% of length)',
                             limits=c(0,1), expand=c(0,0)) +
          annotation_logticks() +
          theme_classic()
      }
    )
  }
  
  
  if (interactive) {
    #------ Test approach on the world in aggregate ---------
    plot_basIRESgam(dt_format = netsub_binir)
    
    #------ Test approach on a single basin ------
    check <- netsub_basbinir[, length(unique(percinter_bin)), by=grouping_var] %>%
      setorder(-V1)
    check[V1>1,]
    #checkpfaf = '18'
    #netsub_binir[get(grouping_var) == checkpfaf,]
    #plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id=checkpfaf)
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id= "35_6")
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id= "45_17")
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id= "14_17")
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id= "91_5")
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id= "12_7")
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id="62_18")
    plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id="14_17")
  }
  
  plot_62_18 <- plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id="62_18")
  plot_14_17 <- plot_basIRESgam(dt_format = netsub_basbinir, pfaf_id="14_17")
  #dt = netsub_basbinir[get(grouping_var) == checkpfaf,]
  
  
  #----------- Implement model across all climates/basins  -----------------------------------
  #Estimate IRES prevalence for smallest log bins by basin using GAM model
  gambasires_bin <- lapply(
    netsub_basbinir[!is.na(get(eval(grouping_var))), unique(get(grouping_var))], function(pfaf_id) {
      print(pfaf_id)
      fit_basIRESgam(dt=netsub_basbinir[get(eval(grouping_var)) == pfaf_id,],
                     verbose = F, onlyextra = T) %>%
        .[, (grouping_var) := pfaf_id]
    }) %>%
    rbindlist(fill=T)
  
  #Compute IRES prevalence based on GAM output for log bins below minimum discharge cutoff
  gambasires_bin_merge <- merge(gambasires_bin, in_extranet,
                                by.x = c(grouping_var, 'dislogbin'),
                                by.y = c(grouping_var, 'dispred'),
                                all.x=T) %>%
    setorderv(cols=c(grouping_var, 'dislogbin'))
  
  gambasires_bin_merge[
    , binL_pred := fifelse(
      dislogbin < 0.09,
      cumL_pred - data.table::shift(cumL_pred, n=1L, type='lead'),
      cumL_pred - cumL_cutoffref)]
  
  gambasires <- gambasires_bin_merge[
    , list(percinter_gam_belowcutoff =
             sum(percinter_bin_adjust*binL_pred)/sum(binL_pred)),
    by = grouping_var]
  
  #Conservatively compute IRES prevalence assuming plateauing below minimum cutoff
  #(based on last log bin
  #e.g. for a cutoff of 0.1, use IRES prevalence from reaches with discharge [0.1,0.2) )
  binlog_alt <- function(x) {
    #Mixed linear-log10 scale NOT USED
    round(x/10^floor(log10(x)))*10^floor(log10(x))
  }
  
  netsub_alt <- in_rivpred %>%
    .[dis_m3_pyr >= min_cutoff & INLAKEPERC < 1,] %>% #Exclude reaches below discharge cutoff and within lakes
    .[, `:=`(PFAF_IDclz = paste0(floor(PFAF_ID05/1000), '_', clz_cl_cmj), #Get PFAF_ID for HydroBASINS level 2
             dislogbin = binlog_alt(dis_m3_pyr), #Get discharge log bins e.g. [0.1,0.2) ; [0.2,0.3) ... [1,2) ; [2,3) ... [10) ...
             LENGTH_KM_NOLAKE = LENGTH_KM*(1-INLAKEPERC)
    )]
  
  netsub_IRESconservative <- netsub_alt[
    dis_m3_pyr >= min_cutoff & INLAKEPERC < 1,
    list(percinter_ref_overcutoff = sum(get(valuevar)*LENGTH_KM_NOLAKE)/sum(LENGTH_KM_NOLAKE),
         percinter_ref_lastbin = .SD[dislogbin == min(dislogbin),
                                     sum(get(valuevar)*LENGTH_KM_NOLAKE)/sum(LENGTH_KM_NOLAKE)]
    ),
    by=grouping_var]
  
  #Merge conservative estimates with GAM estimates
  netsub_IRES <- merge(netsub_IRESconservative, gambasires, by=grouping_var)
  
  #netsub[percinter>0 & percinter < 99, sum(cumL)]/netsub[, sum(cumL)]
  
  extraIRES_pred <- merge(in_extranet[dispred == 0.01],
                          netsub_IRES, by=grouping_var)
  
  #----------- Compute global prevalence according to both methods  ------------
  extraIRES_pred[
    , `:=`(
      percinter_all_conservative =
        (cumL_predextra*percinter_ref_lastbin + cumL_cutoffref*percinter_ref_overcutoff)/
        (cumL_predextra + cumL_cutoffref),
      percinter_all_GAM =
        (cumL_predextra*percinter_gam_belowcutoff + cumL_cutoffref*percinter_ref_overcutoff)/
        (cumL_predextra + cumL_cutoffref)
    )]
  
  check <- extraIRES_pred[(percinter_gam_belowcutoff < percinter_ref_lastbin) &
                            (percinter_gam_belowcutoff < 0.99),]
  
  if (grouping_var != 'PFA_IDclz') {
    extraIRES_pred[, clz_cl_cmj := gsub('(([0-9]*)|(NA))_', '', PFAF_IDclz)]
  }
  
  ggplot(extraIRES_pred, aes(x=percinter_ref_overcutoff)) +
    geom_point(aes(y=percinter_all_conservative), color='black') +
    geom_smooth(aes(y=percinter_all_conservative), color='black') +
    geom_point(aes(y=percinter_all_GAM), color='red') +
    geom_smooth(aes(y=percinter_all_GAM), color='red') +
    geom_abline() +
    scale_x_continuous(expand=c(0,0),
                       name='RF-predicted prevalence of IRES >= 0.1 m3/s') +
    scale_y_continuous(expand = c(0,0),
                       name = 'Extrapolated (<0.1) & RF-predicted (>= 0.1) prevalence of IRES >= 0.01 m3/s')
  
  print(paste0(
    'Assuming that rivers < 0.1 m3/s are as intermittent as those [0.1, 0.2), we predict that ',
    round(100*extraIRES_pred[
      , sum(percinter_all_conservative*cumL_pred, na.rm=T)/sum(cumL_pred, na.rm=T)]),
    '% of rivers >= 0.01 m3/s are intermittent'
  ))
  
  print(paste0(
    'Statistically extrapolating the prevalence of intermittence in rivers < 0.1 m3/s, we predict that ',
    round(100*extraIRES_pred[
      , sum(percinter_all_GAM*cumL_pred, na.rm=T)/sum(cumL_pred, na.rm=T)]),
    '% of rivers >= 0.01 m3/s are intermittent'
  ))
  
  return(list(preds=extraIRES_pred,
              plot_62_18 = plot_62_18,
              plot_14_17 = plot_14_17)
  )
}


#------ analyze_benchmark -----------------
#' Analyze benchmark
#'
#' Evaluate the performance of multiple machine learning models.
#'
#' @param in_bm \link[mlr3]{BenchmarkResult} or file name of one stored on disk.
#' @param in_measure \link[mlr3]{MeasureClassif} or \link[mlr3]{MeasureRegr} to use in model comparison
#' @param inp_resid (character) path to results directory where in_bm is stored
#' 
#' @details this function is useful to optimize the flow intermittence probability 
#' threshold above which to classify reaches as non-perennial (based on a specific 
#' performance measure of interest).
#' 
#' @return list containing:
#' *a ggplot showing sensitivity, specificity, misclassification for different probability thresholds
#' *a ggplot boxplot showing the difference in performance measure between models
#' *a data.table with, for each model, the probability threshold value at which 
#' the performance measure is maximized (and the corresponding performance measure value).
#' 
#' 
#' @export
analyze_benchmark <- function(in_bm, in_measure, inp_resdir=NULL) {
  
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
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
#' Bin misclassification metrics
#'
#' Get misclassification error, sensitivity, specificity, precision, and other
#' classification performance statistics for different subsets (i.e., bins) of 
#' the dataset. 
#'
#' @param in_predictions a data.table of predictions for a set of CV repetitions.
#' @param in_rftuned Output from \link{selecttrain_rf}; list containing inner and 
#' outer resampling results + task.
#' @param in_gaugestats data.table with additional columns to use in binning if 
#' not included as predictor variables in model.
#' @param binvar (character) column that will be used to define bins.
#' @param binfunc (character) binning approach. One of 'manual', 'equal_length', 'equal_freq'. See \link{bin_dt}.
#' @param binarg  (numeric) binning argument, depends on binning approach. See \link{bin_dt}.
#' @param interthresh (numeric) flow intermittence probability threshold above which
#'  to classify records as non-perennial
#' @param spatial_rsp whether to use results from spatial (TRUE) or non-spatial (FALSE) 
#' cross-validation (only used if \code{in_rftuned} is provided, not \code{in_predictions}). 
#' @param rspcol column to use that contains predicted probabilities (only used if
#'  \code{in_predictions} is provided, not \code{in_rftuned}).
#' 
#' @details This function was used for internal assessment of model performance.
#' Right now, the options to change binvar and binfunc are not working 
#' (need some adjustment); only manual binning based on long-term mean annual discharge ('dis_m3_pyr') is possible.
#' Update instances of 'dis_m3_pyr' for changing. 
#' 
#' @return data.table with performance statistics by subset.
#' 
#' 
#' @export
bin_misclass <-  function(in_predictions=NULL, in_resampleresult=NULL,
                          in_gaugestats=NULL, binvar, binfunc, binarg,
                          interthresh=0.5, spatial_rsp=NULL, rspcol=NULL) {

  if (!is.null(in_resampleresult)) {
    if (inherits(in_resampleresult, 'list')) {
      in_resampleresult <- get_outerrsmp(in_resampleresult, spatial_rsp=spatial_rsp)
    }
    
    in_predictions <- in_resampleresult$prediction() %>%
      as.data.table
    
    #Get average predicted probability and response across CV repeats
    predsformat <- in_predictions[, list(prob.1 = mean(prob.1)),
                                  by=.(row_id, truth)] %>%
      .[, response := fifelse(prob.1 >= interthresh, 1, 0)]
    
    #Get discharge data
    rspdat <- in_resampleresult$task$data()
    if (is.null(in_gaugestats) & (binvar %in% names(rspdat))) {
      in_gaugestats <- rspdat
    }
    
    #Get bins
    bin_gauges <- bin_dt(in_dt = in_gaugestats, binvar = 'dis_m3_pyr',
                         binfunc = 'manual', binarg=binarg,
                         bintrans=NULL, na.rm=FALSE)%>%
      .[, row_id := .I]
    
  } else if (is.null(in_resampleresult) & !is.null(in_predictions)) {
    #Get bins
    bin_gauges <- bin_dt(in_dt = in_predictions, binvar = 'dis_m3_pyr',
                         binfunc = 'manual', binarg=binarg,
                         bintrans=NULL, na.rm=FALSE) %>%
      .[, row_id := .I]
    
    if (!is.null(spatial_rsp) & is.null(rspcol)) {
      if (spatial_rsp) {
        rspcol <- 'IRpredcat_CVsp'
      } else {
        rspcol <- 'IRpredcat_CVnosp'
      }
    }
    
    predsformat <- in_predictions[, `:=`(truth = intermittent_o1800,
                                         response = get(rspcol),
                                         row_id = .I)]
  }
  
  #Merge bins with predictions
  gselpreds <- merge(bin_gauges, predsformat, by='row_id', all.y=F)
  
  #Function to get statistics for each bin
  
  
  #Run function for each bin
  performtable <- gselpreds[, compute_confustats(.SD),
                            by=.(bin, paste0(round(bin_lmin, 2), '-',
                                             round(bin_lmax, 2)))] %>%
    setorder(bin) %>%
    rbind(compute_confustats(gselpreds)[, paste0 := 'All'], use.names=T, fill=T)
  
  return(performtable)
}

#------ ggvimp -----------------
#' Plot of variable importance
#'
#' Produce a bar ggplot of variable importance for random forest model.
#' 
#' @param in_rftuned Output from \link{selecttrain_rf}; list containing inner and 
#' outer resampling results + task.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param varnum number of variables whose variable importance to plot.
#' @param spatial_rsp (boolean) whether to use variable importance values derived from
#' spatial (TRUE) or non-spatial (FALSE) cross-validation.
#' 
#' @details this function is used to produce Figure 2 in the Main Text of Messager et al. 2021. 
#' Additional formatting is needed to exactly match Fig. 2.
#' 
#' @return ggplot
#' 
#' 
#' @export
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
  outp <- ggplot(varimp_basic[1:varnum,],aes(x=varname,
                                             color =Category, fill=Category)) +
    geom_bar(aes(y=imp_wmean), stat = 'identity', alpha=0.7) +
    geom_errorbar(aes(ymin=imp_wmean-imp_wsd, ymax=imp_wmean+imp_wsd)) +
    scale_x_discrete(labels = function(x) {
      stringr::str_wrap(tolower(x), width = 27)
    },
    limits=rev) +
    scale_fill_manual(values=c('#fdb462','#80b1d3','#b3de69','#bc80bd','#696868'),
                      drop=FALSE) +
    scale_color_manual(values=c('#fdb462','#80b1d3','#b3de69','#bc80bd','#696868'),
                       drop=FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(size=8),
          axis.title.x  = element_blank(),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.position = c(0.80, 0.5)
    ) +
    scale_y_continuous(expand=c(0,0), position = 'right') +
    coord_flip(ylim=c(0, min(varimp_basic[, max(imp_wmean+imp_wsd)+1], 100)),
               clip='off')
  
  return(outp)
}

#------ ggpartialdep -----------------
#' Plot of partial dependence
#'
#' Produce line or raster ggplots of univariate or bivariate partial dependence.
#' i.e., estimates of the marginal relationship between predictor variables and 
#' the model’s predictions (probability of intermittence) by holding the rest of 
#' the predictors at their respective mean values. Bivariate plots show the co-linearity 
#' in response between two predictors.
#' 
#' @param in_rftuned Output from \link{selecttrain_rf}; list containing inner and 
#' outer resampling results + task.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param colnums number of variables to include, in decreasing order of variable importance.
#' @param ngrid number of predictor variable values to check model's marginally predicted value for.
#' @param nodupli whether to include variable types only once (i.e. if minimum discharge is in, not including mean discharge)
#' @param nvariate (1 or 2) whether to analyze univariate (1) or bivariate (2) 
#' partial dependence.
#' @param parallel (boolean) whether to compute function in parrallel 
#' (as it can be computationally intensive).
#' @param spatial_rsp (boolean) whether to use outputs from spatial (TRUE) or 
#' non-spatial (FALSE) cross-validation.
#' 
#' @details this function is used to produce Figure S5 in the Supplementary Information for Messager et al. 2021. 
#' This function was initially developed for bivariate plots but only univariate plots 
#' were produced for the final manuscript so there may be obsolete snippets left.
#' 
#' @return pages of gridded plots with each page containing 9 plots
#' 
#' 
#' @export
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
#' Plots of Intermittence Prediction Residuals (IPR)
#'
#' Plots relating Intermittence Prediction Residuals (IPR) for gauging stations 
#' to a variety of metadata for the gauging stations including the number of years
#' on record, the average flow intermittence duration, the aridity in the catchment
#' of the station, etc.
#' 
#' @param in_gpredsdt model predictions and environmental characteristics of gauging stations. 
#' Here, output from \link{bind_gaugepreds}.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param spatial_rsp (boolean) whether to use outputs from spatial (TRUE) or 
#' non-spatial (FALSE) cross-validation.
#' @param interthresh (numeric) flow intermittence probability threshold above which
#'  to classify records as non-perennial
#' @param yearthresh  (integer) minimum year from which to analyze/consider discharge record.
#' 
#' @details this function is used to produce Extended Data Fig. 2b-e. Not all
#' plots were used in the final figure though.
#' 
#' @return gtable of plots
#' 
#' 
#' @export
gggaugeIPR <- function(in_gpredsdt, in_predvars, spatial_rsp = FALSE,
                       interthresh = 0.5, yearthresh) {
  
  #Capitlize first letter of string, convert to lower case all other letters.
  SenCapstr <- function(y) {
    paste0(toupper(substring(y, 1,1)),
           tolower(substring(y, 2, 100)))
  }
  
  #Plot numeric variables
  in_gpredsdt[, DApercsdiff := abs(
    (area_correct-UPLAND_SKM)/mean(c(area_correct, UPLAND_SKM))), by=GAUGE_NO]
  
  if (spatial_rsp) {
    predcol <- 'IRpredprob_CVsp'
    uncertcol <- 'preduncert_CVsp'
  } else {
    predcol <- 'IRpredprob_CVnosp'
    uncertcol <- 'preduncert_CVnosp'
  }
  
  predmelt_num <- in_gpredsdt[
    , which(as.vector(unlist(lapply(in_gpredsdt, is.numeric)))), with=F] %>%
    cbind(in_gpredsdt[, c('GAUGE_NO', 'intermittent_o1800'), with=F]) %>%
    melt(id.vars=c('GAUGE_NO', 'intermittent_o1800', predcol, uncertcol))
  
  predmelt_num[, IPRclass := fcase(
    intermittent_o1800 == "0" & get(predcol) < 0.5, 'True: perennial | Pred: perennial',
    intermittent_o1800 == "0" & get(predcol) >= 0.5, 'True: perennial | Pred: non-perennial',
    intermittent_o1800 == "1" & get(predcol) < 0.5, 'True: non-perennial | Pred: perennial',
    intermittent_o1800 == "1" & get(predcol) >= 0.5, 'True: non-perennial | Pred: non-perennial'
  )]
  
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
  
  varlabels[, varname := SenCapstr(varname)]
  
  #levels(predmelt_num$variable) <-  varlabels[varname %in% predmelt_num$variable, varname]
  varstoplot_continuous <- varlabels[varcode %in% c('totalYears_kept_o1800',
                                                    'mDur_o1800',
                                                    'mFreq_o1800',
                                                    'DApercsdiff',
                                                    'ari_ix_uav'),]
  
  varstoplot_log <- varlabels[varcode %in% c('UPLAND_SKM',
                                             'dis_m3_pyr',
                                             'dor_pc_pva'),]
  
  colorpal_continuous = c("#FFAA00", "#a80000", "#004DA8", "#7AB6F5")
  
  
  #plotdt = predmelt_num[variable %in% "DApercdiff",]
  
  ggIPR_util <- function(plotdt, logx = FALSE, legend=FALSE) {
    p <- ggplot(plotdt) +
      # geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
      #                            ymin=ymin, ymax=ymax, fill=fillpal),
      #           alpha=1/4) +
      geom_point(aes(x=value, y=get(uncertcol), color=IPRclass),
                 alpha = 1/3) +
      scale_y_continuous(limits=c(-1, 1), expand=c(0,0)) +
      # scale_fill_manual(values=colorpal,
      #                   name='Predicted regime',
      #                   labels = c('Perennial', 'Intermittent')) +
      
      scale_color_manual(values=colorpal_continuous) +
      labs(x='Value', y='Intermittence Prediction Residuals (IPR)') +
      geom_hline(yintercept=0, alpha=1/2) +
      #new_scale_color() +
      # geom_smooth(aes(x=value, y=get(uncertcol), color=intermittent_o1800),
      #             method='gam', formula = y ~ s(x, k=3)) +
      # annotate("text", x = Inf-5, y = 0.5, angle = 90,
      #          label = "Pred:Int, Obs:Per",
      #          color = colorpal[[1]]) +
      # annotate("text", x = Inf-5, y = -0.5, angle = 90,
      #          label = "Pred:Per, Obs:Int",
      #          color = colorpal[[2]]) +
      #scale_x_sqrt(expand=c(0,0)) +
      coord_cartesian(expand=FALSE, clip='off') +
      facet_wrap(~variable, scales='free', ncol=1,
                 labeller=as_labeller(varlabels)) +
      theme_classic() +
      theme(legend.position = c(0.85, 0.2),
            legend.title = element_blank())
    
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
  predmelt_cat <- in_gpredsdt[, c('GAUGE_NO', 'intermittent_o1800', uncertcol,
                                  'ENDORHEIC', 'clz_cl_cmj'), with=F] %>%
    .[, `:=`(intermittent_o1800 = factor(intermittent_o1800, levels=c('0','1')))] %>%
    melt(id.vars=c('GAUGE_NO', 'intermittent_o1800', uncertcol))
  levels(predmelt_cat$variable) <- c('Endorheic',
                                     'Climate Zone (catchment majority)')
  
  
  colorpal <- c('#1f78b4', '#ff7f00')
  rectdf_cat <- data.table(
    xmin=rep(-Inf, 4),
    xmax=rep(Inf, 4),
    ymin=c(-1, interthresh-1, 0, interthresh),
    ymax=c(interthresh-1, 0, interthresh, 1),
    fillpal = rep(colorpal, 2)
  )
  gaugeIPR_catplot <-
    ggplot(predmelt_cat) +
    geom_rect(data=rectdf_cat, aes(xmin=xmin, xmax=xmax,
                                   ymin=ymin, ymax=ymax, fill=fillpal),
              alpha=1/4) +
    scale_fill_manual(values=colorpal,
                      name='Predicted regime',
                      labels = c('Perennial', 'Intermittent')) +
    new_scale_fill() +
    geom_boxplot(aes(x=as.factor(value), y=get(uncertcol),
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
  
  
  # outp <- grid.arrange(p_continuous, p_log, gaugeIPR_catplot,
  #                      layout_matrix = rbind(c(1,1,1),
  #                                            c(1,1,1),
  #                                            c(2,2,2),
  #                                            c(3,3,3))
  # )
  
  outp <- grid.arrange(p_continuous, ncol=1)
  
  
  return(outp)
}

#------ krige_spgaugeIPR----
#Ignore, not used in final analysis.
krige_spgaugeIPR <- function(in_gpredsdt, in_gaugep,
                             kcutoff=100000, inp_bufrasdir=NULL,
                             overwrite = FALSE) {
  
  predsp_gaugep <- st_as_sf(in_gpredsdt)
  
  # bufras_vec <- file.path(inp_bufrasdir,
  #                         grep('bufras_T.*proj[.]tif$',
  #                              list.files(inp_bufrasdir),
  #                              value=T)
  # )
  
  #Convert to SpatialPointDataFrame to use in gstats and project to Goode Homolosine
  crs_aeqd <- "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  predsp_gaugep_df <- st_transform(predsp_gaugep, crs=crs_aeqd)  %>%
    as_Spatial()
  
  #Possible models
  varmods <- as.character(vgm()[,'short'])
  
  in_sf = predsp_gaugep_df
  ycol = 'intermittent_o1800'
  in_varmods = varmods
  
  customfit_variogram <- function(in_sf, ycol, in_varmods, kcutoff) {
    # Compute the sample variogram; note that the f.1 trend model is one of the
    # parameters passed to variogram(). This tells the function to create the
    # variogram on the de-trended data.
    var_out <- gstat::variogram(intermittent_o1800 ~ 1, in_sf,
                                cutoff=kcutoff, width=500, cressie=TRUE) #cloud=T)
    
    # Compute the variogram model by passing the nugget, sill and range values
    # to fit.variogram() via the vgm() function.
    fit_out  <- fit.variogram(var_out, fit.ranges = TRUE, fit.sills = T,
                              vgm(in_varmods[!(varmods %in% c('Pow', 'Int'))]),
                              fit.kappa = TRUE)
    
    plot <- ggplot(var_out, aes(x=dist, y=gamma)) +
      geom_point(alpha=1/3) +
      geom_line(data = variogramLine(fit_out, maxdist = kcutoff)) +
      coord_cartesian(expand=T) +
      scale_y_log10() +
      scale_x_sqrt(breaks=c(0,100,1000,10000,25000,50000, 100000),
                   labels=c(0,100,1000,10000,25000,50000, 100000)) +
      theme_classic()
    
    return(list(
      plot = plot,
      mod = fit_out
    ))
    
  }
  
  bas = names(table(predsp_gaugep_df$PFAF_ID03) > 20)[8]
  lapply(names(table(predsp_gaugep_df$PFAF_ID03) > 20),
         function(bas) {
           subdat <- subset(predsp_gaugep_df, PFAF_ID03 == bas)
           
           
           fit_out <- customfit_variogram(in_sf = subdat, ycol='intermittent_o1800',
                                          in_varmods = varmods, kcutoff=kcutoff)
           fit_out$plot
           
           
         })
  
  
  
  ################### Kriging on fit data #################################
  # Compute the sample variogram; note that the f.1 trend model is one of the
  # parameters passed to variogram(). This tells the function to create the
  # variogram on the de-trended data.
  var_smpl <- gstat::variogram(IRpredcat_CVsp ~ 1, predsp_gaugep_df,
                               cutoff=kcutoff, width=500, cressie=TRUE) #cloud=T)
  
  ggplot(var_smpl, aes(x=dist, y=gamma)) +
    geom_point(alpha=1/3) +
    coord_cartesian(expand=T) +
    scale_y_log10() +
    theme_classic()
  
  # scale_x_sqrt(breaks=c(0,100,1000,10000,25000,50000, 100000),
  #              labels=c(0,100,1000,10000,25000,50000, 100000)) +
  
  # Compute the variogram model by passing the nugget, sill and range values
  # to fit.variogram() via the vgm() function.
  varmods <- as.character(vgm()[,'short'])
  dat_fit  <- fit.variogram(var_smpl, fit.ranges = TRUE, fit.sills = T,
                            vgm(varmods[!(varmods %in% c('Pow', 'Int'))]),
                            fit.kappa = TRUE)
  # The following plot allows us to assess the fit
  plot(var_smpl, dat_fit)
  attr(dat_fit, "SSErr")
  
  
  ################### Kriging on residuals ####################################
  # Compute the sample variogram; note that the f.1 trend model is one of the
  # parameters passed to variogram(). This tells the function to create the
  # variogram on the de-trended data.
  var_smpl <- gstat::variogram(preduncert_CVsp ~ 1, predsp_gaugep_df,
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
  dat_fit  <- fit.variogram(var_smpl, fit.ranges = TRUE, fit.sills = T,
                            vgm(varmods[!(varmods %in% c('Pow', 'Int'))]),
                            fit.kappa = TRUE)
  # The following plot allows us to assess the fit
  plot(var_smpl, dat_fit)
  
  # Import gauge buffer mask
  
  
  #Predict error
  # kmod <- gstat(formula=IPR~1, locations=predsp_gaugep_df, model=dat_fit,
  #               maxdist = kcutoff)
  # kpl <-  lapply(seq_along(bufras_vec), function(i) {
  #   bufmask <- raster(bufras_vec[[i]]) %>%
  #     as('SpatialGrid')
  #   kp <- round(raster(predict(kmod, bufmask)))
  #
  #   outras = file.path(dirname(bufras_vec[[i]]),
  #                      gsub('bufras', 'krigpred', basename(bufras_vec[[i]])))
  #
  #
  #   print(paste0('Writing ', outras, '...'))
  #   writeRaster(kp, outras, datatype = "INT2S", overwrite=overwrite)
  #   return(outras)
  # }) %>%
  #   do.call(rbind, .)
  
  return(kpl)
}
#------ mosaic_kriging -------------
#Ignore, not used in final analysis.
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


#------ map_basinBACC ------
#' Map prediction performance at the basin level
#'
#' Compute the number of gauging stations, bias, residual spatial autocorrelation
#' and balanced classification accuracy for each HydroBASIN level 3.
#'
#' @param in_gaugepred data.table of predictions for gauging stations. Output from \link{bind_gaugepreds}.
#' @param in_rivernetwork river network with model predictions and other attributes. Output from \link{netpredformat}.
#' @param inp_basin (character) absolute path to feature class of HydroBASINS level 3 polygons
#' @param outp_basinerror (character) absolute path of output basin polygons with computed statistics
#' @param spatial_rsp (boolean) whether to use outputs from spatial (TRUE) or 
#' non-spatial (FALSE) cross-validation.
#'
#' @details This function was used to produce the data underlying Figure 3 of the Main Text of Messager et al. 2021.
#' see \link{test_joincount} for more details on the procedure for testing spatial autocorrelation.
#' 
#' @return list containing a data.table and two plots.
#' The data.table contains an estimate 
#' 
#' 
#' @export
map_basinBACC <- function(in_gaugepred, #rfpreds_gauges,
                          in_rivernetwork,
                          inp_basin,
                          outp_basinerror,
                          spatial_rsp) {
  bas03 <- st_read(dsn = dirname(inp_basin),
                   layer = basename(inp_basin))
  
  if (spatial_rsp) {
    predcol <- 'IRpredcat_CVsp'
  } else {
    predcol <- 'IRpredcat_CVnosp'
  }
  
  classcol='intermittent_o1800'
  dat <- in_gaugepred[, .N, by=classcol]%>%
    setnames(old=classcol, new='classcol')
  minoratio <- dat %>%
    setorder(N) %>%
    .[, list(minoclass=classcol[1], ratio=N[2]/N[1])]
  
  in_gaugepred[, ':='(classweights = fifelse(intermittent_o1800 == minoratio$minoclass,
                                             minoratio$ratio, 1)
  )]
  
  in_rivernetwork[, PFAF_ID03 := substr(PFAF_ID05, 1, 3)]
  gnetjoin <- merge(in_gaugepred,
                    in_rivernetwork[,.(HYRIV_ID, PFAF_ID03)],
                    by='HYRIV_ID', all.x=T, all.y=F)
  
  
  
  #Join count test for each basin with more than 20 gauges
  jclist <- test_joincount(gnetjoin)
  
  #Compute statistics by basin
  basbacc <- gnetjoin[, list(
    basbacc =  sum((intermittent_o1800==get(predcol)) *
                     classweights) / sum(classweights),
    acc = sum(intermittent_o1800==get(predcol))/.N,
    gnum = .N,
    percinter = round((100*sum(intermittent_o1800=='1')/.N), 2),
    predbias = (sum(get(predcol)==1)/.N) - (sum(intermittent_o1800=='1')/.N)
  ), by=PFAF_ID03] %>%
    merge(jclist, by='PFAF_ID03', all.x=T) %>%
    .[, `:=`(stdeviate_diff = stdeviate_pred - stdeviate_ref,
             stdeviate_ratio = stdeviate_pred/stdeviate_ref
    )]
  
  #Investigate how join count statistics of predictions relate to ref data
  plot_jcdeviate <- ggplot(basbacc, aes(x=stdeviate_ref, y=stdeviate_pred)) +
    geom_point() +
    geom_abline(size=1) +
    labs(x='Standard deviate - reference flow intermittency class',
         y='Standard deviate - predicted flow intermittency class') +
    coord_fixed() +
    theme_bw()
  
  #Check relationship between accuracy and number of gauges
  palette_acc <- c('#314D8F', '#3C8A9E', '#3EB591', '#2BD93D',
                   '#8EF026', '#FCF228', '#F2C035', '#E08E46', '#CC6D5A')
  basbacc[, colacc := factor(
    fcase(
      acc >= 0.95, 1L,
      ((0.90 <= acc) & (acc < 0.95)), 2L,
      ((0.85 <= acc) & (acc < 0.90)), 3L,
      ((0.80 <= acc) & (acc < 0.85)), 4L,
      ((0.70 <= acc) & (acc < 0.80)), 5L,
      ((0.60 <= acc) & (acc < 0.70)), 6L,
      ((0.50 <= acc) & (acc < 0.60)), 7L,
      ((0.005 <= acc) & (acc < 0.50)), 8L,
      acc < 0.005, 9L
    )
  ), by=.I]
  
  plot_acc <- ggplot(basbacc[, .N, by=.(gnum, acc, colacc)],
                     aes(x=gnum, y=100*acc, size=as.numeric(N), fill=colacc)) +
    geom_point(colour="black",pch=21) +
    #geom_quantile() +
    scale_size_continuous(name='# of basins', breaks=c(1, 2, 5, 10)) +
    scale_fill_manual(values=palette_acc, guide = FALSE) +
    scale_x_log10(breaks=c(1, 10, 100, 500)) +
    labs(x='# of gauging stations', y='Accuracy (%)') +
    theme_classic() +
    theme(legend.position = c(0.85, 0.2),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.height=unit(0.75,"line")
    )
  
  
  #Check relationship between bias and number of gauges
  palette_bias <- c('#a8273a', '#d46e4c', '#ffac70',
                    '#76c284', '#00ffff', '#3d7294', '#374559')
  basbacc[, colbias := factor(
    fcase(
      predbias >= 0.50, 1L,
      ((0.20 <= predbias) & (predbias < 0.50)), 2L,
      ((0.10 <= predbias) & (predbias < 0.20)), 3L,
      ((-0.10 <= predbias) & (predbias < 0.10)), 4L,
      ((-0.20 <= predbias) & (predbias < -0.10)), 5L,
      ((-0.50 <= predbias) & (predbias < -0.20)), 6L,
      predbias < -0.50, 7L
    )
  ), by=.I]
  
  plot_bias <-ggplot(basbacc[, .N, by=.(gnum, predbias, colbias)],
                     aes(x=gnum, y=100*predbias, size=as.numeric(N), fill=colbias)) +
    geom_hline(yintercept=0, alpha=1/2) +
    geom_point(colour="black",pch=21) +
    scale_size_continuous(name='# of basins', breaks=c(1, 2, 5, 10)) +
    scale_fill_manual(values=palette_bias, guide = FALSE) +
    scale_x_log10(breaks=c(1, 10, 100, 500)) +
    labs(x='# of gauging stations', y='Bias (%)') +
    theme_classic() +
    theme(legend.position = c(0.85, 0.8),
          legend.background = element_blank())
  
  
  #Output stats to .gpkg
  basbacc_format <- base::merge(
    bas03, basbacc,
    by.x='PFAF_ID', by.y='PFAF_ID03',
    all.x.=T, all.y=T)
  basbacc_format$gdens <- with(basbacc_format, gnum/UP_AREA)
  
  st_write(obj=basbacc_format,
           dsn=outp_basinerror,
           driver = 'gpkg',
           delete_dsn=F)
  
  return(list(
    plot_acc = plot_acc,
    plot_bias = plot_bias,
    plot_jcdeviate = plot_jcdeviate
  ))
}

##### -------------------- Report functions -----------------------------------
#------ get_basemapswintri ------------
#' Get and project basemaps
#'
#' Get and project basemaps into Winkel Tripel for subsequent mapping of gauging stations.
#' 
#' @return list containing basemaps and graticule in Winkel Tripel projection
#' 
#' @export
get_basemapswintri <- function() {
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"
  
  wcountries <- rnaturalearth::ne_countries(
    scale = "medium", returnclass = "sf") %>%
    sfformat_wintri
  wland <- rnaturalearth::ne_download(
    scale = 110, type = 'land', category = 'physical', returnclass = "sf") %>%
    sfformat_wintri
  wlakes <- rnaturalearth::ne_download(
    scale = 110, type = 'lakes', category = 'physical', returnclass = "sf") %>%
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
#' Plot river gauging stations basemaps
#' 
#' Utility function to format plot of basic basemap to show the distribution of river gauging stations.
#' 
#' @param in_basemaps list of basemaps \link[sf][sf] objects.
#' 
#' @details does not include rivers in basemaps. Too long to draw, too messy/
#' 
#' @return ggplot
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
#' Map river gauging stations
#'
#' Make two maps of streamflow gauging stations used in model training and testing.
#' One map of perennial stations and the other one of non-perennial stations. 
#'
#' @param in_gaugepred data.table of predictions for gauging stations. Output from \link{bind_gaugepreds}.
#' @param in_basemaps list of basemaps \link[sf][sf] objects.
#' @param binarg (vector) inner bin breaks to plot histogram for gauging stations.
#' @param binvar (character) column to plot histogram of gauging stations with
#'
#' @details 
#' The maps include a histogram of streamflow record duration in each map.
#' binarg c(30, 60, 100) is supplied to show four bins in the histogram (10-30, 30-60, 60-100, >100).
#' 
#' This function was used to produce the data underlying Figure 3 of the Main Text of Messager et al. 2021.
#' see \link{test_joincount} for more details on the procedure for testing spatial autocorrelation.
#' 
#' @return ggplot (patchwork)
#' 
#' 
#' @export
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
      geom_vline(xintercept = meanpos - 0.5) +
      annotate(geom='text', angle=90,
               x=meanpos-0.6,
               y=0.5*max(as.data.table(in_gdf)[, .N, by=bin]$N),
               label=paste(round(permean), 'y')) +
      scale_fill_manual(values=cs) +
      scale_x_continuous(name = 'Years of data',
                         breaks = c(binarg_tick, x_tick),
                         labels = c(rep(c(""), len), c(l_freq, u_freq[bins]))) +
      scale_y_continuous(name = 'Count (gauges)')+
      coord_cartesian(clip='off', expand=c(0,0)) +
      theme_classic() +
      theme(legend.position = 'none',
            text = element_text(size=9),
            plot.background = element_blank(),
            panel.background = element_blank(),
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
  outp <- p_patch + plot_annotation(tag_levels = 'a')
  return(outp)
}

#------ comp_GRDCqstats ----------
#' Compute GRDC discharge statistics.
#'
#' Compute number of years with discharge data,  qmean, mean minimum flow,
#' mean minimum 3-day average flow, q10, q90 and q99 for a GRDC gauging station
#'
#' @param path (character) full path to the formatted daily streamflow record for a GRDC gauging station.
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param minyear (integer) start year to include in computing statistics.
#' @param maxyear (integer) last year to include in computing statistics.
#' @param verbose (boolean) whether to print the input path upon executing the function.
#' 
#' @return single-row data.table 
#' 
#' 
#' @export
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
#' Compute modeled discharge statistics for GRDC gauging stations
#'
#' Compute performance statistics (regression R2 and sMAPE) for WaterGAP v2.2 discharge 
#' predictions of long-term mean discharge and Q90 against observed discharge 
#' at GRDC gauging stations with ≥20 years of streamflow data.
#' 
#' @param in_qstats data.table of discharge statistics for GRDC streamflow gauging stations. 
#' Each row corresponds to a single station. Output from \link{comp_GRDCqstats}.
#' @param in_selgauges data.table or data.frame of gauging stations to analyze. 
#' @param binarg discharge bin limits to divide performance statistics asssessment by 
#' (to create size classes based on long-term mean annual flow in m3/s).
#' 
#' 
#' @details gauging stations to analyze \code{in_selgauges}  are further subsetted
#' to keep only those with at least 20 years of daily discharge data.
#' 
#' This function performs a log-log regression across all gauges and non-log regressions
#' for each size class to compute R2.   
#'  
#'  This function was used to produce Table S1 in the Supplementary Information of Messager et al. 2021
#'  
#' @return data.table of performance statistics by streamflow size class for 
#' mean annual flow and Q90.
#' 
#' @export
eval_watergap <- function(in_qstats, in_selgauges, binarg) {
  qsub <- in_qstats[!is.na(GRDC_NO) &
                      nyears >= 20 &
                      GRDC_NO %in% in_selgauges$GRDC_NO,] %>%
    .[!(GRDC_NO %in% c('5204010')),]
  #5204010 - issue with units which does not affect zero flow assessment
  #Checked 1160520 - mixed but not clear enough
  
  qsubp <- merge(qsub, in_selgauges, by='GRDC_NO', all.x=T, all.y=F) %>%
    .[dor_pc_pva < 100,]
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
#' Format plot scales
#'
#' Utility function to format plot scales for density distribution plots of environmental variables.
#' 
#' @param in_df data.frame with all records for environmental variables. 
#' Used to determine the appropriate range of values for each variable. In this case, 
#' a data.table of the river network hydro-environmental attributes.
#' @param varstopplot vector of variable names that will be plots and for which 
#' to return a list of scales
#'  
#' @return list of x and y scale objects + cartesian coordinates for ggplot
#' 
#' 
#' @export
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
#' Plot of environmental histogram
#'
#' Utility function to create an individual density plot of the distribution of a 
#' given environmental variables across gauges and the whole global river network.
#' 
#' @param vartoplot (column) variable for which to produce a density plot.
#' @param in_gaugedt data.table of gauging stations' environmental attributes.
#' @param in_rivdt data.table of global river network's environmental attributes 
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param scalesenvhist list of scale objects to format plot. From \link{formatscales}.
#' 
#' @return ggplot with two density distributions of the environmental variable,
#' one for the gauging stations and one for the global river network. 
#' 
#' 
#' @export
ggenvhist <- function(vartoplot, in_gaugedt, in_rivdt, in_predvars,
                      scalesenvhist) {
  print(vartoplot)
  
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
#' Layout plots of environmental histograms
#'
#' Run plotting functions across predictor variables and arrange plots
#' 
#' @param in_rivernetwork data.table of global river network's environmental attributes. Here, output from \link{rformat_network}.
#' @param in_gaugepred selected gauging stations' environmental attributes. Here, output from \link{write_gaugepreds}.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' 
#' @details function used to produce Extended Data Fig. 8  c-p in Messager et al. 2021.
#' 
#' @return ggplots with two density distributions of the environmental variable,
#' one for the gauging stations and one for the global river network. 
#' 
#' 
#' @export
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
#' Tabulate benchmarks
#'
#' Function to compile hyperparameters and performance statistics for a 
#' given random forest model into two summary tables.
#
#' @param in_bm \link[mlr3]{BenchmarkResult} or file name of one stored on disk.
#' @param in_bmid (character)Type of model. one of 'classif1' (selection of classification algorithm), 
#' 'regr1' (selection of regression model), 'classif2' (comparison between without and without variable selection). 
#' @param inp_resdir (character) path to results directory where in_bm is stored.
#' @param interthresh (numeric) flow intermittence probability threshold above which
#'  to classify records as non-perennial.
#'  
#' @details Function used to produce Tables S2 and Tables S3.
#' 
#' @return list of two one-row table, one with hyperparameter values and
#' the other with performance metrics.
#' 
#' 
#' @export
tabulate_benchmarks <- function(in_bm, in_bmid, inp_resdir=NULL, interthresh=NULL) {
  
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  print('Getting table content...')
  tbbm <- bm_msrtab(in_bm, interthresh=interthresh) %>%
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
                                   spe = round(spe, 3), sen = round(sen, 3),
                                   pre = round(pre, 3), bbrier=round(bbrier, 3), 
                                   auc=round(auc, 3)
  )]
  
  return(list(setup=setup_table, results=results_table))
}

#------ compute_IRpop -------
#' Compute population living nearest to an intermittent river
#'
#' Compute total population and percentage of the world's population which live
#' nearest to a river or stream reach predicted to be non-perennial.
#
#' @param in_rivpred data.table of model predictions for global river network. output from \link{netpredformat}.
#' @param inp_linkpop absolute path to table identifying number of people living 
#' nearest to each river reach in the global river network.
#' @param valuevar name of the column in \code{in_rivpred} by which to compute 
#' percentage of the population (i.e., column name that contains the predicted
#' binary flow intermittence class for each river reach)
#' 
#' @return two-row data.table with the number and percentage of people living nearest
#' to non-perennial vs. perennial river reaches.
#' 
#' @export
compute_IRpop <- function(in_rivpred, inp_linkpop, valuevar) {
  linkpop = fread(inp_linkpop) %>%
    merge(in_rivpred, by.x = 'VALUE', by.y = 'HYRIV_ID', all.x=T, all.y=T)
  linkpop[!is.na(VALUE) & is.na(dis_m3_pyr), .N]
  linkpop[is.na(VALUE) & !is.na(dis_m3_pyr), .N]
  
  #Compute number of people living nearest to a IRES
  globalestimate <- linkpop[, list(popsum=sum(SUM, na.rm=T)), by=valuevar] %>%
    .[, popperc:=popsum/sum(.$popsum)]
  
  return(globalestimate)
}

#------ formatmisclass_bm -------------
#Not used in analysis in the end
formatmisclass_bm <- function(in_bm, in_bmid, inp_resdir=NULL) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(file.path(inp_resdir, in_bm))
  }
  
  print('Getting table content...')
  thresh_dt <- threshold_dat(bmres =in_bm) %>%
    format_modelcompdat(in_bmid)
  return(thresh_dt)
}


#------ ggmisclass_bm -------------
#Not used in analysis in the end
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


#------ test_thresholdsensitivity ---------
#' Test threshold sensitivity
#'
#' Test sensitivity of model predictions to flow intermittence probability threshold
#' used to classify global river reaches into non-perennial and perennial classes.
#' 
#' @param in_gpredsdt data.table of model predictions for gauging stations. Here, output from \link{bind_gaugepreds}.
#' @param in_rivpred data.table of model predictions for global river network. Here, output from \link{netpredformat}.
#' @param threshrange_gauges numerical vector of probability threshold values for
#'  which to assess and plot model predictions and performance for gauging stations.
#' @param threshrange_network numerical vector of probability threshold values for which to produce predictions
#' of flow intermittence for the global network (set narrower than \code{thresrange_gauges} as each one takes
#' much longer to compute).
#' @param mincutoff minimum long-term mean annual flow (MAF) to include in assessment of gauges and network 
#' (i.e., mincutoff == 0.1 means that only gauging stations and river reaches with a WaterGAP estimated MAF >= 0.1 
#' are included in the sensitivity analysis)
#' @param gaugescol name of the column in \code{in_gpredsdt} containing the predicted probability of flow intermittence.
#' @param netcol name of the column in \code{in_rivpred} containing the predicted probability of flow intermittence.
#' 
#' @return list with three elements:
#' * A ggplot of the probability threshold ranges for gauges >= 10m3/s and gauges 
#' < 10 m3/s that maximize balanced accuracy, bias, raw accuracy, and |sensitivity-specificity|.
#' * A ggplot of the predicted global prevalence of IRES as a function of probability thresholds for two gauge size classes (<10 and >=10 m3/s)
#' * data.table containing the predicted global prevalence of IRES based on the range of probability threshold \code{threshrange_network}.
#'
#' @export
test_thresholdsensitivity <- function(in_gpredsdt,
                                      in_rivpred,
                                      threshrange_gauges = seq(0.30, 0.70, 0.01),
                                      threshrange_network = seq(0.45, 0.55, 0.01),
                                      mincutoff = 0.1,
                                      gaugescol = 'IRpredprob_CVnosp',
                                      netcol = 'predbasic800') {
  
  #--------------- Compute sensitivity for gauges -------------------
  gpreds_format <- in_gpredsdt[dis_m3_pyr >= mincutoff,] %>%
    .[,`:=`(truth = intermittent_o1800,
            row_id = .I)]
  
  threshgrid <- expand.grid(threshrange_gauges, threshrange_gauges) %>%
    setnames(c('thresh_u10', 'thresh_o10'))
  
  gstats <- mapply(function(tu10, to10) {
    #print(paste('thresh_u10:', tu10, 'thresh_o10:', to10))
    gpreds_format[dis_m3_pyr < 10,
                  response := as.numeric(get(gaugescol) >= tu10)]
    gpreds_format[dis_m3_pyr >= 10,
                  response := as.numeric(get(gaugescol) >= to10)]
    
    confustats <- compute_confustats(gpreds_format, ndigits=4) %>%
      .[, predratio := (as.numeric(gsub("[|][0-9]+", "", predtrue_inter)) -
                          as.numeric(gsub("[0-9]+[|]", "", predtrue_inter)))/
          as.numeric(gsub("[0-9]+[|]", "", predtrue_inter))] %>%
      cbind(gpreds_format[, list(
        bacc = mlr3measures::bacc(truth, as.factor(response)),
        threshold_u10 = tu10,
        threshold_o10 = to10
      )]
      )
    
    return(confustats)
  },
  tu10 = threshgrid$thresh_u10,
  to10 = threshgrid$thresh_o10,
  SIMPLIFY = F
  ) %>%
    do.call(rbind, .)
  
  #--------------- Plot sensitivity for gauges -------------------
  collimit_bias <- max(abs(gstats$predratio)) * c(-100, 100)
  minmax_gpoints <- gstats[c(which.min(abs(predratio)),
                             which.min(misclas),
                             which.max(bacc),
                             which.min(abs(sens-spec))),] %>%
    .[, perfstat := c('Minimum bias', 'Maximum accuracy',
                      'Maximum BACC', 'Sensitivity = Specificity')]
  
  threshpoly_gpredratio <- gstats[
    which(abs(predratio)<=min(abs(predratio))+0.01),] %>%
    .[, perfstat := 'Bias']
  threshpoly_gpredratio[, .(min(predratio), max(predratio))]
  
  threshpoly_gmisclass <- copy(gstats)[
    which(misclas<=min(misclas)+0.01),] %>%
    .[, perfstat := 'Raw accuracy']
  
  threshpoly_gmisclass[, .(1-max(misclas), 1-min(misclas))]
  
  threshpoly_gbacc <- copy(gstats)[
    which(bacc>=max(bacc)-0.01),] %>%
    .[, perfstat := 'Balanced accuracy (BACC)']
  
  threshpoly_gbacc[, .(min(bacc), max(bacc))]
  
  threshpoly_gsenspec <- gstats[
    which(abs(sens-spec)<=min(abs(sens-spec)+0.01)),] %>%
    .[,perfstat := 'Sensitivity = Specificity']
  threshpoly_gsenspec[, .(min(abs(sens-spec)), max(abs(sens-spec)))]
  
  thresh_ghexdt <- rbindlist(list(
    threshpoly_gpredratio,
    threshpoly_gmisclass,
    threshpoly_gbacc,
    threshpoly_gsenspec
  ))
  
  thresh_gpolydt <-thresh_ghexdt[, .SD[chull(threshold_u10, threshold_o10)],
                                 by=perfstat]
  
  minnetrhesh <- min(threshrange_network)
  maxnethresh <- max(threshrange_network)
  boxnet <- data.table(xstart=c(minnetrhesh, minnetrhesh,
                                maxnethresh , maxnethresh ),
                       xend=c(minnetrhesh, maxnethresh ,
                              minnetrhesh, maxnethresh ),
                       ystart=c(minnetrhesh, minnetrhesh,
                                maxnethresh , maxnethresh ),
                       yend=c(maxnethresh , minnetrhesh,
                              maxnethresh , minnetrhesh))
  
  sensitivityplot_gperf <- ggplot(data=threshpoly_gpredratio,
                                  aes(x=threshold_u10, y=threshold_o10)) +
    geom_vline(xintercept=0.5, color='black') +
    geom_hline(yintercept=0.5, color='black') +
    geom_segment(data=boxnet, aes(x=xstart, y=ystart, xend=xend, yend=yend),
                 color='darkgrey') +
    stat_bin2d(aes(fill=perfstat),
               stat='identity', alpha=1/2) +
    stat_bin2d(data=threshpoly_gbacc, aes(fill=perfstat),
               stat='identity', alpha=1/2) +
    stat_bin2d(data=threshpoly_gmisclass, aes(fill=perfstat),
               stat='identity', alpha=1/2) +
    stat_bin2d(data=threshpoly_gsenspec, aes(fill=perfstat),
               stat='identity', alpha=1/2) +
    scale_fill_manual(
      name='Zone within 1% of optimum',
      values=c("#999999", "#E69F00", "#56B4E9", "#D55E00")) +
    # geom_point(data=minmax_gpoints, size=2) +
    # geom_text(data=minmax_gpoints, aes(label=perfstat),
    #   size=3, hjust=0.5, vjust=-0.8) +
    # geom_text(data=minmax_gpoints, aes(
    #   label=paste0(' x=', threshold_u10, ', y=', threshold_o10)),
    #   size=3, hjust=0.5, vjust=1.1) +
    scale_x_continuous(
      name=bquote('Probability threshold for gauges with MAF < 10'~m^3~s^-1)) +
    scale_y_continuous(
      name=expression(Probability~threshold~"for"~gauges~with~MAF >= 10~m^3~s^-1)) +
    coord_fixed(clip='off') +
    theme_classic() +
    theme(legend.position = c(0.20, 0.12),
          legend.background = element_blank())
  
  #--------------- Compute sensitivity for global network -------------------
  in_rivpred[, LENGTH_KM_NOLAKE := LENGTH_KM*(1-INLAKEPERC)]
  
  threshgrid_net <- expand.grid(threshrange_network, threshrange_network) %>%
    setnames(c('thresh_u10', 'thresh_o10'))
  
  netstats <- mapply(function(tu10, to10) {
    print(paste(tu10, to10))
    in_rivpred[(dis_m3_pyr>=mincutoff) & (INLAKEPERC < 1),
               list(
                 IRESperc = (
                   .SD[dis_m3_pyr<10,
                       sum(LENGTH_KM_NOLAKE*as.numeric(get(netcol)>=tu10))] +
                     .SD[dis_m3_pyr>=10,
                         sum(LENGTH_KM_NOLAKE*as.numeric(get(netcol)>=to10))]
                 )/
                   sum(LENGTH_KM_NOLAKE),
                 threshold_u10 = tu10,
                 threshold_o10 = to10
               )
    ]
  },
  tu10 = threshgrid_net$thresh_u10,
  to1 = threshgrid_net$thresh_o10,
  SIMPLIFY = F
  ) %>%
    do.call(rbind, .)
  
  predbounds <- netstats[
    (threshold_u10 %in% c(minnetrhesh, maxnethresh ))
    &
      (threshold_o10 %in%  c(minnetrhesh, maxnethresh )),
  ]
  
  sensitivityplot_netpred <- ggplot(
    netstats, aes(x=threshold_u10, y=threshold_o10)) +
    geom_vline(xintercept=0.5, color='black') +
    geom_hline(yintercept=0.5, color='black') +
    geom_bin2d(aes(fill=IRESperc), stat='identity', alpha=0.8) +
    scale_fill_distiller(name='Predicted % of IRES',
                         palette='Spectral') +
    scale_x_continuous(
      name=bquote('Probability threshold for gauges with MAF < 10'~m^3~s^-1)) +
    scale_y_continuous(
      name=expression(Probability~threshold~"for"~gauges~with~MAF >= 10~m^3~s^-1)) +
    coord_fixed(clip='off', expand=c(0,0)) +
    theme_classic()
  
  return(list(gperf = sensitivityplot_gperf,
              netpred = sensitivityplot_netpred,
              predbounds = predbounds
  ))
}


#------ tabulate_globalsummary -----
#' Tabulate global summary
#'
#' Create a summary table of the predicted prevalence of IRES across continuous size classes
#' and categorical classes.
#' 
#' @param outp_riveratlaspred (character) full path to .csv table containing flow
#'  intermittence predictions for global river network
#' @param inp_riveratlas (character) full path to .csv. table containing attributes
#'  for global river network
#' @param inp_riveratlas_legends (character) full path to .xlsx table containing 
#' metadata on river atlas attributes (meaning of category acronyms and codes e.g., for land cover)
#' @param idvars (character) name of the variable (usually categorical) whose 
#' categories to use for producing summary statistics of flow prevalence (e.g., climate or country)
#' @param interthresh (numeric) flow intermittence probability threshold above which
#'  to classify records as non-perennial.
#' @param castvar (character) name of the column that will be used to define bins/table columns (e.g., discharge) —
#' the equivalent of \code{binvar} in \link{bin_dt}.
#' @param castvar_num (boolean) whether casting column (castvar) is a numeric (used for sorting properly).
#' @param weightvar (character) variable to weigh proportion by (e.g. river reach length).
#' @param valuevar (character) variable to summarize (e.g. intermittency class).
#' @param valuevarsub (character) value of \code{valuevar} to summarize (e.g. '1' for intermittent rivers)
#' @param binfunc (character) binning approach. One of 'manual', 'equal_length', 'equal_freq'. See \link{bin_dt}.
#' @param binarg  (numeric) binning argument, depends on binning approach. See \link{bin_dt}.
#' @param bintrans (character or numeric) transformation of \code{binvar}, default is NULL.
#' @param na.rm (logical) whether to include NAs.
#' @param tidy (boolean) whether to keep in tidy long form (TRUE) or to cast into
#' formatted shape (FALSE; i.e., idvars as rows and castvar as columns + a world row)
#' @param nolake (boolean) whether to include sections of river reaches that intersect lakes (FALSE) or not (TRUE)
#' @param mincutoff (numeric) smallest river reaches to include in computation, in terms of long-term mean annual flow (m3/s; e.g., 0.1)
#' 
#' @details this function is used to produce Table 1 in Main Text, Extended Data Figure 1c, and 
#' Supplementary Data tables from Messager et al. (2021). 
#' 
#' bintrans can either be 'log' (for natural log) or a numeric exponent to transform
#' according to x^bintrans.
#' 
#' @return data.table
#'
#' @export
tabulate_globalsummary <- function(outp_riveratlaspred,
                                   inp_riveratlas,
                                   inp_riveratlas_legends,
                                   idvars,
                                   interthresh,
                                   castvar, castvar_num=TRUE,
                                   weightvar,
                                   valuevar, valuevarsub,
                                   binfunc=NULL, binarg=NULL, bintrans=NULL,
                                   na.rm=T, tidy=FALSE,
                                   nolake = TRUE,
                                   mincutoff = 0) {
  
  
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
    riveratlas <- riveratlas[, LENGTH_KM_NOLAKE := LENGTH_KM*(1-INLAKEPERC)]
    weightvar = 'LENGTH_KM_NOLAKE'
  }
  
  riveratlas <- riveratlas[dis_m3_pyr >= mincutoff,]
  
  riveratlas[, sum(LENGTH_KM_NOLAKE)]
  riveratlas[, mean(LENGTH_KM_NOLAKE)]
  riveratlas[, .N]
  
  
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
  
  format_tidyperc <- function(interthresh,
                              riveratlas, castvar, newvar, binfunc,
                              binarg, bintrans, na.rm, cast, tidy) {
    newvar <- paste0(valuevar, 'cat')
    riveratlas[, (newvar) := fifelse(get(valuevar)>=interthresh, valuevarsub, 0)]
    
    #Bin castvar if needed
    if (!is.null(binfunc) & !is.null(binarg)) {
      riveratlas_format <- bin_dt(in_dt = riveratlas, binvar = castvar, valuevar = newvar,
                                  binfunc = binfunc, binarg = binarg, bintrans = bintrans,
                                  na.rm = na.rm)
      castvar = 'bin_lformat'
      castvar_num = FALSE
    }
    
    #Compute overall number or weight for all combinations of castvar, valuevar and idvars
    #e.g. total number or river length for each flow state category for each country and river order
    statall <- riveratlas_format[if (na.rm) {!is.na(eval(newvar))},
                                 if (!is.null(weightvar)) sum(get(weightvar)) else .N,
                                 by=c(eval(castvar), eval(newvar), eval(idvars))] %>%
      rbind(
        .[, list('World', sum(V1)), by=c(eval(castvar), eval(newvar))] %>%
          setnames(c('V1', 'V2'), c(idvars, 'V1'))
      )
    
    
    #Compute for each cast and id var, the percentage of the newvar that is of the category of interest
    #(e.g. for each country and river order, percentage of river length that is intermittent)
    tidyperc <- statall[, 100*.SD[get(newvar)==valuevarsub, sum(V1)]/sum(V1),
                        by=c(eval(castvar), eval(idvars))]
    
    #If cast variable is a numeric (e.g. Strahler Order), sort table in increasing order
    if (castvar_num) {
      tidyperc[, eval(castvar) := as.numeric(as.character(get(eval(castvar))))]
    }
    setorderv(tidyperc, cols=c(castvar, idvars))
    
    
    #Compute totals for each newvar and idvar summed across castvar
    #e.g. (total length of rivers for each country, and of each flow regime)
    tidytotal <- statall[, list('Total intermittency (%)',
                                100*.SD[get(newvar)==valuevarsub, sum(V1)]/sum(V1)),
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
  
  tidyperc_list <- lapply(interthresh, format_tidyperc,
                          riveratlas, castvar, newvar, binfunc,
                          binarg, bintrans, na.rm, cast, tidy)
  
  names(tidyperc_list) <- interthresh
  
  return(tidyperc_list)
}

#------ extend_globalsummary ---------
#' Extend global summary
#'
#' Add statistics for extrapolated discharge size class to table of summary statistics
#' on global prevalence
#' 
#' @param in_IRESextra estimate of prevalence of flow intermittence in global river 
#' network for river reaches with 0.01 < MAF 0.099 
#' @param in_globaltable summary data.table of the predicted prevalence of IRES 
#' across discharge size classes and climate zones
#' @param inp_riveratlas_legends (character) full path to .xlsx table containing 
#' metadata on river atlas attributes (meaning of category acronyms and codes e.g., for land cover)  
#' 
#' @details this function is used to produce Table 1 in Main Text, Extended Data Figure 1c, and 
#' Supplementary Data tables from Messager et al. (2021). 
#' 
#' @return data.table
#'
#' @export
extend_globalsummary_clz <- function(in_IRESextra, in_globaltable, 
                                     inp_riveratlas_legends) {
  
  in_IRESextra <- in_IRESextra$preds
  
  worldpred_belowcutoff <- in_IRESextra[
    , list(`0.01-0.099` = 100*sum(percinter_gam_belowcutoff*cumL_predextra, na.rm=T)/sum(cumL_predextra, na.rm=T),
           clz_cl_cmj = 'World')]
  
  clzpred_belowcutoff <- in_IRESextra[
    , list(`0.01-0.099` = 100*sum(percinter_gam_belowcutoff*cumL_predextra, na.rm=T)/sum(cumL_predextra, na.rm=T)),
    by=clz_cl_cmj] %>%
    rbind(worldpred_belowcutoff) %>%
    setorder(clz_cl_cmj)
  
  worldpred <- in_IRESextra[
    , list(`Total intermittence extra (%)` = 100*sum(percinter_all_GAM*cumL_pred, na.rm=T)/sum(cumL_pred, na.rm=T),
           `Total stream length extra` = sum(cumL_pred)/1000,
           clz_cl_cmj = 'World')]
  
  clzpred <- in_IRESextra[
    , list(`Total intermittence extra (%)` = 100*sum(percinter_all_GAM*cumL_pred, na.rm=T)/sum(cumL_pred, na.rm=T),
           `Total stream length extra` = sum(cumL_pred, na.rm=T)/1000
    ),
    by=clz_cl_cmj] %>%
    rbind(worldpred) %>%
    setorder(clz_cl_cmj)
  
  clznames <- readxl::read_xlsx(inp_riveratlas_legends, sheet='clz_cl') %>%
    setDT
  
  clzpred_bind <- cbind(clzpred_belowcutoff, clzpred[, -c('clz_cl_cmj'), with=F]) %>%
    .[, clz_cl_cmj := as.numeric(as.character(clz_cl_cmj))] %>%
    merge(clznames, by.x = 'clz_cl_cmj', by.y='GEnZ_ID', all.x=T) %>%
    .[is.na(clz_cl_cmj), GEnZ_Name := 'World']
  
  tableextend <- merge(in_globaltable[[0.5]],
                       clzpred_bind[, -c('GEnZ_Code', 'clz_cl_cmj'), with=F],
                       by.x='clz_cl_cmj', by.y = 'GEnZ_Name') %>%
    setorder(-'Total stream length extra')
  
  
  numcols <-  names(tableextend)[unlist(tableextend[, lapply(.SD, is.numeric)])]
  
  tableextend[, (numcols) := lapply(.SD, function(x) {
    fifelse(is.na(x), '-', as.character(round(x)))}),
    .SDcols=numcols] %>%
    .[, `:=`(`Total intermittence (without extrapolation) - %` = paste0(
      `Total intermittence extra (%)`, ' (', `Total intermittency (%)`, ')'),
      `Total intermittence extra (%)` = NULL,
      `Total intermittency (%)` = NULL,
      `Total river length (without extrapolation) - 10^3 km` = paste0(
        `Total stream length extra`, ' (', `Total stream length (10^3 km)`, ')'),
      `Total stream length (10^3 km)` = NULL,
      `Total stream length extra` = NULL
    )]
  
  return(tableextend)
}

#------ compare_fr --------------------------------------
#' Compare model estimates for France
#'
#' Compare model estimates generated for mainland France by this study to model
#' estimates generated by Snelder et al. (2013) by producing comparative histograms
#' and statistics.
#' 
#' @param inp_frdir (character) full path to directory containing formatted
#'  network and basins for France. 
#' @param in_rivpred data.table of model predictions for global river network. Here, output from \link{netpredformat}.
#' @param predcol (character) name of the column for the predicted probability 
#' of flow intermittence generated by this study (Messager et al. 2021).
#' @param binarg  (numeric) limits of the mean annual flow bins across which to compare estimates (and for plotting histogram). See \link{bin_dt}.
#' @param mincutoff (numeric) smallest river reaches to include in comparison, in terms of long-term mean annual flow (m3/s; e.g., 0.1)
#' 
#' @details this function is used to produce Extended Data Figure 5 in Main Text from Messager et al. (2021). 
#' 
#' @return list with two elements:
#' *"plot": histogram comparing estimated prevalence of flow intermittence and 
#' total river length across discharge size classes.
#' *"data": summary comparison statistics.
#' @export
compare_fr <- function(inp_frdir, in_rivpred, predcol, binarg,
                       mincutoff) {
  in_netpath <- file.path(inp_frdir, 'network')
  in_baspath <- file.path(inp_frdir, 'hydrobasins12')
  valuevarsub <- "1"
  
  net_all <- st_read(dsn = dirname(in_netpath),
                     layer = basename(in_netpath))
  net <- net_all[((net_all$rhtvs2_all_phi_qclass_MODULE >= mincutoff) |
                    (net_all$rhtvs2_all_phi_qclass_SURF_BV >= 10)),]
  
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
                              datname = 'France') %>%
    .[, binsumlength := binsumlength/1000]
  
  tidyperc_riv  <- formathistab(in_dt = rivpredsub,
                                castvar = 'dis_m3_pyr',
                                valuevar = predcol,
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
  print(datmerge[bin > 1, weighted.mean(perc, binsumlength), by=dat])
  
  return(
    list(
      plot = ggcompare(datmerge, binarg) +
        geom_rect(aes(xmin=-Inf, xmax=1.5, ymin=0, ymax=Inf),
                  fill='grey', alpha=0.09) +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank()),
      data = datmerge
    )
  )
}

#------ compare_us ----------------
#' Compare model estimates for the U.S.
#'
#' Compare model estimates generated for the contiguous U.S. by this study to the
#' prevalence of intermittence as depicted in the National Hydrography Dataset
#' Plus (NHDPlus) medium-resolution (100k) by producing comparative histograms
#' and statistics.
#' 
#' @param inp_usresdir (character) full path to directory containing formatted
#'  network and basins data for the U.S. 
#'  @param inp_usdatdir (character) full path to directory containing raw data
#'  for the U.S.
#' @param in_rivpred data.table of model predictions for global river network. Here, output from \link{netpredformat}.
#' @param predcol (character) name of the column for the predicted probability 
#' of flow intermittence generated by this study (Messager et al. 2021).
#' @param binarg  (numeric) limits of the mean annual flow bins across which to compare estimates (and for plotting histogram). See \link{bin_dt}.
#' @param mincutoff (numeric) smallest river reaches to include in comparison, in terms of long-term mean annual flow (m3/s; e.g., 0.1)
#' 
#' @details this function is used to produce Extended Data Figure 5 in Main Text from Messager et al. (2021). 
#' 
#' @return list with two elements:
#' *"plot": histogram comparing estimated prevalence of flow intermittence and 
#' total river length across discharge size classes.
#' *"data_all": summary comparison statistics assuming that line segments in the 
#' NHDPlus classified as "Artificial path" are all perennial.
#' *"data_noartificial": summary comparison statistics excluding line segments 
#' classified as "Artificial path" are all perennial.
#' @export
compare_us <- function(inp_usresdir, inp_usdatdir, in_rivpred, predcol,
                       binarg, mincutoff) {
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
  rivpredsubmr <- rivpredbas[HUC8 %in% unique(netmr$HUC8),]
  
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
  
  
  netmr_o01 <- netmr[QE_MAm3 >= mincutoff,]
  
  print(paste0(
    'Range of prevalence of intermittency in NHDplus medium resolution >= 1 m3/s,',
    'depending on whether unclassified reaches are counted as fully intermittent or perennial: ',
    round(100*netmr_o01[FCODE %in% c(46003, 46007), sum(LENGTHKM)]/
            netmr_o01[FCODE %in% keepfcodes, sum(LENGTHKM)]),
    '-',
    round(100*netmr_o01[FCODE %in% c(46003, 46007, 46000, 55800), sum(LENGTHKM)]/
            netmr_o01[FCODE %in% keepfcodes, sum(LENGTHKM)]),
    '%'))
  
  print(paste0(
    'Total estimated prevalence in HydroSHEDS: ',
    rivpredsubmr[dis_m3_pyr >= mincutoff & get(predcol)==1, sum(LENGTH_KM)]/
      rivpredsubmr[dis_m3_pyr >= mincutoff, sum(LENGTH_KM)]
  ))
  
  
  #Check the percentage of each type of line in the NHD by HUC8 and Drainage area
  netmr_fsub <- copy(netmr[FCODE %in% keepfcodes,])
  netmr_fsub[, HUC8sum := sum(LENGTHKM), by=HUC8]
  FCode_HUC8 <- netmr_fsub[FCODE %in% keepfcodes,
                           as.integer(100*.SD[,sum(LENGTHKM)]/max(HUC8sum)),
                           by=.(HUC8, FCODE)]
  
  netmrbinned <- bin_dt(in_dt=netmr,
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
  # tidyperc_nhdHUC <- netmr_fsub[!is.na(intermittent),][
  #   ,
  #   formathistab(in_dt = .SD,
  #                castvar = "QE_MAm3",
  #                valuevar = "intermittent",
  #                valuevarsub = as.numeric(as.character(valuevarsub)),
  #                weightvar = "LENGTHKM",
  #                binfunc = 'manual',
  #                binarg =  binarg,
  #                binlabels = binlabels,
  #                datname = 'U.S. NHD'),
  #   by=HUC8]
  #
  # tidyperc_rivHUC  <- rivpredsubmr[
  #   , formathistab(in_dt = .SD,
  #                  castvar = 'dis_m3_pyr',
  #                  valuevar = 'get(predcol)',
  #                  valuevarsub = valuevarsub,
  #                  weightvar = 'LENGTH_KM',
  #                  binfunc = 'manual',
  #                  binarg =  binarg,
  #                  binlabels = binlabels,
  #                  datname = 'Global predictions'),
  #   by=HUC8]
  #
  # hucmajcl <- rivpredsubmr[, .N,
  #                          by=.(clz_cl_cmj, HUC8)][
  #                            , .SD[which.max(N), clz_cl_cmj], by=HUC8] %>%
  #   setnames('V1', 'majcl')
  #
  # setnames(tidyper_nhdHUC, paste0(names(tidyper_nhdHUC), 'NHD'))
  # tidyHUCmerge <- merge(tidyperc_rivHUC, tidyperc_nhdHUC,
  #                       by.x=c('HUC8','bin'), by.y=c('HUC8NHD','binNHD'),
  #                       all.x = T, all.y=T) %>%
  # ggplot(tidyHUCmerge, aes(x=perc, y=percNHD, color=factor(majcl))) +
  #   geom_point() +
  #   facet_wrap(~binformat)
  #
  # tidyperc_hucbind <- rbind(tidyperc_nhdHUC, tidyperc_rivHUC) %>%
  #   merge(hucmajcl, by='HUC8', all.x=T)
  #
  # clhucbox <- ggplot(tidyperc_hucbind, aes(x=factor(majcl), y=perc, fill=dat)) +
  #   geom_boxplot() +
  #   facet_wrap(~binformat)
  #
  # Check what difference is due to differences in drainage density
  # compare_dens <- dcast(tidyperc_hucbind, HUC8+bin~dat,
  #                       value.var=c('binsumdens', 'perc', 'binsumlength')
  # )
  # compcols <- names(compare_dens[
  #   , which(as.vector(unlist(lapply(compare_dens, is.numeric)))), with=F])
  # compare_dens[, (compcols) := lapply(.SD, function(x)
  #   {fifelse(is.na(x), 0, x)}), .SDcols=compcols] %>%
  #   .[,`:=`(densdiff = (`binsumdens_Global predictions` - `binsumdens_U.S. NHD`)/
  #             ((`binsumdens_Global predictions` + `binsumdens_U.S. NHD`)/2),
  #           percdiff = (`perc_Global predictions` - `perc_U.S. NHD`)/
  #             ((`perc_Global predictions` + `perc_U.S. NHD`)/2))]
  #
  # ggplot(compare_dens, aes(x=percdiff, y=densdiff)) +
  #   geom_point() +
  #   geom_smooth() +
  #   facet_wrap(~bin, scales='free')
  
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
                                 datname = 'U.S.') %>%
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
                                          datname = 'U.S.') %>%
    .[, binsumlength := binsumlength]
  
  #Compute bin statistics for river networks
  tidyperc_riv  <- formathistab(in_dt = rivpredsubmr,
                                castvar = 'dis_m3_pyr',
                                valuevar = predcol,
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
    list(
      plot = ggcompare(datmerge_all, binarg=binarg[]) +
        geom_rect(aes(xmin=-Inf, xmax=1.5, ymin=0, ymax=Inf),
                  fill='grey', alpha=0.09) +
        scale_y_continuous(breaks=c(0,25,75,100)) +
        coord_cartesian(ylim=c(0,120), expand=c(0,0)) +
        theme(axis.title.y = element_blank(),
              legend.title = element_blank(),
              legend.position = c(0.8, 0.2)),
      data_all = datmerge_all,
      data_noartificial = datmerge_noartificial
    )
  )
}


#------ compare_au ----------------
#' Compare model estimates for Australia
#'
#' Compare model estimates generated for Australia by this study to model
#' estimates in the Australian National Geofabric by producing comparative histograms
#' and statistics.
#' 
#' @param inp_resdir (character) full path to directory containing formatted
#'  network and basins data for Australia. 
#' @param in_rivpred data.table of model predictions for global river network. Here, output from \link{netpredformat}.
#' @param predcol (character) name of the column for the predicted probability 
#' of flow intermittence generated by this study (Messager et al. 2021).
#' @param binarg  (numeric) limits of the mean annual flow bins across which to compare estimates (and for plotting histogram). See \link{bin_dt}.
#' 
#' @details this function is used to produce Extended Data Figure 5 in Main Text from Messager et al. (2021). 
#' Only rivers and stream segments with a drainage area >= 10 km2 are included in the comparison.
#' 
#' @return list with two elements:
#' *"plot": histogram comparing estimated prevalence of flow intermittence and 
#' total river length across discharge size classes.
#' *"data": summary comparison statistics.
#' @export
compare_au <- function(inp_resdir, in_rivpred, predcol, binarg) {
  valuevarsub <- "1"
  
  #Get Australian network
  net <- fread(file.path(inp_resdir, 'netbas12_inters_australia.csv'))%>%
    .[, UpstrDArea := UpstrDArea/(10^6)]  %>%
    .[UpstrDArea >=10,]
  
  #Join HydroSHEDS basins with RiverATLAS network and subselect network to match NHD selection
  rivpredbas <- merge(in_rivpred, unique(net[, .(HYBAS_ID)]),
                      by.x="HYBAS_L12", by.y="HYBAS_ID", all.x=F) #To match full US
  rivpredsub <- rivpredbas[dis_m3_pyr > 0 & UPLAND_SKM >= 10,] #Remove segments with 0 discharge
  
  #Comparing by basin (adjusting prevalence of intermittence predicted by RF based
  # on length of rivers in Geofabric)
  # bastats_au <- net[,list(sumlen_IRES_au = .SD[Perennial == "Non Perennial", sum(LENGTH_GEO)],
  #                         sumlen_all_au = sum(LENGTH_GEO)),
  #                   by=HYBAS_ID]
  # bastats_riv <- rivpredsub[,list(sumlen_IRES_riv = sum(get(predcol)*LENGTH_KM),
  #                                 sumlen_all_riv = sum(LENGTH_KM)),
  #                           by=HYBAS_L12]
  #
  # bastats <- merge(bastats_au, bastats_riv, by.x='HYBAS_ID', by.y='HYBAS_L12')
  # bastats[, sum((sumlen_IRES_riv/sumlen_all_riv)*sumlen_all_au)/sum(sumlen_all_au)]
  # bastats[, sum((sumlen_IRES_au/sumlen_all_au)*sumlen_all_riv)/sum(sumlen_all_riv)]
  
  #Get binned prevalence of IRES and length
  binlabels <- label_manualbins(binarg=binarg,
                                minval=10)
  
  tidyperc_au <- formathistab(in_dt = net,
                              castvar = "UpstrDArea",
                              valuevar = "Perennial",
                              valuevarsub = "Non Perennial",
                              weightvar = "LENGTH_GEO",
                              binfunc = 'manual',
                              binarg =  binarg,
                              binlabels = binlabels,
                              datname = 'Australia')
  
  tidyperc_riv  <- formathistab(in_dt = rivpredsub,
                                castvar = 'UPLAND_SKM',
                                valuevar = predcol,
                                valuevarsub = "1",
                                weightvar = 'LENGTH_KM',
                                binfunc = 'manual',
                                binarg =  binarg,
                                binlabels = binlabels,
                                datname = 'Global')
  
  datmerge <- rbind(tidyperc_au, tidyperc_riv) %>%
    .[, dat:=factor(dat, levels=c('Global', 'Australia'))] %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]
  print('Percentage intermittence for Australia')
  print(datmerge[, weighted.mean(perc, binsumlength), by=dat])
  
  # dcast(datmerge, formula=binformat~dat, value.var = c('perc', 'binsumlength')) %>%
  #   .[, weighted.mean(perc_Global, binsumlength_Geofabric)]
  
  return(
    list(
      plot = ggcompare(datmerge, binarg, insetx = 0.85, insety = 0.9),
      data = datmerge
    )
  )
  
}

#------ qc_pnw ------------
#' Quality check model estimates in U.S. Pacific Northwest
#'
#' Analyze performance of model predictions of flow intermittence against PROSPER field
#' observations of flow state across the U.S. Pacific Northwest.
#' 
#' @param inp_pnwresdir (character) full path to directory containing formatted
#'  data for Pacific Northwest. 
#' @param in_rivpred data.table of model predictions for global river network. Here, output from \link{netpredformat}.
#' @param predcol (character) name of the column for the predicted probability 
#' of flow intermittence generated by this study (Messager et al. 2021).
#' @param interthresh (numeric) flow intermittence probability threshold above which
#'  to classify records as non-perennial.
#' @param mincutoff (numeric) smallest river reaches to include in comparison, in terms of long-term mean annual flow (m3/s; e.g., 0.1)

#' 
#' @details this function is used to format data used in producing Extended Data Figure 6.
#' A side effect of this function is the writing of a .gpkg of river reaches for
#' which there are flow observations in the U.S. Pacific Northwest with quality-checking
#' attributes.
#' 
#' @return list with two elements:
#' *"plot": scatterplots analyzing Intermittence Prediction Residuals (IPR) as a function of
#'  modeled mean annual flow, the position of the field observation along the RiverATLAS river reach,
#'  the number of field observations at the site, and the population density in the entire upstream area
#'  of the corresponding reach.
#' *"stats": summary comparison statistics.
#' @export
qc_pnw <- function(inp_pnwresdir, in_rivpred, predcol,
                   interthresh=0.5, mincutoff = 0) {
  in_refpts <- file.path(inp_pnwresdir, 'StreamflowPermObs_final')
  in_fulldat <- file.path(inp_pnwresdir, 'StreamflowPermObs_sub')
  
  valuevarsub <- "1"
  probcol <- gsub('cat', '', predcol)
  
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
  
  ################# for peer-review ############################################
  check_nobsinter <- fulldat[, list(nobs=.N,
                                    refinter=as.numeric(any(Category == "Non-perennial"))),
                             by=.(POINT_X, POINT_Y)]
  
  nobs_inter_ratio <- check_nobsinter[, .SD[refinter == 1, .N]/
                                        .SD[refinter == 0, .N],
                                      by = nobs] %>%
    setorder(nobs)
  
  summary(glm(refinter~nobs, data=check_nobsinter, family=binomial(link='logit')))
  
  
  ##############################################################################
  
  
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
                       in_rivpred[, c('HYRIV_ID', 'HYBAS_L12', 'dis_m3_pyr',
                                      probcol, predcol),
                                  with = F],
                       by='HYRIV_ID', all.y=F) %>%
    .[, refinter := fifelse(Category == 'Non-perennial', 1, 0)] #Assign ephemeral and intermittent categories to 1, perennial to 0 for PNW obs
  
  
  #Consider that if one location is non-perennial, that reach is not perennial
  refpts_reach <- refpts_join[
    , `:=`(refinter_reach = as.numeric(any(refinter == 1)),
           nobs_reach = sum(nobs)
    ), by=HYRIV_ID] %>%
    setorder(HYRIV_ID, -refinter) %>%
    .[!duplicated(HYRIV_ID),] %>%
    .[dis_m3_pyr >= mincutoff,] %>%
    .[, IPR := get(probcol) - refinter] #Compute prediction error
  
  
  #Compute BACC
  refpts_reach[, `:=`(truth = factor(as.character(refinter_reach),levels=c('0', '1')),
                      response = factor(get(predcol),levels=c('0','1'))
  )]
  
  pnw_measures <- refpts_reach[, list(
    bacc = mlr3measures::bacc(truth, response),
    auc = mlr3measures::auc(truth, get(probcol), positive='1'),
    ce = mlr3measures::ce(truth, response),
    sen = mlr3measures::sensitivity(truth, response, positive='1'),
    spe = mlr3measures::specificity(truth, response, positive='1'),
    pre = mlr3measures::precision(truth, response, positive='1'),
    nobs_total = sum(nobs),
    nreaches = .N,
    nreaches_perennial = .SD[truth=='0', .N],
    nreaches_nonperennial = .SD[truth=='1', .N]
  )]
  
  #Compute relative position of observation point on HydroSHEDS reach line
  refpts_reach[, RATIOLENGTH := min(c(1, (fromM*95)/LENGTH_KM)), by=OBJECTID]
  
  #Write data out to points
  predtowrite <- merge(
    st_read(dsn = inp_pnwresdir,
            layer = basename(in_refpts))[,c('OBJECTID', 'Shape')],
    refpts_reach[,-c('Shape','i.Shape'), with=F],
    by='OBJECTID', all.x=F, all.y=T) %>%
    setnames(old=c("HYDROSHEDSdis", "HYDROSHEDSDA",
                   predcol, probcol),
             new=c("LINEdis", "LINEDA", "predcat", "predp")) %>%
    setnames(new=unlist(lapply(names(.), function(x) gsub('[.]','_', x))))
  
  st_write(obj=predtowrite[,c('OBJECTID', 'HYRIV_ID',
                              'IPR',
                              "LINEdis", "LINEDA",
                              "predcat", "predp",
                              'nobs', 'refinter_reach',
                              'pop_ct_usu')],
           dsn = dirname(inp_pnwresdir),
           layer=paste0('pnwobs_IPR_',
                        predcol,
                        format(Sys.time(), '%Y%m%d%H%M'), '.shp'),
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
  
  return(list(plot = pnw_qcplot,
              stats = pnw_measures)
  )
}

#------ qc_onde ------
#' Quality check model estimates in France
#'
#' Analyze performance of model predictions of flow intermittence against ONDE field
#' observations of flow state across mainland France.
#' 
#' @param inp_ondedatdir (character) full path to directory containing raw
#'  data for France.
#' @param inp_onderesdir (character) full path to directory containing formatted
#'  data for France.
#' @param inp_riveratlas (character) full path to .csv. table containing attributes
#'  for global river network
#' @param in_rivpred data.table of model predictions for global river network. Here, output from \link{netpredformat}.
#' @param predcol (character) name of the column for the predicted probability 
#' of flow intermittence generated by this study (Messager et al. 2021).
#' @param interthresh (numeric) flow intermittence probability threshold above which
#'  to classify records as non-perennial.
#' @param mincutoff (numeric) smallest river reaches to include in comparison, in terms of long-term mean annual flow (m3/s; e.g., 0.1)
#' 
#' @details this function is used to format data used in producing Extended Data Figure 6.
#' A side effect of this function is the writing of a .gpkg of river reaches for
#' which there are flow observations in France with quality-checking
#' attributes.
#' 
#' @return list with two elements:
#' *"plot": scatterplots analyzing Intermittence Prediction Residuals (IPR) as a function of
#'  modeled mean annual flow, the position of the field observation along the RiverATLAS river reach,
#'  the ratio of the drainage area of the watercourse at the location of the field 
#'  observation to that for the pourpoint of the corresponding RiverATLAS river reach,
#'  the percentage of field observations that are either no-flow or dry,
#'  the number of field observations at the site, and the population density in the entire upstream area
#'  of the corresponding reach.
#' *"stats": summary comparison statistics.
#' @export
qc_onde <- function(inp_ondedatdir, inp_onderesdir, inp_riveratlas,
                    in_rivpred, predcol, interthresh=0.5, mincutoff=0) {
  in_refpts <- file.path(inp_onderesdir, 'obs_finalwgs')
  in_fulldat <- file.path(inp_ondedatdir, 'onde_france_merge.csv')
  
  valuevarsub <- "1"
  probcol <- gsub('cat', '', predcol)
  
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
                       in_rivpred[, c('HYRIV_ID', 'HYBAS_L12', 'LENGTH_KM',
                                      probcol, predcol, 'dis_m3_pyr'),
                                  with = F],
                       by.x='HYRIVIDjoinedit', by.y='HYRIV_ID', all.y=F) %>%
    merge(riveratlas, by.x='HYRIVIDjoinedit', by.y='HYRIV_ID', all.y=F) %>%
    .[, IPR := get(probcol) - refinter] #Compute prediction error
  
  #Consider that if one location is non-perennial, that reach is not perennial
  refpts_reach <- refpts_join[
    , refinter_reach := as.numeric(any(refinter == 1)), by=HYRIVIDjoinedit] %>%
    setorder(HYRIVIDjoinedit, -refinter) %>%
    .[!duplicated(HYRIVIDjoinedit),] %>%
    .[dis_m3_pyr >= mincutoff,] %>%
    .[, IPR := get(probcol) - refinter] #Compute prediction error
  
  
  #Compute BACC
  refpts_reach[, `:=`(truth = factor(as.character(refinter_reach),levels=c('0', '1')),
                      response = factor(get(predcol),levels=c('0','1'))
  )]
  onde_measures <- refpts_reach[, list(
    bacc = mlr3measures::bacc(truth, response),
    auc = mlr3measures::auc(truth, get(probcol), positive='1'),
    ce = mlr3measures::ce(truth, response),
    sen = mlr3measures::sensitivity(truth, response, positive='1'),
    spe = mlr3measures::specificity(truth, response, positive='1'),
    pre = mlr3measures::precision(truth, response, positive='1'),
    nobs_total = sum(nobs),
    nreaches = .N,
    nreaches_perennial = .SD[truth=='0', .N],
    nreaches_nonperennial = .SD[truth=='1', .N]
  )]
  
  #Compute relative position of observation point on HydroSHEDS reach line
  refpts_join[, RATIOLENGTH := (fromM/1000)/LENGTH_KM]
  
  #Write data out to points
  predtowrite <- merge(st_read(dsn = inp_onderesdir,
                               layer = basename(in_refpts))[, "F_CdSiteHydro_"],
                       refpts_join,
                       by.x = "F_CdSiteHydro_", by.y = 'CdSiteHydro') %>%
    setnames(old=c("HYDROSHEDSdis", "HYDROSHEDSDA",
                   predcol, probcol),
             new=c("LINEdis", "LINEDA", "predcat", "predp"))
  
  write_sf(predtowrite, file.path(dirname(inp_onderesdir),
                                  paste0(
                                    'ondeobs_IPR_',
                                    predcol,
                                    format(Sys.time(), '%Y%m%d%H%M'), '.shp')
  ))
  
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
  
  return(list(plot=onde_qcplot,
              stats = onde_measures)
  )
}



#################################### END #######################################