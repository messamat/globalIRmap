library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

in_gaugestats <- readd(gaugestats_format)
in_predvars <- readd(predvars)
in_rftuned <- readd(rftuned)
plot(readd(vimp_plot))

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