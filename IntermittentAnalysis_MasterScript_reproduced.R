library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

in_filestructure <- readd(filestructure)
in_gaugep <- readd(gaugep)
in_gaugestats <- readd(gaugestats_format)
in_rftuned <- readd(rftuned)
in_predvars <- readd(predvars)

loadd(rfbm_classif)
in_rf = rfbm_classif$bm_classif$clone()$filter(learner_ids = "oversample.classif.ranger.tuned")
in_task = rfbm_classif$bm_tasks$task_classif
in_measure = rfbm_classif$measure_classif
featimpfilt = 0.01
insamp_nfolds =  2
insamp_nevals = 1
outsamp_nrep = 2
outsamp_nfolds =  10
imp = TRUE
perf = TRUE
pvalue = TRUE
pvalue_permutn = 10


rfresamp <-  in_rf$resample_result(
  uhash=unique(as.data.table(in_rf)$uhash))

bmcheck_classif <- readd(bm_checked_rfbm_classif.bm_classif_rfbm_classif.measure_classif)
bmcheck_regr <- readd(bm_checked_rfbm_regr.bm_regr_rfbm_regr.meassure_regr)
loadd(rfbm_classif)
loadd(rfbm_regr)
loadd(rfeval_featsel)



bmres <- rfbm_classif$bm_classif
check <- lapply(1:bmres$n_resample_results, function(i) {
  print(i)
  rsmp_preds <- bmres$resample_result(rsmp_i)$prediction()

  baccthresh <- ldply(seq(0,1,0.01), function(threshold_class) {
    print(threshold_class)
    cbind(
      rsmp_bacc(bmres, rsmp_i=i, threshold_class=threshold_class),
      threshold_misclass(i=threshold_class, in_preds=rsmp_preds)
    ) %>%
      .[, i:=NULL]
  }) %>%
    setDT
  return(baccthresh)
})

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



















################### KEEP TO ANALYZE SPATIAL RESAMPLING OUTPUT ###################
check <- nestedresamp_bmrout_classif$resample_result(2)
coordfolds <- cbind(nestedresamp_bmrout_classif$tasks$task[[1]]$coordinates(), check$resampling$instance)

library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world)+
  geom_sf(data=sf::st_as_sf(x =coordfolds, coords = c("X", "Y"),
                            crs = "+proj=longlat +datum=WGS84"),
          aes(color=as.factor(fold)))



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
