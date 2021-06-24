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

loadd(rfbm_classif)
loadd(rfbm_regr)
loadd(rfeval_featsel)
loadd(gaugep)
loadd(rftuned)
gaugepred <- readd(rfpreds)
loadd(gaugestats_format)

in_filestructure <- filestructure

#Import global predictions
rivpred <- fread(in_filestructure['out_riveratlas']) %>%
  .[rivernetwork[, c('HYRIV_ID', 'HYBAS_L12', 'LENGTH_KM', 'dis_m3_pyr', 'UPLAND_SKM'),
                 with=F], on='HYRIV_ID']
rivpredsub <- merge(rivpred, bas, by.x="HYBAS_L12", by.y="HYBAS_ID", all.x=F) %>%
  .[, UPLAND_SKM := round(UPLAND_SKM)]


in_netpath <- filestructure[['in_netfr']]
in_baspath <- filestructure[['in_basfr']]
valuevarsub <- "1"
castvar <- "rhtvs2_all_phi_qclass_SURF_BV"
valuevar <- "INT_RF_txt_V1"
weightvar <- "rhtvs2_all_phi_qclass_LONG_"
binarg <- c(10,20,50,100,200,500,1000,2000,5000,10000,
            50000,100000,150000)

label_manualbins <- function(binarg, minval) {
  minlabel <- paste(minval, binarg[1], sep="-")
  otherlabels <- mapply(function(x, y) {paste(x, y-1, sep="-")},
                        binarg[1:(length(binarg)-1)], binarg[2:length(binarg)])
  return(c(minlabel, otherlabels))
}

net <- st_read(dsn = dirname(in_netpath),
               layer = basename(in_netpath))

bas <- st_read(dsn = dirname(in_baspath),
               layer = basename(in_baspath)) %>%
  .[, 'HYBAS_ID', with=F]

#"rhtvs2_all_phi_qclass_MODULE" average discharge
#"rhtvs2_all_phi_qclass_SURF_BV" drainage area
#"rhtvs2_all_phi_qclass_LONG_" reach length
#"INT_RF_txt_V1"

binlabels <- label_manualbins(binarg=binarg, minval=min(net$rhtvs2_all_phi_qclass_SURF_BV))

netbin <- bin_dt(in_dt = as.data.table(net),
                 binvar = castvar,
                 valuevar = valuevar,
                 binfunc = "manual",
                 binarg = binarg)

netstat <- as.data.table(netbin)[,sum(get(weightvar)),
                              by=c(eval(valuevar), "bin")]
tidyperc <- netstat[, list(perc = 100*.SD[get(valuevar)==valuevarsub,
                                          sum(V1)]/sum(V1),
                           binsumlength = sum(V1)/1000),
                    by=c("bin")] %>%
  .[, `:=`(dat='Snelder et al. (2013) predictions',
           binformat = binlabels[bin])]


#Format RiverATLAS
castvar <- 'UPLAND_SKM'
valuevar <- 'predbasic800cat'
weightvar <- 'LENGTH_KM'
rivbin <- bin_dt(in_dt = as.data.table(rivpredsub),
                 binvar = castvar,
                 valuevar = valuevar,
                 binfunc = "manual",
                 binarg = binarg)

netstat <- as.data.table(rivbin)[,sum(get(weightvar)),
                                 by=c(eval(valuevar), "bin")]
tidyperc_riv <- netstat[, list(perc = 100*.SD[get(valuevar)==valuevarsub,
                                              sum(V1)]/sum(V1),
                               binsumlength = sum(V1)),
                        by=c('bin')] %>%
  setorder(bin) %>%
  .[, `:=`(dat='Global predictions',
           binformat = binlabels[bin])]

datmerge <- rbind(tidyperc, tidyperc_riv) %>%
  setorder(bin) %>%
  .[, binformat := factor(binformat, levels=unique(binformat))]

ggplot(datmerge, aes(x=binformat, y=perc, fill=dat)) +
  geom_bar(stat='identity', position='dodge') +
  #geom_smooth() +
  #scale_x_log10() +
  coord_cartesian(ylim=c(0,35), expand=FALSE, clip="off") +
  theme_classic()

ggplot(datmerge, aes(x=binformat, y=binsumlength, fill=dat)) +
  geom_bar(stat='identity', position='dodge') +
  #geom_smooth() +
  #scale_x_log10() +
  coord_cartesian(expand=FALSE, clip="off") +
  theme_classic()





######################################################
ggplot(rivpredsub, aes(x=dis_m3_pyr, y=predbasic800)) +
  geom_point(alpha=1/5) +
  geom_smooth(aes(weight=LENGTH_KM)) +
  geom_quantile(aes(weight=LENGTH_KM), method = "rqss") +
  scale_x_sqrt() +
  theme_classic()

ggplot(rivpredsub, aes(x=UPLAND_SKM, y=predbasic800)) +
  geom_point(alpha=1/5) +
  geom_smooth(aes(weight=LENGTH_KM)) +
  scale_x_sqrt() +
  theme_classic()




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
