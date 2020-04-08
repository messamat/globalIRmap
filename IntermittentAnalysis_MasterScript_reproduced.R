


#Import GRDC stations and only keep those < 200 m from a HydroSHEDS reach
GRDCpjoin <- st_read(dsn=outgdb, layer='GRDCstations_riverjoin') %>%
  .[.$station_river_distance<200,]

#Get data paths of daily records for GRDC stations
fileNames <- file.path(datdir, 'GRDCdat_day', paste(GRDCpjoin$GRDC_NO,  ".txt", sep=""))
which(do.call(rbind, lapply(fileNames, file.exists))) #Check whether any file does not exist


#### -------------------------- Define functions ------------------------------------------

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
rf_formula <- as.formula(paste0('intermittent~', 
                                paste(rfpredcols, collapse="+"), 
                                collapse=""))

task_inter <- mlr3::TaskClassif$new(id ='inter_basic',
                                    backend = gaugestats_join[!is.na(cly_pc_cav), c('intermittent', rfpredcols),with=F],
                                    target = "intermittent")

#Set default ranger learner with explicit parameter set
rangerlrn <- mlr3::lrn('classif.ranger', id = 'ranger', 
                       num.trees = 500,
                       #mtry = , #Default is square root of number of variables
                       min.node.size = 1, #Default is 1 for classification, 5 for regression and 10 for probability
                       replace = FALSE, #in ranger package, default replace is True but Boulesteix et al. 2012 show that less biased
                       sample.fraction = 0.632, #Default for sampling without replacement
                       split.select.weights = ,
                       always.split.variables= ,
                       respect.unordered.factors = 'ignore',
                       importance='permutation',
                       write.forest=TRUE, 
                       scale.permutation.importance = FALSE,
                       num.threads = bigstatsr::nb_cores(),
                       save.memory = FALSE,
                       verbose = TRUE,
                       splitrule = "gini",
                       num.random.splits = 1L,
                       keep.inbag = FALSE,
                       predict_type = 'prob',
                       max.depth = NULL #Unlimited depth default
                       )

mlr3::lrn('classif.ranger')$param_set

mlr_measures
measure = lapply(c("classif.sensitivity", 'classif.specificity'), msr)
basic800pred$score(measure)


tune_ps = ParamSet$new(list(
  ParamDbl$new("cp", lower = 0.001, upper = 0.1),
  ParamInt$new("minsplit", lower = 1, upper = 10)
))
tune_ps



#### -------------------------- Diagnose initial model -------------------------------------------------
# ---- Variable importance ----
varimp_basic <- data.table(varcode=names(mod_basic800$variable.importance), 
                           importance=mod_basic800$variable.importance)[
                             rfpredcols_dt, on='varcode' ]%>%
  .[, varname := factor(varname, levels=.[,varname[order(-importance)]])] %>%
  setorder(-importance)

ggplot(varimp_basic[1:30,],aes(x=varname, y=importance)) + 
  geom_bar(stat = 'identity') +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=8))

# ---- Check misclassification rate with different threshold probabilities ----
gaugestats_join[!is.na(cly_pc_cav), intermittent_predprob := mod_basic800$predictions[,'1']]

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

threshold_misclass[i==0.5,]
gaugestats_join[, intermittent_predcat := ifelse(intermittent_predprob>=0.3, 1, 0)]

# ---- Partial dependence plots ----
pdtiles_grid(mod = mod_basic800, dt = gaugestats_join, colnums = 5)

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
             predbasic800 := predict(mod_basic800, type='response', data=.SD)$predictions[,'1'], 
             .SDcols = c("HYRIV_ID", rfpredcols)]
  toc()
}

riveratlas[, predbasic800cat := ifelse(predbasic800>=0.3, 1, 0)]
fwrite(riveratlas[, c('HYRIV_ID', 'predbasic800cat'), with=F], file.path(resdir, 'RiverATLAS_predbasic800.csv'))










############# TO DO ##############
## Map uncertainty in predictions for each gauge — relate to length of record and environmental characteristics
## Compare gauge environmental characteristics compared to full river network, looking at confusion matrix results
## Understand variable associations with gauges

######## Model improvements

#Bootstrap without replacement to avoid bias
#Hyper-parameter tuning
#Implement conditional inference forest
#Use CAST package for variable selection and Leave one location out CV
#Use mlr3? (read mlr3 book) — test spatial and aspatial CV









##################################################################################################################################################
############################ MORE ALREADY SPED UP/TRANSCRIBED ####################################################################################
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
############################## MORE CHARLOTTE ####################################################################################################
### Leave-one-out cross-validation (run overnight)
jack_results <- jackRF(MyYVar = class, MyXPreds = predictors)
Sys.time()
write.table(jack_results, "jack_results.csv", sep = ",", row.name = FALSE)

setwd("D:/Charlotte/Current_datasets")
jack_results <- read.csv("jack_results.csv")

#Finding agreement per climate zone
rounded <- as.numeric(ifelse(jack_results[,3] < 0.40, 0, 1))
jack_rounded <- cbind(jack_results[,1:2], rounded, jack_results[,4])
colnames(jack_rounded) <- colnames(jack_results)

zonal_agreement <- integer(18)
zonal_disagreement <- integer(18)
for(i in 1:length(jack_results$observed == jack_results$predicted)){
  if(jack_rounded$observed[i] == jack_rounded$predicted[i]){
    zonal_agreement[jack_rounded$clz_cl_cmj[i]] <- as.numeric(zonal_agreement[jack_rounded$clz_cl_cmj[i]]) +1
  }else{
    zonal_disagreement[jack_rounded$clz_cl_cmj[i]] <- as.numeric(zonal_disagreement[jack_rounded$clz_cl_cmj[i]]) +1
  }
}

indices <- vector()
for(i in 1:18){
  if(zonal_agreement[i] ==0 & zonal_disagreement[i]==0){
    indices[i] <- FALSE
  }else{
    indices[i] <- TRUE
  }
}

zonal <- cbind(c(1:18), zonal_agreement, zonal_disagreement)
zonal <- zonal[indices,]
percent_agreement <- zonal[,2]/(zonal[,2] + zonal[,3])
zonal <- cbind(zonal, percent_agreement)

numPerZone <- zonal[,2] + zonal[,3]
names(numPerZone) <- zonal[,1]

names(percent_agreement) <- zonal[,1]

barplot(numPerZone,  
        xlab="Climate zone",
        ylab = "Number of GRDC stations",
        cex.lab = 1.5
)

barplot(percent_agreement,  
        xlab="Climate zone",
        ylab = "Percent correctly classified",
        cex.lab = 1.5
)

### Funding auc and plotting 
rounded <- as.numeric(ifelse(results < 0.40, 0, 1))
plot(roc(rounded, as.factor(jack_results[,1])))

DATA <- jack_results[,1:3]
par(mar= c(4, 2, 3, 1)+0.1)
auc.roc.plot(DATA= DATA, find.auc= TRUE, threshold = 40, which.model=1, xlab = "False positive rate", ylab = "True positive rate",
             main = "ROC curve")
auc <- auc(DATA = DATA, which.model = 1)


### Finding Cohen's Kappa and percent correctly classified
setwd("D:/Charlotte/Current_datasets")
jack_results <- read.csv("jack_results.csv")
rounded_jack <- round(jack_results[,3])
kappa_data <- cbind(rounded_jack, jack_results[,2])
kappa <- kappa2(kappa_data)
agreement <- agree(kappa_data, tolerance = 0)

#Find threshold that maximizes AUC and Cohen's Kappa
threshold_matrix <- matrix(ncol = 3, nrow = 1001)

count <- 1
pb <- progress_bar$new(total = length(seq(0.0, 1, by = 0.001)))

for(j in seq(0.0, 1, by = 0.001)){
  results <- jack_results[,3]
  rounded <- as.numeric(ifelse(results < j, 0, 1))
  kappa_data <- cbind(rounded, jack_results[,2])
  kappa <- kappa2(kappa_data)
  agreement <- agree(kappa_data, tolerance = 0)
  threshold_matrix[count,] <- c(j, kappa[[5]], agreement[[5]])
  count<- count +1
  #Progress bar
  pb$tick()
  Sys.sleep(1 / length(seq(0.0, 1, by = 0.001)))
}

colnames(threshold_matrix) <- c("Threshold", "Cohen's kappa", "% Agreement")

max_kap_loc <- which.max(threshold_matrix[,2])
max_ag_loc <- which.max(threshold_matrix[,3])
max_kappa <- threshold_matrix[max_kap_loc, 1:2]
max_agreement <- threshold_matrix[max_ag_loc, c(1,3)]

kap_agreement <- threshold_matrix[max_kap_loc, c(1,3)]

plot(threshold_matrix[,1],threshold_matrix[,2],type="l",col="red", ylim = c(0, 0.9),
     ylab = "Accuracy measure", xlab = "Threshold", cex.lab = 1.5)
points(max_kappa[1], max_kappa[2])
points(max_agreement[1], max_agreement[2]/100)
lines(threshold_matrix[,1],threshold_matrix[,3]/100,lty = 2, col="blue")
par(mar= c(6, 6, 3, 1)+0.1)

legend("topleft", c("Kappa", "PCC"),
       col=c("red", "blue"),lty=1:2)

my.legend.size <-legend("topright",c("Series1","Series2","Series3"),plot = FALSE)



### Predictions 
geospatial_analyst <- function(continent){
  setwd('D:/Charlotte/Reach_R_data')
  newdata <- readRDS(paste(continent,  ".rds", sep=""))
  continent_reaches <- readRDS(paste(continent, "_reaches.rds", sep = ""))
  climate_zones <- continent_reaches$clz_cl_cmj
  continent_predictions <- predict(RF_model, newdata =newdata, type="prob", progress= "window")
  continent_preds <- as.data.frame(cbind(newdata[,1], ifelse(continent_predictions[,2]<0.4,0,1), continent_predictions[,2], climate_zones))
  colnames(continent_preds) <- c("Reach_ID", "Class", "Prob_of_Intermittence", "Climate zone")
  st_geometry(continent_preds) <- continent_reaches$geometry
  setwd("D:/Charlotte/GIS")
  #st_write(continent_preds, dsn = paste(continent,"_preds", sep =""), driver ='ESRI Shapefile')
  percent_intermittent <- sum(continent_preds$Class)/length(continent_preds$Class)
  return(percent_intermittent)
}

percent_Af <-geospatial_analyst("Africa")
percent_Asia <-geospatial_analyst("Asia")
percent_Aus <-geospatial_analyst("Australia")
percent_Eur <-geospatial_analyst("Europe")
percent_Na <-geospatial_analyst("NorthAmerica")
percent_North <- geospatial_analyst("Northern")
percent_Sa <- geospatial_analyst("SouthAmerica")
percent_Sib <-geospatial_analyst("Siberia")

probabilities <- c(percent_Af, percent_Asia, percent_Aus, percent_Eur, 
                   percent_Na, percent_North, percent_Sa, percent_Sib)*100
names(probabilities) <- c("Africa", "Asia", "Australia", "Europe", "N America",
                          "Northern", "S America", "Siberia")
setwd("D:/Charlotte/Current_datasets")
write.table(probabilities, "probabilities.csv", row.names = FALSE, sep = ",")

barplot(probabilities,  
        xlab="Region",
        ylab = "% of rivers predicted to be intermittent",
)
par(mar= c(10, 5, 3, 1)+0.1)


#Getting total percent intermittent

setwd("D:/Charlotte/GIS")
Africa_preds<- read_sf(dsn="Africa_preds", layer="Africa_preds")
Asia_preds<- read_sf(dsn="Asia_preds", layer="Asia_preds")
Australia_preds<- read_sf(dsn="Australia_preds", layer="Australia_preds")
Europe_preds<- read_sf(dsn="Europe_preds", layer="Europe_preds")
NorthAmerica_preds<- read_sf(dsn="NorthAmerica_preds", layer="NorthAmerica_preds")
Northern_preds<- read_sf(dsn="Northern_preds", layer="Northern_preds")
SouthAmerica_preds<- read_sf(dsn="SouthAmerica_preds", layer="SouthAmerica_preds")
Siberia_preds<- read_sf(dsn="Siberia_preds", layer="Siberia_preds")

setwd('D:/Charlotte/Reach_R_data')
Africa_reaches <- readRDS("Africa_reaches.rds")
Asia_reaches <- readRDS("Asia_reaches.rds")
Australia_reaches <- readRDS("Australia_reaches.rds")
Europe_reaches <- readRDS("Europe_reaches.rds")
NorthAmerica_reaches <- readRDS("NorthAmerica_reaches.rds")
Northern_reaches <- readRDS("Northern_reaches.rds")
SouthAmerica_reaches <- readRDS("SouthAmerica_reaches.rds")
Siberia_reaches <- readRDS("Siberia_reaches.rds")


#Getting percent intermittent per continent and globally
total_percent <- sum(sum(Africa_preds$Class), sum(Asia_preds$Class),sum(Australia_preds$Class),
                     sum(Europe_preds$Class), sum(NorthAmerica_preds$Class), 
                     sum(Northern_preds$Class), sum(SouthAmerica_preds$Class), 
                     sum(Siberia_preds$Class))/sum(length(Africa_preds$Class), 
                                                   length(Asia_preds$Class),
                                                   length(Australia_preds$Class),
                                                   length(Europe_preds$Class), 
                                                   length(NorthAmerica_preds$Class), 
                                                   length(Northern_preds$Class), 
                                                   length(SouthAmerica_preds$Class), 
                                                   length(Siberia_preds$Class))




#Getting percent intermittent per climate zone 

climate_zone_analyst <- function(continent){
  setwd("D:/Charlotte/GIS")
  name <- paste(continent, "_preds", sep ="")
  continent_preds<- read_sf(dsn=name, layer=name)
  zonal_intermittent <- integer(18)
  zonal_total <- integer(18)
  pb <- progress_bar$new(total = length(continent_preds$Class))
  for(i in 1:length(continent_preds$Class)){
    if(continent_preds$Class[i] ==1){
      zonal_intermittent[continent_preds$Climtzn[i]] <- as.numeric(zonal_intermittent[continent_preds$Climtzn[i]]) +1
    }
    zonal_total[continent_preds$Climtzn[i]] <- as.numeric(zonal_total[continent_preds$Climtzn[i]]) +1
    pb$tick()
    Sys.sleep(1 / length(continent_preds$Class))
  }
  zonal <- rbind(zonal_intermittent, zonal_total)
  rownames(zonal) <- c("Intermittent", "Total")
  colnames(zonal) <- c(1:18)
  return(zonal)
}

Af_zonal <- climate_zone_analyst("Africa")
Asia_zonal <- climate_zone_analyst("Asia")
Aus_zonal <- climate_zone_analyst("Australia")
Eur_zonal <- climate_zone_analyst("Europe")
Na_zonal <- climate_zone_analyst("NorthAmerica")
North_zonal <- climate_zone_analyst("Northern")
Sa_zonal <- climate_zone_analyst("SouthAmerica")
Sib_zonal <- climate_zone_analyst("Siberia")

total_zonal <- (Af_zonal[2,]+Asia_zonal[2,]+Aus_zonal[2,]+Eur_zonal[2,]+
                  Na_zonal[2,]+North_zonal[2,]+Sa_zonal[2,]+Sib_zonal[2,])


zonal_percent_intermittent <- ((Af_zonal[1,]+Asia_zonal[1,]+Aus_zonal[1,]+Eur_zonal[1,]+
                                  Na_zonal[1,]+North_zonal[1,]+Sa_zonal[1,]+Sib_zonal[1,])/
                                 (Af_zonal[2,]+ Asia_zonal[2,]+Aus_zonal[2,]+Eur_zonal[2,]+
                                    Na_zonal[2,]+North_zonal[2,]+Sa_zonal[2,]+Sib_zonal[2,]))*100
par(mfrow=c(1,2)) 
barplot(total_zonal,  
        xlab="Climate zone",
        ylab = "Number of river reaches",
        cex.lab = 1.4
)
barplot(zonal_percent_intermittent, xlab = "Climate zone", ylab = "Percent intermittent", ylim = c(0, 100),  cex.lab = 1.4)


### Predicting France 
setwd('D:/Charlotte/IntermittenStudyFrance')
French_stations <- read.csv("French_RF_data.csv")


setwd('D:/Charlotte/Current_datasets/Reach_R_data')
France <- readRDS("France.rds")
setwd('D:/Charlotte/')
France_reaches <- read_sf(dsn="GIS", layer="France_lines")
newData <- France[,2:32]
French_predictions <- predict(RF_model, newdata =newData, type="prob", progress= "window")
French_preds <- as.data.frame(cbind(France[,1], ifelse(French_predictions[,2]<0.4,0,1), French_predictions[,2]))
colnames(French_preds)<- c("Reach_ID", "Class", "Probability")
st_geometry(French_preds) <- France_reaches$geometry
setwd("D:/Charlotte/GIS")
st_write(French_preds, dsn = "French_preds", driver ='ESRI Shapefile')


our_percent <- sum(French_preds$Class)/length(French_preds$Class)


setwd('D:/Charlotte/IntermittentStudyFrance/transfer_1487122_files_8a23b49b/')
French_reaches <- read_sf(dsn="rhtvs2_all", layer="rhtvs2_all_phi_qclass")
setwd('D:/Charlotte/IntermittentStudyFrance/')
french_classes <- read.delim("INT_RF.txt")
french_classes <- french_classes[french_classes$AllPred.ID_DRAIN %in% French_reaches$ID_DRAIN,]
French_reaches <- French_reaches[French_reaches$ID_DRAIN %in% french_classes$AllPred.ID_DRAIN,]
french_classes <- french_classes[!duplicated(french_classes$AllPred.ID_DRAIN),c(1:2)]
colnames(french_classes)<- c("Drain_ID", "Class")
st_geometry(french_classes) <- French_reaches$geometry
setwd("D:/Charlotte/GIS")
st_write(french_classes, dsn = "France_study_classes", driver ='ESRI Shapefile')

their_percent <- sum(french_classes$Class)/length(french_classes$Class)










