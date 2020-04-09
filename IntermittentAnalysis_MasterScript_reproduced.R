#### -------------------------- Run RF model -------------------------------------------------
#---- Tune model ----
# rf_formula <- as.formula(paste0('intermittent~', 
#                                 paste(rfpredcols, collapse="+"), 
#                                 collapse=""))

#Create taxk
task_inter <- mlr3::TaskClassif$new(id ='inter_basic',
                                    backend = gaugestats_join[!is.na(cly_pc_cav), c('intermittent', rfpredcols),with=F],
                                    target = "intermittent")

#Create learner
lrn_ranger <- mlr3::lrn('classif.ranger', num.trees=500, replace=F, predict_type="prob", importance = "permutation")
lrn_ranger$param_set

#Define parameter ranges to tune on
tune_ranger <- ParamSet$new(list(
  ParamInt$new("mtry", lower = 1, upper = 11),
  ParamInt$new("min.node.size", lower = 1, upper = 10),
  ParamDbl$new("sample.fraction", lower = 0.1, upper = 0.9)
))

#Define resampling strategy
rcv_ranger = rsmp("cv", folds=5) #5-fold aspatial CV repeated 10 times
#Define performance measure
measure_ranger = msr("classif.bacc") #use balanced accuracy to account for imbalanced dataset
#Define termination rule 
evals20 = term("evals", n_evals = 20) #termine tuning after 20 rounds


#Define hyperparameter tuner wrapper for inner sampling
at_ranger <- AutoTuner$new(learner= lrn_ranger,
                           resampling = rcv_ranger, 
                           measures = measure_ranger,
                           tune_ps = tune_ranger, 
                           terminator = evals20,
                           tuner =  tnr("random_search", batch_size = 4L))

#Launch nested resampling (outer resampling)
future::plan("multiprocess")
tic()
nestedresamp_ranger <- resample(task = task_inter, 
                                learner = at_ranger, 
                                resampling = rsmp("repeated_cv", repeats = 10, folds = 5),
                                store_models = TRUE)
toc()
nestedresamp_ranger$aggregate()


lapply(nestedresamp_ranger$data$learner, function(x) x$variable.importance)

check <- nestedresamp_ranger$data$learner[[1]]$model$learner$model



#### -------------------------- Diagnose initial model -------------------------------------------------
# ---- Variable importance ----

varimp_basic <- data.table(varcode=names(ranger_basictrained$variable.importance), 
                           importance=ranger_basictrained$variable.importance)[
                             rfpredcols_dt, on='varcode' ]%>%
  .[, varname := factor(varname, levels=.[,varname[order(-importance)]])] %>%
  setorder(-importance)

ggplot(varimp_basic[1:30,],aes(x=varname, y=importance)) + 
  geom_bar(stat = 'identity') +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=8))

# ---- Check misclassification rate with different threshold probabilities ----
gaugestats_join[!is.na(cly_pc_cav), intermittent_predprob := 
                  as.data.table(nestedresamp_ranger$prediction())[order(row_id), mean(prob.1), by=row_id]$V1]


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

threshold_misclass[i %in% c(0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.35, 0.5, 0.45, 0.5),]
gaugestats_join[, intermittent_predcat := ifelse(intermittent_predprob>=0.28, 1, 0)]

# ---- Partial dependence plots ----
pdtiles_grid(mod = ranger_basictrained, dt = gaugestats_join, colnums = 5)

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
           dis_pc_pvar=ifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx), #min/max monthly watershed discharge
           dis_pc_pvaryr=ifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),#min monthly/average yearly watershed discharge
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
             predbasic800 := predict(ranger_basictrained, type='response', data=.SD)$predictions[,'1'], 
             .SDcols = c("HYRIV_ID", rfpredcols)]
  toc()
}

riveratlas[, predbasic800cat := ifelse(predbasic800>=0.3, 1, 0)]
fwrite(riveratlas[, c('HYRIV_ID', 'predbasic800cat'), with=F], file.path(resdir, 'RiverATLAS_predbasic800.csv'))




















##################################################################################################################################################
############################ MORE ALREADY SPED UP/TRANSCRIBED ####################################################################################
#Set default ranger learner with explicit parameter set
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