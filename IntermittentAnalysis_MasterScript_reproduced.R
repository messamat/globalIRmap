library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

in_gaugestats= readd(gaugestats_format)
in_predvars = readd(predvars)
in_rftuned = readd(rftuned)
insamp_nfolds = 2
insamp_neval = 20
insamp_nbatch = parallel::detectCores(logical=FALSE)
outsamp_nrep = 1
outsamp_nfolds = 5

datsel <- in_gaugestats[!is.na(cly_pc_cav), 
                        c('intermittent',in_predvars$varcode),
                        with=F]
task_inter <- mlr3::TaskClassif$new(id ='inter_basic',
                                    backend = datsel,
                                    target = "intermittent")

#Create learner
lrn_ranger <- mlr3::lrn('classif.ranger', 
                        num.trees=500, 
                        replace=F, 
                        predict_type="prob", 
                        importance = "permutation")
#print(lrn_ranger$param_set)

tune_ranger <- ParamSet$new(list(
  ParamInt$new("mtry", lower = 1, upper = 11),
  ParamInt$new("min.node.size", lower = 1, upper = 10),
  ParamDbl$new("sample.fraction", lower = 0.1, upper = 0.95)
))

#Define inner resampling strategy
rcv_ranger = rsmp("cv", folds=insamp_nfolds) #5-fold aspatial CV repeated 10 times
#Define performance measure
measure_ranger = msr("classif.bacc") #use balanced accuracy to account for imbalanced dataset
#Define termination rule 
evals20 = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds

at_ranger <- AutoTuner$new(learner= lrn_ranger,
                           resampling = rcv_ranger, 
                           measures = measure_ranger,
                           tune_ps = tune_ranger, 
                           terminator = evals20,
                           tuner =  tnr("random_search", 
                                        batch_size = insamp_nbatch))

nestedresamp_ranger <- mlr3::resample(task = task_inter, 
                                      learner = at_ranger, 
                                      resampling = rsmp("repeated_cv", 
                                                        repeats = 1, 
                                                        folds = 2),
                                      store_models=TRUE)

loadd(rftuned)
in_mod <- rftuned$rf_outer$learners[[1]] 
in_mod <- nestedresamp_ranger$learners[[1]]


#in_mod <- nestedresamp_ranger$learners[[1]] 

if (inherits(in_mod$learner, "GraphLearner")) {
  in_fit <- in_mod$learner$model$classif.ranger$model
} else {
  in_fit <- in_mod$learner$model 
}
class(in_fit)
in_fit

foldperf <- extract_impperf_nestedrf(in_mod, imp=F, perf=T)

# selcols <- in_vimp_plot$data %>% #Can use that if extracting from tunredrf is expensive
#   setorder(-imp_wmean) %>%
#   .[colnums, variable]
ngrid <- c(10,10)
selcols <- c('run_mm_cyr', 'dis_m3_pmn')
datdf <- as.data.frame(task_inter$data())
rftuned$rf_outer$task$data()
nestedresamp_ranger$task$data()

pdout <- edarf::partial_dependence(in_fit, vars = selcols, n = ngrid,
                                   interaction = TRUE, data = datdf) #Warning: does not work with data_table
  
  
  

rftuned$rf_outer$aggregate(msr('time_both'))

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