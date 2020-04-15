library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

in_gaugestats= readd(gaugestats_format)
in_predvars = readd(predvars)
insamp_nfolds = 2
insamp_neval = 20
insamp_nbatch = parallel::detectCores(logical=FALSE)
outsamp_nrep = 1
outsamp_nfolds = 5

#Create task
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

#Create mlr3 pipe operator to oversample minority class based on major/minor ratio
#https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
#https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
#Sampling happens only during training phase.
po_over <- po("classbalancing", id = "oversample", adjust = "minor", 
              reference = "minor", shuffle = TRUE, 
              ratio = get_oversamp_ratio(task_inter)$ratio)
#table(po_over$train(list(task_inter))$output$truth()) #Make sure that oversampling worked

#Create a graph learner so that oversampling happens systematically upstream of all training
lrn_ranger_over <- GraphLearner$new(po_over %>>% lrn_ranger)

#Define parameter ranges to tune on
regex_tuneset <- function(in_lrn) {
  prmset <- names(in_lrn$param_set$tags)
  
  tune_ranger <- ParamSet$new(list(
    ParamInt$new(prmset[grep(".*mtry", prmset)], 
                 lower = 1, upper = 11),
    ParamInt$new(prmset[grep(".*min.node.size", prmset)], 
                 lower = 1, upper = 10),
    ParamDbl$new(prmset[grep(".*sample.fraction", prmset)], 
                 lower = 0.1, upper = 0.95)
  ))
}

#Define inner resampling strategy
rcv_ranger = rsmp("cv", folds=insamp_nfolds) #5-fold aspatial CV repeated 10 times
#Define performance measure
measure_ranger = msr("classif.bacc") #use balanced accuracy to account for imbalanced dataset
#Define termination rule 
evals20 = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds

#Define hyperparameter tuner wrapper for inner sampling
learns = list(
  AutoTuner$new(learner= lrn_ranger,
                resampling = rcv_ranger, 
                measures = measure_ranger,
                tune_ps = regex_tuneset(lrn_ranger), 
                terminator = evals20,
                tuner =  tnr("random_search", 
                             batch_size = insamp_nbatch)), #batch_size determines level of parallelism
  
  AutoTuner$new(learner= lrn_ranger_over,
                resampling = rcv_ranger, 
                measures = measure_ranger,
                tune_ps = regex_tuneset(lrn_ranger_over), 
                terminator = evals20,
                tuner =  tnr("random_search", 
                             batch_size = insamp_nbatch)) 
)
names(learns) <-mlr3misc::map(learns, "id")

#Perform outer resampling, keeping models for diagnostics later
outer_resampling = rsmp("repeated_cv", 
                        repeats = outsamp_nrep, 
                        folds = outsamp_nfolds)

nestedresamp_bmrdesign <- benchmark_grid(tasks = task_inter, 
                                         learners = learns, 
                                         resamplings = outer_resampling)

nestedresamp_bmrout <- benchmark(nestedresamp_bmrdesign,
                                store_models = TRUE)

print(nestedresamp_bmrout$aggregate(measure_ranger))
mlr3viz::autoplot(nestedresamp_bmrout, measure = measure_ranger)

#Train learners
learns$classif.ranger.tuned$train(task_inter) 
learns$oversample.classif.ranger.tuned$train(task_inter)

#Return outer sampling object for selected model
nestedresamp_bmrout$resample_result(
  which(names(learns)=="oversample.classif.ranger.tuned"))

rfresamp <- rftuned$rf_outer


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