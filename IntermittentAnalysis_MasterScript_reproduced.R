library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

in_gaugestats <- readd(gaugestats_format)
in_predvars <- readd(predvars)
insamp_nfolds = 2
insamp_neval = 5
insamp_nbatch = parallel::detectCores(logical=FALSE)
outsamp_nrep = 1
outsamp_nfolds = 2

benchmark_rf <- function(in_gaugestats, in_predvars, 
                         insamp_nfolds, insamp_neval, insamp_nbatch,
                         outsamp_nrep, outsamp_nfolds) {
  
  #---------- Create tasks -----------------------------------------------------
  datsel <- in_gaugestats[!is.na(cly_pc_cav), 
                          c('intermittent',in_predvars$varcode),
                          with=F]
  
  task_classif <- mlr3::TaskClassif$new(id ='inter_basic',
                                        backend = datsel,
                                        target = "intermittent")
  
  task_regr <- convert_clastoregrtask(in_task = task_classif,
                                      in_id = 'inter_regr',
                                      oversample=FALSE) 
  
  task_regrover <- convert_clastoregrtask(in_task = task_classif,
                                          in_id = 'inter_regrover',
                                          oversample=TRUE) 
  
  #---------- Create learners --------------------------------------------------
  #Create basic learner
  lrn_ranger <- mlr3::lrn('classif.ranger', 
                          num.trees = 500, 
                          replace = FALSE, 
                          splitrule = 'gini',
                          predict_type = "prob", 
                          importance = "permutation",
                          respect.unordered.factors = 'order')
  
  #print(lrn_ranger$param_set)
  
  #Create learner with 
  lrn_ranger_maxstat <- mlr3::lrn('regr.ranger', 
                                  num.trees=500, 
                                  replace=FALSE, 
                                  splitrule = 'maxstat',
                                  importance = "permutation",
                                  respect.unordered.factors = 'order')
  
  #Create mlr3 pipe operator to oversample minority class based on major/minor ratio
  #https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
  #https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
  #Sampling happens only during training phase.
  po_over <- po("classbalancing", id = "oversample", adjust = "minor", 
                reference = "minor", shuffle = TRUE, 
                ratio = get_oversamp_ratio(task_classif)$ratio)
  #table(po_over$train(list(task_classif))$output$truth()) #Make sure that oversampling worked
  
  #Create a graph learner so that oversampling happens systematically upstream of all training
  lrn_ranger_overp <- GraphLearner$new(po_over %>>% lrn_ranger)
  
  #---------- Set up inner resampling ------------------------------------------
  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)
    
    tune_ranger <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T), 
                   lower = 1, upper = 11),
      ParamDbl$new(grep(".*sample.fraction", prmset, value=T), 
                   lower = 0.2, upper = 0.8)
    ))
    
    in_splitrule =in_lrn$param_set$get_values()[
      grep(".*splitrule", prmset, value=T)]
    
    if (in_splitrule == 'maxstat') {
      tune_ranger$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)], 
                     lower = 0.01, upper = 0.1)
      )
      
    } else {
      tune_ranger$add(
        ParamInt$new(prmset[grep(".*min.node.size", prmset)], 
                     lower = 1, upper = 10)
      )
    }
  }
  
  #Define inner resampling strategy
  rcv_ranger = rsmp("cv", folds=insamp_nfolds) #5-fold aspatial CV repeated 10 times
  
  #Define performance measure
  measure_ranger_class = msr("classif.bacc") #use balanced accuracy as objective function
  measure_ranger_reg = msr("regr.mae") 
  
  #Define termination rule 
  evalsn = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds
  
  #Define hyperparameter tuner wrapper for inner sampling
  learns_classif = list(
    AutoTuner$new(learner= lrn_ranger,
                  resampling = rcv_ranger, 
                  measures = measure_ranger_class,
                  tune_ps = regex_tuneset(lrn_ranger), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch)), #batch_size determines level of parallelism
    
    AutoTuner$new(learner= lrn_ranger_overp,
                  resampling = rcv_ranger, 
                  measures = measure_ranger_class,
                  tune_ps = regex_tuneset(lrn_ranger_overp), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch))
  )
  names(learns_classif) <-mlr3misc::map(learns_classif, "id")
  
  learns_regr = list(
    AutoTuner$new(learner= lrn_ranger_maxstat,
                  resampling = rcv_ranger, 
                  measures = measure_ranger_reg,
                  tune_ps = regex_tuneset(lrn_ranger_maxstat), 
                  terminator = evalsn,
                  tuner =  tnr("random_search", 
                               batch_size = insamp_nbatch)))
  
  #---------- Set up outer resampling benchmarking -----------------------------
  
  #Perform outer resampling, keeping models for diagnostics later
  outer_resampling = rsmp("repeated_cv", 
                          repeats = outsamp_nrep, 
                          folds = outsamp_nfolds)
  
  #Run outer resampling and benchmarking on classification learners
  nestedresamp_bmrdesign_classif <- benchmark_grid(
    tasks = task_classif, 
    learners = learns_classif, 
    resamplings = outer_resampling)
  
  nestedresamp_bmrout_classif <- benchmark(
    nestedresamp_bmrdesign_classif, store_models = TRUE)
  
  #Run outer resampling and benchmarking on regression learners
  nestedresamp_bmrdesign_regr <- benchmark_grid(
    tasks = list(task_regr, task_regrover), 
    learners = learns_regr, 
    resamplings = outer_resampling)
  
  nestedresamp_bmrout_regr <- benchmark(
    nestedresamp_bmrdesign_regr, store_models = TRUE)
  
  return(
    list(
      bm_classif = nestedresamp_bmrout_classif,
      bm_regr = nestedresamp_bmrout_regr,
      bm_tasks = list(task_classif=task_classif, 
                      task_regr=task_regr,
                      task_regrover=task_regrover),
      measure_classif = measure_ranger_class,
      measure_regr = measure_ranger_reg
    )
  )
}
################################################
rfbm <- list(
  bm_classif = nestedresamp_bmrout_classif,
  bm_regr = nestedresamp_bmrout_regr,
  bm_tasks = list(task_classif, task_regr, task_regrover),
  measure_classif = measure_ranger_class,
  measure_regr = measure_ranger_reg
)
in_bm <- rfbm$bm_classif
in_measure <- rfbm$measure_classif

###################################################

################################################################################
#in_rf <- 
in_rf <- rfbm$bm_classif$filter(learner_ids = "oversample.classif.ranger.tuned")
check <- in_rf$learners$learner[[1]]
in_rf$resamplings$resampling


selecttrain_rf <- function(in_rf, in_folds =  NULL, in_nevals = NULL) {
  lrn_autotuner <- in_rf$learners$learner[[1]]
  
  if (!is.null(in_folds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- in_folds
  }
  
  if (!is.null(in_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- in_nevals
  }
  
  #Train learners
  in_rf$learners$learner[[1]]$train(rfbm$bm_tasks$task_classif)
  
  #Return outer sampling object for selected model
  outer_resampling_output <- in_rf$resample_result(1)
  
  
  
  return(list(rf_outer = outer_resampling_output, #Resampling results
              rf_inner = in_rf$learners$learner[[1]], #Core learner (with hyperparameter tuning)
              task_basic = in_rf$tasks$task[[1]])) #Task
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