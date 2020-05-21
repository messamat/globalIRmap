plan <- drake_plan(
  filestructure = def_filestructure(),
  
  gaugep = read_gaugep(in_filestructure = filestructure, dist = 200),
  
  gauged_filenames = read_gauged_paths(filestructure, gaugep),
  
  gaugestats = future_map(gauged_filenames, #To map
                          comp_durfreq, #Function to run on each file name
                          maxgap = 20, monthsel = NULL, #Other arguments in function
                          .progress = TRUE),
  
  gaugestats_format = format_gaugestats(gaugestats, gaugep),
  
  predvars = selectformat_predvars(filestructure, in_gaugestats = gaugestats_format),
  
  tasks = create_tasks(in_gaugestats = gaugestats_format, in_predvars = predvars), 
  
  rfbm_classif = benchmark_classif(in_tasks = tasks, 
                                   insamp_nfolds = 2, insamp_neval = 5, 
                                   insamp_nbatch = parallel::detectCores(logical=FALSE),
                                   outsamp_nrep = 1, outsamp_nfolds = 2),
  
  rfbm_regr = benchmark_regr(in_tasks = tasks, 
                             insamp_nfolds = 2, insamp_neval = 5, 
                             insamp_nbatch = parallel::detectCores(logical=FALSE),
                             outsamp_nrep = 1, outsamp_nfolds = 2),
  
  bm_checked = target(
    analyze_benchmark(in_bm, in_measure),
    transform = map(in_bm = c(rfbm_classif$bm_classif, rfbm_regr$bm_regr),
                    in_measure = c(rfbm_classif$measure_classif, rfbm_regr$meassure_regr))
  ),
  
  rfeval_featsel = target(
    benchmark_featsel(in_rf = rfbm_classif$bm_classif$clone()$filter(learner_ids = "oversample.classif.ranger.tuned"),
                      in_task = rfbm_classif$bm_tasks$task_classif,
                      in_measure = rfbm_classif$measure_classif,
                      pcutoff = 0.1,
                      insamp_nfolds =  2, insamp_nevals = 1,
                      outsamp_nrep = 2, outsamp_nfolds =  10, 
                      outsamp_nfolds_sp = 10) 
  ),
  
  #  Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
  rftuned = target(
    selecttrain_rf(in_rf = rfeval_featsel$bm_classif$clone()$filter(learner_ids = "oversample.classif.ranger.tuned"),
                   in_task = rfeval_featsel$bm_tasks$task_classif,
                   insamp_nfolds = 2,
                   insamp_nevals = 2),
    trigger =  trigger(condition = TRUE)),
  
  misclass_plot = ggmisclass(in_predictions=rftuned$rf_outer$prediction()),
  
  vimp_plot = ggvimp(rftuned, predvars),
  
  pd_plot = ggpd(rftuned, predvars, colnums=1:5, ngrid=c(10,10), parallel=T),
  
  rfpreds = write_preds(filestructure, gaugep, gaugestats_format, rftuned, predvars)
)
                    