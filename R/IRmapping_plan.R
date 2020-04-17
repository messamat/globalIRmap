plan <- drake_plan(
  filestructure = def_filestructure(),
  
  gaugep = read_gaugep(in_filestructure = filestructure, dist = 200),
  
  gauged_filenames = read_gauged_paths(filestructure, gaugep),
  
  gaugestats = future_map(gauged_filenames, #To map
                          comp_durfreq, #Function to run on each file name
                          maxgap = 20, monthsel = NULL, #Other arguments in function
                          .progress = TRUE),
  
  gaugestats_format = format_gaugestats(gaugestats, gaugep),
  
  predvars = selectformat_predvars(filestructure, gaugestats_format),
  
  rfbm = benchmark_rf(in_gaugestats= gaugestats_format, in_predvars = predvars, 
                      insamp_nfolds = 2, insamp_neval = 20, 
                      insamp_nbatch = parallel::detectCores(logical=FALSE),
                      outsamp_nrep = 2, outsamp_nfolds = 5),
  
  #bm_ = target() , #Add analyze_benchmark with dynamic variables with in_bm, in_tasks, in_measure
  
  rftuned = selecttrain_rf(
    in_rf = rfbm$bm_classif$filter(learner_ids = "oversample.classif.ranger.tuned"),
    in_task = rfbm$bm_tasks$task_classif,
    insamp_nfolds = 5,
    insamp_nevals = 50), 
  
  misclass_plot = ggmisclass(gaugestats_format, rftuned$rf_outer$prediction()),
  
  vimp_plot = ggvimp(rftuned, predvars),
  
  pd_plot = ggpd(rftuned, predvars, colnums=1:5, ngrid=c(10,10), parallel=T),
  
  rfpreds = write_preds(filestructure, gaugep, gaugestats_format, rftuned, predvars)
)
                    