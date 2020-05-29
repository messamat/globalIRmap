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

  baselearners = target(
    create_baselearners(in_task = taskbase),
    transform = map(taskbase = list(tasks$classif, tasks$regr),
                    .names= c('baselearners.classif', 'baselearners.regr'))
  ),

  measures = target(list(classif = msr("classif.bacc"),
                         regr = msr("regr.mae"))
  ),

  ############## NOT WORKING YET #########################
  tuningset = target(
    set_tuning(in_learner = learnertune,
               in_measures = measures,
               nfeatures = length(tasks$classif$feature_names),
               insamp_nfolds = 5, insamp_neval = 100,
               insamp_nbatch = parallel::detectCores(logical=FALSE)),
    transform = map(learnertune= unlist(baselearners, recursive=F))
  )
)

,
#
#   resamplingset = set_cvresampling(rsmp_id = 'repeated_cv',
#                                    in_task = tasks$classif,
#                                    outsamp_nrep = 1,
#                                    outsamp_nfolds = 2),
#
#   rfresampled_classif = target(
#     resample(task = tasks$classif,
#              learner =learners_classif,
#              resampling = resamplingset,
#              store_models = TRUE),
#     transform = map(learners_classif = baselearners[[1]])
#   ),
#
#   rfresampled_regr = target(
#     resample(task = tasks$regr, learner = learners_regr, resampling = resamplingset,
#              store_models = TRUE),
#     transform = map(learners_regr = baselearners[[2]])
#   ),
#
#   rfbm_classif = target(
#     combine_bm(in_resampleresults = rsmps_classif),
#     transform = map(rsmps_classif = rfresampled_classif)
#   ),
#
#   rfbm_regr = target(
#     combine_bm(in_resampleresults = rsmps_regr),
#     transform = map(rsmps_regr = rfresampled_regr)
#   )
# )


#   bm_checked = target(
#     analyze_benchmark(in_bm = in_bm, in_measure = measures$classif),
#     dynamic = map(in_bm = list(rfbm_classif, rfbm_regr))#,
#                   #in_measure = list(measures$classif, measures$regr))
#   ),
#
#   rfeval_featsel = target(
#     benchmark_featsel(in_rf = rfbm_classif$clone()$filter(learner_ids = "oversample.classif.ranger.tuned"),
#                       in_task = tasks$task_classif,
#                       in_measure = measures$classif,
#                       pcutoff = 0.05,
#                       insamp_nfolds =  2, insamp_nevals = 1,
#                       outsamp_nrep = 1, outsamp_nfolds =  2,
#                       outsamp_nfolds_sp = 2)
#   ),
#
#   #  Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
#   rftuned = target(
#     selecttrain_rf(in_rf = rfeval_featsel$bm_classif$clone()$filter(learner_ids = "oversample.classif.ranger.tuned"),
#                    in_task = tasks$task_classif,
#                    insamp_nfolds = 2,
#                    insamp_nevals = 1)),
#
#   misclass_plot = ggmisclass(in_rftuned = rftuned, spatial_rsp = FALSE),
#
#   vimp_plot = ggvimp(rftuned, predvars, varnum=20, spatial_rsp = FALSE),
#
#   pd_plot = ggpd(in_rftuned=rftuned, in_predvars=predvars, colnums=1:10,
#                  nodupli = TRUE, ngrid = 20, parallel = T, spatial_rsp = FALSE),
#
#   uncertainty_plot = gguncertainty(in_rftuned = rftuned,
#                                    in_gaugestats = gaugestats_format,
#                                    in_predvars = predvars,
#                                    spatial_rsp = FALSE),
#
#   rfpreds = write_preds(filestructure, gaugep, gaugestats_format,
#                         rftuned, predvars)
# )

