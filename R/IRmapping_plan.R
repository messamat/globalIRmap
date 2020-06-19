plan <- drake_plan(
  filestructure = def_filestructure(),

  monthlydischarge = read_monthlydis(in_filestructure = filestructure),

  gaugep = read_gaugep(in_filestructure = filestructure, dist = 200,
                       in_monthlydischarge = monthlydischarge),

  gauged_filenames = read_gauged_paths(filestructure, gaugep),

  gaugestats = future_map(gauged_filenames, #To map
                          comp_durfreq, #Function to run on each file name
                          maxgap = 20, monthsel = NULL, #Other arguments in function
                          .progress = TRUE),

  gaugestats_format = format_gaugestats(gaugestats, gaugep),

  predvars = selectformat_predvars(filestructure, in_gaugestats = gaugestats_format),

  tasks = create_tasks(in_gaugestats = gaugestats_format, in_predvars = predvars),

  baselearners = target(
    create_baselearners(tasks),
    dynamic = map(tasks)),

  #Subdivide the baselearners dynamic targets from 3 to 8
  seplearners = target(readd(baselearners, subtarget_list = FALSE)),

  measures = target(list(classif = msr("classif.bacc"),
                         regr = msr("regr.mae"))
  ),

  autotuningset = target(
    set_tuning(in_learner = seplearners,
               in_measures = measures,
               nfeatures = length(tasks$classif$feature_names),
               insamp_nfolds = 5, insamp_neval = 100,
               insamp_nbatch = parallel::detectCores(logical=FALSE)-1
    ),
    dynamic = map(seplearners)
  ),

  resamplingset = set_cvresampling(rsmp_id = 'repeated_cv',
                                   in_task = tasks$classif,
                                   outsamp_nrep = 2,
                                   outsamp_nfolds = 3),

  rfresampled_classif = target(
    dynamic_resample(in_task = tasks$classif,
                     in_learner = autotuningset,
                     in_resampling = resamplingset,
                     store_models = TRUE,
                     type = 'classif'),
    dynamic = map(autotuningset)
  ),

  rfresampled_regr = target(
    dynamic_resample(in_task = in_tsk,
                     in_learner = in_lrn,
                     in_resampling = resamplingset,
                     store_models = TRUE,
                     type = 'regr'),
    transform = map(in_tsk = c(tasks$regr, tasks$regover),
                    in_lrn = c(autotuningset[[7]], autotuningset[[8]]),
                    .names = c('rfresampled_regr_res', 'rfresampled_regover_res'))

  ),

  rfbm_classif = target(
    combine_bm(readd(rfresampled_classif, subtarget_list = TRUE),
               out_qs = file_out(!!file.path(readd(filestructure)[['resdir']],
                                             'rfbm_classif.qs')))
  ),

  # bm_checked = target(
  #   analyze_benchmark(in_bm, in_measure = measures),
  #   transform = map(in_bm = list(file_in(!!file.path(readd(filestructure)[['resdir']],
  #                                                     'rfbm_classif.qs')),
  #                                c(rfresampled_regr_res, rfresampled_regover_res)))#,
  #                 #in_measure = list(measures$classif, measures$regr))
  # ),

  resamplingset_featsel = target(
    set_cvresampling(rsmp_id = in_strategy,
                     in_task = tasks$classif,
                     outsamp_nrep = in_outrep,
                     outsamp_nfolds = in_outfolds),
    transform = map(in_strategy = c('repeated_cv', "repeated-spcv-coords"),
                    in_outrep = c(1, 1),
                    in_outfolds = c(2, 2),
                    .names = c('featsel_cv', 'featsel_spcv'))
  ),

  tasks_featsel = select_features(
    in_bm = file_in(!!file.path(readd(filestructure)[['resdir']],
                                'rfbm_classif.qs')),
    in_lrnid =  "oversample.classif.ranger",
    in_task = tasks$classif,
    pcutoff = 0.05
  ),

  rfresampled_featsel = target(
    dynamic_resamplebm(in_task = in_taskfeatsel,
                       in_bm = file_in(!!file.path(readd(filestructure)[['resdir']],
                                                   'rfbm_classif.qs')),
                       in_lrnid =  "oversample.classif.ranger",
                       in_resampling = in_resampling,
                       store_models = TRUE,
                       type = 'classif'),
    transform= cross(in_taskfeatsel = c(tasks_featsel[[1]], tasks_featsel[[2]]),
                     in_resampling = c(featsel_cv, featsel_spcv),
                     .names = c('res_all_cv', 'res_all_spcv',
                                'res_featsel_cv', 'res_featsel_spcv'))
  ),

  rfeval_featsel = target(
    c(res_all_cv, res_all_spcv,
      res_featsel_cv, res_featsel_spcv), #Cannot use combine as lead to BenchmarkResult directly in the branching
  ),

  rfbm_featsel = analyze_benchmark(in_bm = rfeval_featsel,
                                   in_measure = measures$classif),

  #  Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
  rftuned = target(
    selecttrain_rf(in_rf = rfeval_featsel,
                   in_learnerid ="oversample.classif.ranger",
                   in_taskid = "inter_basicsp_featsel",
                   insamp_nfolds = 2,
                   insamp_nevals = 1)),

  rivernetwork = rformat_network(in_filestructure = filestructure,
                                 in_predvars = predvars,
                                 in_monthlydischarge = monthlydischarge),

  uncertainty_plot = gguncertainty(in_rftuned = rftuned,
                                   in_gaugestats = gaugestats_format,
                                   in_predvars = predvars,
                                   spatial_rsp = TRUE),

  rfpreds = write_preds(in_filestructure = filestructure, in_gaugep = gaugep,
                        in_gaugestats = gaugestats_format, in_network = rivernetwork,
                        in_rftuned = rftuned, in_predvars = predvars,
                        in_uncertainty = uncertainty_plot[['out_uncertainty']]),

  vimp_plot = ggvimp(rftuned, predvars, varnum=20, spatial_rsp = FALSE),

  pd_plot = ggpd(in_rftuned=rftuned, in_predvars=predvars, colnums=1:10,
                 nodupli = TRUE, ngrid = 20, parallel = T, spatial_rsp = FALSE),

  basemaps = get_basemapswintri(in_filestructure = filestructure),

  gauges_plot = gggauges(in_gaugepred = rfpreds, in_basemaps = basemaps),

  envhist = layout_ggenvhist(in_rivernetwork = rivernetwork,
                             in_gaugepred = rfpreds,
                             in_predvars = predvars),

  table_allbm = target(
    tabulate_benchmarks(in_bm, in_bmid),
    transform = map(
      in_bm = list(
        file_in(!!file.path(readd(filestructure)[['resdir']], 'rfbm_classif.qs')),
        c(rfresampled_regr_res, rfresampled_regover_res),
        rfeval_featsel),
      in_bmid = list('classif1', 'regr1', 'classif2'),
      .names = c('tablebm_classif1', 'tablebm_regr1', 'tablebm_classif2'))
  ),

  misclass_format = target(
    formatmisclass_bm(in_bm = in_bm, in_bmid = in_bmid),
    transform = map(
      in_bm = list(
        file_in(!!file.path(readd(filestructure)[['resdir']], 'rfbm_classif.qs')),
        c(rfresampled_regr_res, rfresampled_regover_res),
        rfeval_featsel),
      in_bmid = list('classif1', 'regr1', 'classif2'),
      .names = c('misclass_classif1', 'misclass_regr1', 'misclass_classif2'))
  ),

  misclass_plot = ggmisclass_bm(list(misclass_classif1,
                                     misclass_regr1,
                                     misclass_classif2)),

  krigepreds = krige_spuncertainty(in_filestructure = filestructure, in_rftuned = rftuned,
                                   in_gaugep = gaugep, in_gaugestats = gaugestats_format,
                                   kcutoff=50000, overwrite=T),

  krigepreds_mosaic = mosaic_kriging(in_filestructure = filestructure,
                                     in_kpathlist = krigepreds,
                                     overwrite = TRUE),

  globaltables = target(
    tabulate_globalsummary(in_filestructure = filestructure,
                           idvars = in_idvars,
                           castvar = 'ORD_STRA', castvar_num = TRUE,
                           weightvar = 'LENGTH_KM',
                           valuevar = 'predbasic800cat',
                           valuevarsub = 1,
                           na.rm=T, tidy = FALSE),
    transform = map(in_idvars = c('gad_id_cmj', 'fmh_cl_cmj',
                                  'tbi_cl_cmj', 'clz_cl_cmj'))
  ),

  fr_plot = compare_fr(in_filestructure = filestructure,
                       in_rivernetwork = rivernetwork,
                       binarg = c(10,20,50,100,200,500,1000,
                                  2000,5000,10000,50000,100000,150000)),

  us_plot = compare_us(in_filestructure = filestructure,
                       in_rivernetwork = rivernetwork,
                       binarg = c(10,20,50,100,200,500,1000,
                                  2000,5000,10000,50000,100000,2000000, 3200000)),

  pnw_plot = qc_pnw(in_filestructure = filestructure,
                    in_rivernetwork = rivernetwork)
)

