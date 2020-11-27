source('R/IRmapping_functions.R')
source('R/planutil_functions.R')

rootdir <- find_root(has_dir("src")) #Set root working directory

########################### plan_preprocess ####################################
plan_preprocess <- drake_plan(
  #-------------------- Define input and output files ##########################
  #Note: for dependency tracking, drake only reads string literals inside file_in and file_out; hence the use of absolute paths
  #Tidy evaluation (!!) didn't work for some reason -- can look into it later
  path_resdir = file.path(!!rootdir, "results"),
  path_GRDCgaugep = file.path(path_resdir, "spatialoutputs.gdb\\grdcstations_cleanjoin"),
  path_GRDCgaugedir =  file_in(!!file.path(rootdir, "data\\GRDCdat_day")),
  path_GSIMgaugep = file.path(path_resdir, "GSIM\\GSIM.gdb\\GSIMstations_cleanriverjoin"),
  path_GSIMindicesdir = file.path(!!rootdir, "data\\GSIM\\GSIM_indices"),
  path_riveratlas_meta = file_in(
    !!file.path(rootdir, "data\\HydroATLAS\\HydroATLAS_metadata_MLMv11.xlsx")),
  path_riveratlas_legends = file_in(!!file.path(rootdir, "data\\HydroATLAS\\HydroATLAS_v10_Legends.xlsx")),
  path_monthlynetdischarge = file.path(!!rootdir,
                                       "data\\HydroSHEDS\\HS_discharge_monthly.gdb\\Hydrosheds_discharge_monthly"),
  path_riveratlas = file_in(!!file.path(rootdir, 'results\\RiverATLAS_v10tab.csv')),
  path_riveratlas2 = file_in(!!file.path(rootdir, 'results\\RiverATLAS_v11tab.csv')),
  path_bas03 = file.path(!!rootdir, "data\\HydroATLAS\\BasinATLAS_v10.gdb\\BasinATLAS_v10_lev03"),
  path_compresdir = file_in(!!file.path(rootdir, 'results\\Comparison_databases')),
  path_frresdir = file.path(path_compresdir, 'france.gdb'),
  path_usdatdir = file.path(rootdir, 'data\\Comparison_databases', 'US'),
  path_usresdir = file.path(path_compresdir, 'us.gdb'),
  path_insitudatdir = file.path(rootdir, 'data\\Insitu_databases'),
  path_insituresdir = file.path(path_resdir, 'Insitu_databases'),
  path_pnwresdir =  file.path(path_insituresdir, 'pnw.gdb'),
  path_ondedatdir = file.path(path_insitudatdir, 'OndeEau'),
  path_onderesdir = file.path(path_insituresdir, 'ondeeau.gdb'),

  outpath_gaugep = file_out(!!file.path(rootdir, 'results\\GRDCstations_predbasic800.gpkg')),
  outpath_riveratlaspred = file_out(!!file.path(rootdir, 'results\\RiverATLAS_predbasic800.csv')),
  outpath_bas03error = file_out(!!file.path(rootdir, 'results\\BasinATLAS_v10_lev03_errors.gpkg')),
  path_bufrasdir = file.path(path_resdir, 'bufrasdir')
  ,

  outpath_krigingtif = file_out("E:\\Mathis\\results\\prederror_krigingtest.tif"),

  #-------------------- Pre-analysis ------------------------------------------
  #monthlydischarge = read_monthlydis(in_path = path_monthlynetdischarge),

  gaugep = target(
    read_gaugep(inp_GRDCgaugep = path_GRDCgaugep,
                inp_GSIMgaugep = path_GSIMgaugep,
                #in_monthlydischarge = monthlydischarge,
                inp_riveratlas2 = path_riveratlas2
    )
  )
  ,

  GRDCgauged_filenames = read_GRDCgauged_paths(
    inp_GRDCgaugedir = path_GRDCgaugedir,
    in_gaugep = gaugep),

  GSIMgaugedmo_filenames = read_GSIMgauged_paths(
    inp_GSIMindicesdir = path_GSIMindicesdir,
    in_gaugep = gaugep,
    timestep = 'month'),

  GSIMgaugedsea_filenames = read_GSIMgauged_paths(
    inp_GSIMindicesdir = path_GSIMindicesdir,
    in_gaugep = gaugep,
    timestep = 'season'),

  GRDCqstats = rbindlist(future_map(GRDCgauged_filenames,
                                    comp_GRDCqstats,
                                    maxgap = 20,
                                    minyear = 1971,
                                    maxyear = 2000,
                                    verbose = FALSE,
                                    .progress = TRUE)),

  basemaps = get_basemapswintri()
)

########################### plan_setupdata ####################################
plan_setupdata <- drake_plan(
  GRDCgaugestats = future_map(GRDCgauged_filenames, #To map
                              comp_GRDCdurfreq, #Function to run on each file name
                              maxgap = 20,
                              in_gaugep = gaugep,
                              windowsize = 20,
                              fullwindow = FALSE,
                              monthsel = NULL, #Other arguments in function
                              mdurthresh = 1,
                              verbose = FALSE,
                              .progress = TRUE),

  GSIMgaugestats = future_map2(GSIMgaugedmo_filenames, #To map
                               GSIMgaugedsea_filenames,
                               comp_GSIMdurfreq, #Function to run on each file name
                               maxgap = 20,
                               in_gaugep = gaugep,
                               windowsize = 20,
                               fullwindow = FALSE,
                               monthsel = NULL, #Other arguments in function
                               mdurthresh = 1,
                               .progress = TRUE),

  # GRDCplots = plot_GRDCflags(in_GRDCgaugestats = GRDCgaugestats,
  #                            yearthresh = 1800,
  #                            inp_resdir = path_resdir,
  #                            maxgap = 20),
  # GSIMplots = plot_GSIM(in_GSIMgaugestats = GSIMgaugestats,
  #                       yearthresh = 1800,
  #                       inp_resdir = path_resdir,
  #                       maxgap = 20,
  #                       showmissing = TRUE),


  gaugestats_analyzed = analyzemerge_gaugeir(in_GRDCgaugestats = GRDCgaugestats,
                                             in_GSIMgaugestats = GSIMgaugestats,
                                             in_gaugep = gaugep,
                                             inp_resdir = path_resdir,
                                             plotseries = FALSE),

  gaugestats_format = target(format_gaugestats(in_gaugestats = gaugestats_analyzed$data,
                                               in_gaugep = gaugep,
                                               yearthresh = 1800),
                             trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,
  
  
  watergap_stats = eval_watergap(in_qstats = GRDCqstats,
                                 in_selgauges = gaugestats_format,
                                 binarg = c(1, 10, 100, 1000, 10000, 1000000)
  ),
  
  predvars = target(
    selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta,
                          in_gaugestats = gaugestats_format),
    trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,
  
  
  measures = target(list(classif = msr("classif.bacc"),
                         regr = msr("regr.mae"))
  ),
  
  rivernetwork = rformat_network(in_predvars = predvars,
                                 #in_monthlydischarge = monthlydischarge,
                                 inp_riveratlasmeta = path_riveratlas_meta,
                                 inp_riveratlas = path_riveratlas,
                                 inp_riveratlas2 = path_riveratlas2
  )
  ,

  #-------------------- set-up tasks -------------------------------------
  gauges_div = target(
    gaugestats_format[(dis_m3_pyr >= discharge_interval[1]) &
                        (dis_m3_pyr < discharge_interval[2]),]
    ,
    transform = map(discharge_interval = list(c(0, 10), c(1, Inf)),
                    .names = c('gaugestats_format_u10',
                               'gaugestats_format_o1')
    ),
    trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,

  tasks = target(
    create_tasks(in_gaugestats = in_gauges,
                 in_predvars = predvars,
                 id_suffix = in_id,
                 include_discharge = map_includedis),
    transform = map(
      in_gauges = c(gaugestats_format_u10,
                    gaugestats_format_o1),
      in_id = c('_u10', '_o1'),
      map_includedis = c(TRUE, TRUE),
      .names = c('tasks_u10', 'tasks_o1'),
      tag_in = task,
      tag_out = size
    ),
    trigger  = trigger(mode = "condition", condition =FALSE)
  )
)



########################### plan_runmodels  #########################################
plan_runmodels <- drake_plan(
  baselearners = target(
    create_baselearners(tasks),
    dynamic = map(tasks),
    trigger  = trigger(mode = "condition", condition =FALSE)
    )
  ,

  #Subdivide the baselearners dynamic targets
  seplearners = target(readd(baselearners, subtarget_list = FALSE),
                       trigger  = trigger(mode = "condition", condition =FALSE)),

  autotuningset = target(
    set_tuning(in_learner = seplearners,
               in_measure = measures,
               nfeatures = length(tasks$classif$feature_names),
               insamp_nfolds = 4, insamp_neval = 100,
               insamp_nbatch = parallel::detectCores(logical=FALSE)-2
    ),
    dynamic = map(seplearners),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  resamplingset = target(
    set_cvresampling(rsmp_id = 'repeated_cv',
                     in_task = tasks$classif,
                     outsamp_nrep = 2,
                     outsamp_nfolds = 3),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),
  
  rfresampled_classif = target(
    dynamic_resample(in_task = tasks$classif,
                     in_learner = autotuningset,
                     in_resampling = resamplingset,
                     store_models = FALSE,
                     type = 'classif'),
    dynamic = map(autotuningset),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  rfresampled_regr = target(
    dynamic_resample(in_task = in_tsk,
                     in_learner = in_lrn,
                     in_resampling = resamplingset,
                     store_models = FALSE,
                     type = 'regr'),
    transform = map(in_tsk = c(tasks$regr, tasks$regover),
                    in_lrn = c(autotuningset[[7]], autotuningset[[8]]),
                    .names = c('rfresampled_regr_res', 'rfresampled_regover_res')),
    trigger  = trigger(mode = "condition", condition =FALSE)

  ),

  rfbm_classif = target(
    combine_bm(in_resampleresults = readd(rfresampled_classif,
                                          subtarget_list = TRUE),
               write_qs = T, inp_resdir = path_resdir),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  # bm_checked = target(
  #   analyze_benchmark(in_bm, in_measure = measures),
  #   transform = map(in_bm = list(file_in(!!file.path(resdir,'rfbm_classif.qs')),
  #                                c(rfresampled_regr_res, rfresampled_regover_res)))#,
  #                 #in_measure = list(measures$classif, measures$regr))
  # ),

  selected_learner = target("oversample.classif.ranger"),

  tasks_featsel = target(
    select_features(
      in_bm = rfbm_classif,
      in_lrnid =  selected_learner,
      in_task = tasks$classif,
      pcutoff = 0.05
    ),
    trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,

  resamplingset_featsel = target(
    set_cvresampling(rsmp_id = in_strategy,
                     in_task = tasks$classif,
                     outsamp_nrep = in_outrep,
                     outsamp_nfolds = in_outfolds),
    transform = map(in_strategy = c('repeated_cv', "repeated-spcv-coords"),
                    in_outrep = c(2, 1),
                    in_outfolds = c(3, 40),
                    .names = c('featsel_cv', 'featsel_spcv'))
  ),

  # rfresampled_featsel = target(
  #   dynamic_resamplebm(in_task = in_taskfeatsel,
  #                      in_bm = rfbm_classif,
  #                      in_lrnid =  selected_learner,
  #                      in_resampling = in_resampling,
  #                      store_models = FALSE,
  #                      type = 'classif'),
  #   transform= cross(in_taskfeatsel = c(tasks_featsel[[1]], tasks_featsel[[2]]),
  #                    in_resampling = c(featsel_cv, featsel_spcv),
  #                    .names = c('res_all_cv', 'res_featsel_cv',
  #                               'res_all_spcv', 'res_featsel_spcv')),
  #   trigger  = trigger(mode = "condition", condition =TRUE)
  # ),
  
  res_featsel_spcv = dynamic_resamplebm(in_task = tasks_featsel[[2]],
                                        in_bm = rfbm_classif,
                                        in_lrnid =  selected_learner,
                                        in_resampling = featsel_spcv,
                                        store_models = FALSE,
                                        type = 'classif'),
  
                                        
  rfeval_featall = target(c(res_all_cv, res_all_spcv),
                          trigger  = trigger(mode = "condition", condition =FALSE)),

  rfeval_featsel = target(c(res_featsel_cv, res_featsel_spcv), #Cannot use combine as lead to BenchmarkResult directly in the branching
                          trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,

  rfbm_featall= target(
    analyze_benchmark(in_bm = rfeval_featall,
                      in_measure = measures),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  rfbm_featsel = target(
    analyze_benchmark(in_bm = rfeval_featsel,
                      in_measure = measures),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  interthresh = target(0.5), #target(rfbm_featsel$interthresh_dt[learner == selected_learner, max(thresh)]),

  # Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
  rftuned = target(
    selecttrain_rf(in_rf = res_featsel_cv,
                   in_learnerid = selected_learner,
                   in_task = "inter_class_featsel"),
    trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  # vimp_plot = ggvimp(in_rftuned = rftuned, in_predvars = predvars,
  #                    varnum=20, spatial_rsp = FALSE),
  #
  # pd_plot = ggpartialdep(in_rftuned=rftuned,
  #                        in_predvars=predvars,
  #                        colnums=1:27,
  #                        nvariate=1,  nodupli = FALSE, ngrid = 20, parallel = F,
  #                        spatial_rsp = FALSE),
  #
  table_allbm = target(
    tabulate_benchmarks(in_bm, in_bmid, interthresh=0.5),
    transform = map(
      in_bm = list(
        rfbm_classif,
        c(rfresampled_regr_res, rfresampled_regover_res),
        c(rfeval_featall, rfeval_featsel)),
      in_bmid = list('classif1', 'regr1', 'classif2'),
      .names = c('tablebm_classif1', 'tablebm_regr1', 'tablebm_classif2'))
  ),
  #
  # bin_rftunedmisclass = bin_misclass(in_resampleresult = res_featsel_spcv,
  #                                    binvar = 'dis_m3_pyr',
  #                                    binfunc = 'manual',
  #                                    binarg = c(0.1, 1, 10, 100, 1000, 10000, 1000000),
  #                                    interthresh=interthresh,
  #                                    spatial_rsp=TRUE
  # ),

  misclass_format = target(
    formatmisclass_bm(in_bm = in_bm, in_bmid = in_bmid),
    transform = map(
      in_bm = list(
        rfbm_classif,
        c(rfresampled_regr_res, rfresampled_regover_res),
        c(rfeval_featall, rfeval_featsel)),
      in_bmid = list('classif1', 'regr1', 'classif2'),
      .names = c('misclass_classif1', 'misclass_regr1', 'misclass_classif2')),
    trigger = trigger(mode = "condition", condition =FALSE)
  ),

  misclass_plot = target(
    ggmisclass_bm(list(misclass_classif1,
                       misclass_regr1,
                       misclass_classif2)),
    trigger = trigger(mode = "condition", condition =FALSE)
  )
  ,

  gpredsdt = target(
    make_gaugepreds(in_rftuned = rftuned,
                    in_res_spcv = res_featsel_spcv,
                    in_gaugestats = gaugestats_format,
                    in_predvars = predvars,
                    interthresh = interthresh)
  )
)

########################### plan_runsimplemodels_30d ###########################
plan_runsimplemodels_30d <- drake_plan(
  tasks_featsel_mdur30 = change_tasktarget(in_task = tasks_featsel[[2]],
                                           in_gaugestats = gaugestats_format,
                                           newmdurthresh = 30),

  autotunerdur10 = change_oversamp(in_oversamprf = rftuned$rf_inner,
                                   in_task = tasks_featsel_mdur30),

  rftuned_mdur30 = selecttrain_rf(in_rf = autotunerdur10,
                                  in_task = tasks_featsel_mdur30,
                                  insamp_nfolds =  4, insamp_nevals = 100),

  gpredsdt_mdur30 = make_gaugepreds(in_rftuned = rftuned_mdur30,
                                    in_gaugestats = gaugestats_format,
                                    in_predvars = predvars,
                                    interthresh = interthresh,
                                    simple = TRUE)
)

########################### plan_getpreds #####################################
plan_getpreds <- drake_plan(
  gpredsdt = bind_gaugepreds(list(gpredsdt_u10, gpredsdt_o1),
                             interthresh = 0.5
  )
  ,

  rfpreds_gauges = write_gaugepreds(in_gaugep = gaugep,
                                    in_gpredsdt = gpredsdt,
                                    outp_gaugep = outpath_gaugep
  )
  ,

  rfpreds_network = write_netpreds(
    in_network = rivernetwork,
    in_rftuned = list(rftuned_u10, rftuned_o1),
    predcol = 'predbasic800',
    discharge_interval = list(c(0, 10), c(1, Inf)),
    interthresh = list(interthresh_u10, interthresh_o1),
    in_predvars = predvars,
    outp_riveratlaspred = outpath_riveratlaspred
  ),

  gauges_plot = target(
    gggauges(in_gaugepred =rfpreds_gauges,
             in_basemaps = basemaps,
             binarg = c(30, 60, 100),
             binvar = 'totalYears_kept_o1800'
    ),
    trigger = trigger(mode = "condition", condition =FALSE)
  )
  ,

  envhist = target(
    layout_ggenvhist(in_rivernetwork = rivernetwork,
                     in_gaugepred = rfpreds_gauges,
                     in_predvars = predvars),
    trigger = trigger(mode = "condition", condition =FALSE)
  )
)

########################### plan_getpreds_30d #####################################
plan_getpreds_30d <- drake_plan(
  gpredsdt_mdur30 = bind_gaugepreds(list(gpredsdt_mdur30_u10,
                                         gpredsdt_mdur30_o1),
                                    interthresh = 0.5
  )
  ,

  rfpreds_gauges_mdur30 = write_gaugepreds(
    in_gaugep = gaugep,
    in_gpredsdt = gpredsdt_mdur30,
    outp_gaugep = gsub('[.]gpkg',
                       '_mdur30.gpkg',
                       outpath_gaugep)
  )
  ,

  rfpreds_network_mdur30 = write_netpreds(
    in_network = rivernetwork,
    in_rftuned = list(rftuned_mdur30_u10, rftuned_mdur30_o1),
    predcol = 'predbasic800_mdur30',
    discharge_interval = list(c(0, 10), c(1, Inf)),
    interthresh = list(interthresh_u10, interthresh_o1),
    in_predvars = predvars,
    outp_riveratlaspred = gsub('[.]csv',
                               '_mdur30.csv',
                               outpath_riveratlaspred)
  )
)
########################### plan_getoutputs #####################################
plan_getoutputs <- drake_plan(
  bin_finalmisclass = bin_misclass(in_predictions = gpredsdt,
                                   binvar = 'dis_m3_pyr',
                                   binfunc = 'manual',
                                   binarg = c(0.1, 1, 10, 100, 1000, 10000, 1000000),
                                   interthresh=0.5,
                                   spatial = TRUE
  ),

  gaugeIPR_plot = gggaugeIPR(in_gpredsdt = gpredsdt,
                             in_predvars = predvars,
                             spatial_rsp = TRUE,
                             yearthresh = 1800),

  rivpred = netpredformat(in_rivernetwork = rivernetwork,
                          outp_riveratlaspred = rfpreds_network),

  basinBACC = map_basinBACC(in_gaugepred = gpredsdt,
                            in_rivernetwork = rivernetwork,
                            inp_basin = path_bas03,
                            outp_basinerror = outpath_bas03error,
                            spatial_rsp = TRUE
  )
  

  # NOT SUPER USEFUL - maybe only to reply to reviewers around spatial auto-correlation
  # krigepreds = krige_spgaugeIPR(in_rftuned = rftuned,
  #                               in_gaugep = gaugep,
  #                               in_gaugestats = gaugestats_format,
  #                               inp_bufrasdir = path_bufrasdir,
  #                               kcutoff=50000, overwrite=T),
  #
  # krigepreds_mosaic = mosaic_kriging(in_filestructure = filestructure,
  #                                    in_kpathlist = krigepreds,
  #                                    outp_krigingtif = outpath_krigingtif,
  #                                    overwrite = TRUE),

  # globaltables = target(
  #   tabulate_globalsummary(outp_riveratlaspred = rfpreds_network,
  #                          inp_riveratlas = path_riveratlas,
  #                          inp_riveratlas_legends = path_riveratlas_legends,
  #                          idvars = in_idvars,
  #                          castvar = 'dis_m3_pyr',
  #                          castvar_num = FALSE,
  #                          weightvar = 'LENGTH_KM',
  #                          valuevar = 'predbasic800cat',
  #                          valuevarsub = 1,
  #                          binfunc = 'manual',
  #                          binarg = c(1, 10, 100, 1000, 10000, 1000000),
  #                          na.rm=T,
  #                          tidy = FALSE,
  #                          nolake = TRUE,
  #                          nozerodis = FALSE),
  #   transform = map(in_idvars = c('gad_id_cmj', 'fmh_cl_cmj',
  #                                 'tbi_cl_cmj', 'clz_cl_cmj'))
  # )
)

########################### plan_compareresults ################################
plan_compareresults <- drake_plan(
  fr_plot = compare_fr(inp_frdir = path_frresdir,
                       in_rivpred = rivpred,
                       binarg = c(0.1, 1, 10, 100, 1000, 10000, 100000)
  )
  ,
                       # binarg = c(10,20,50,100,200,500,1000,
                       #            2000,5000,10000,50000,100000,150000)),

  us_plot = compare_us(inp_usresdir = path_usresdir,
                       inp_usdatdir = path_usdatdir,
                       in_rivpred = rivpred,
                       binarg = c(0.1, 1, 10, 100, 1000, 10000, 100000)
                       #c(10,20,50,100,200,500,1000,
                       #        2000,5000,10000,100000,2000000, 3200000)
  )
  ,

  pnw_plot = qc_pnw(inp_pnwresdir = path_pnwresdir,
                    in_rivpred = rivpred,
                    interthresh = 0.5),

  onde_plot = qc_onde(inp_ondedatdir = path_ondedatdir,
                      inp_onderesdir = path_onderesdir,
                      inp_riveratlas = path_riveratlas,
                      in_rivpred = rivpred,
                      interthresh=0.5
  )
)


########################### BIND AND BRANCH PLANS ##############################
plan_runmodels_branches_default <- lapply(
  c('_u10', '_o1'), function(suffix) {
    return(
      branch_plan(
        plan = plan_runmodels,
        branch_suffix = suffix,
        external_arguments_to_modify = c('tasks', 'gaugestats_format'),
        verbose = FALSE)
    )
  }
) %>%
  do.call(bind_plans, .)

plan_runsimplemodels_branches_30d <- lapply(
  c('_u10', '_o1'), function(suffix) {
    return(
      branch_plan(
        plan = plan_runsimplemodels_30d,
        branch_suffix = suffix,
        external_arguments_to_modify = c('tasks_featsel',
                                        'gaugestats_format',
                                        'rftuned', 'interthresh'),
        verbose = FALSE)
    )
  }
) %>%
  do.call(bind_plans, .)

plan_getoutputs_30d <- branch_plan(
  plan = plan_getoutputs,
  branch_suffix = '_mdur30',
  external_arguments_to_modify = c('gpredsdt',
                                   'rfpreds_gauges',
                                   'rfpreds_network'),
  verbose = FALSE) %>%
  .[
    substr(.$target, 1, 12) == 'globaltables',] %>%
  mutate(command = lapply(command, function(call) {
    rlang::call_modify(call, valuevar = "predbasic800_mdur30cat")
  })
  )

plan <- bind_plans(plan_preprocess,
                   plan_setupdata,
                   plan_runmodels_branches_default,
                   #plan_runsimplemodels_branches_30d,
                   plan_getpreds,
                   #plan_getpreds_30d,
                   plan_getoutputs,
                   #plan_getoutputs_30d,
                   #plan_compareresults
)

