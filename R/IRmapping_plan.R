source('R/IRmapping_packages.R')
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
  path_linkpop = file_in(!!file.path(rootdir, 'results\\worldpop_link_20201220.csv')),
  path_frresdir = file.path(path_compresdir, 'france.gdb'),
  path_usdatdir = file.path(rootdir, 'data\\Comparison_databases', 'US'),
  path_usresdir = file.path(path_compresdir, 'us.gdb'),
  path_auresdir = file.path(path_compresdir, 'australia.gdb'),
  path_insitudatdir = file.path(rootdir, 'data\\Insitu_databases'),
  path_insituresdir = file.path(path_resdir, 'Insitu_databases'),
  path_pnwresdir =  file.path(path_insituresdir, 'pnw.gdb'),
  path_ondedatdir = file.path(path_insitudatdir, 'OndeEau'),
  path_onderesdir = file.path(path_insituresdir, 'ondeeau.gdb'),

  outpath_riveratlaspred = file.path(rootdir, 'results\\RiverATLAS_predbasic800.csv'),
  outpath_gaugep = file_out(!!file.path(rootdir, 'results\\GRDCstations_predbasic800.gpkg')),
  outpath_bas03error = file_out(!!file.path(rootdir, 'results\\BasinATLAS_v10_lev03_errors.gpkg')),
  path_bufrasdir = file.path(path_resdir, 'bufrasdir'),

  #outpath_krigingtif = file_out("E:\\Mathis\\results\\prederror_krigingtest.tif"),

  #-------------------- Pre-analysis ------------------------------------------
  #monthlydischarge = read_monthlydis(in_path = path_monthlynetdischarge),

  gaugep = target(
    read_gaugep(inp_GRDCgaugep = path_GRDCgaugep,
                inp_GSIMgaugep = path_GSIMgaugep,
                #in_monthlydischarge = monthlydischarge,
                inp_riveratlas2 = path_riveratlas2
    )
    ,
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
                                               yearthresh = 1800)
                             # ,
                             # trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,


  watergap_stats = eval_watergap(in_qstats = GRDCqstats,
                                 in_selgauges = gaugestats_format,
                                 binarg = c(1, 10, 100, 1000, 10000, 1000000)
  ),

  predvars = target(
    selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta,
                          in_gaugestats = gaugestats_format)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
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

  netlength_extra = extrapolate_networklength(
    inp_riveratlas = path_riveratlas,
    min_cutoff = 0.1,
    dispred =  seq(0.01, 0.09, 0.01),
    interactive = F,
    grouping_var = 'PFAF_IDclz'
  ),

  #-------------------- set-up tasks -------------------------------------
  gauges_div = target(
    gaugestats_format[(dis_m3_pyr >= discharge_interval[1]) &
                        (dis_m3_pyr < discharge_interval[2]),]
    ,
    transform = map(discharge_interval = list(c(0, Inf), c(0, 10), c(1, Inf)),
                    .names = c('gaugestats_format_all',
                               'gaugestats_format_u10',
                               'gaugestats_format_o1')
    )
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  )
  ,

  tasks = target(
    create_tasks(in_gaugestats = in_gauges,
                 in_predvars = predvars,
                 id_suffix = in_id,
                 include_discharge = map_includedis),
    transform = map(
      in_gauges = c(gaugestats_format_all,
                    gaugestats_format_u10,
                    gaugestats_format_o1),
      in_id = c('_all', '_u10', '_o1'),
      map_includedis = c(TRUE, TRUE, TRUE),
      .names = c('tasks_all', 'tasks_u10', 'tasks_o1')
    )
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  )
)



########################### plan_runmodels  #########################################
plan_runmodels <- drake_plan(
  baselearners = target(
    create_baselearners(tasks),
    dynamic = map(tasks)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
    )
  ,

  #Subdivide the baselearners dynamic targets
  seplearners = target(readd(baselearners, subtarget_list = FALSE)
                       # ,
                       # trigger  = trigger(mode = "condition", condition =FALSE)
                       ),

  ###when update mlr3, uncomment store_tuning_instance
  autotuningset = target(
    set_tuning(in_learner = seplearners,
               in_measure = measures,
               nfeatures = length(tasks$classif$feature_names),
               insamp_nfolds = 4, insamp_neval = 100,
               insamp_nbatch = parallel::detectCores(logical=FALSE)-1
    ),
    dynamic = map(seplearners)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  resamplingset = target(
    set_cvresampling(rsmp_id = 'repeated_cv',
                     in_task = tasks$classif,
                     outsamp_nrep = 2,
                     outsamp_nfolds = 3)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  rfresampled_classif = target(
    dynamic_resample(in_task = tasks$classif,
                     in_learner = autotuningset,
                     in_resampling = resamplingset,
                     store_models = TRUE,
                     type = 'classif'),
    dynamic = map(autotuningset)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  rfresampled_regr = target(
    dynamic_resample(in_task = in_tsk,
                     in_learner = in_lrn,
                     in_resampling = resamplingset,
                     store_models = FALSE,
                     type = 'regr'),
    transform = map(in_tsk = c(tasks$regr, tasks$regover),
                    in_lrn = c(autotuningset[[7]], autotuningset[[8]]),
                    .names = c('rfresampled_regr_res', 'rfresampled_regover_res'))
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)

  ),

  rfbm_classif = target(
    combine_bm(in_resampleresults = readd(rfresampled_classif,
                                          subtarget_list = TRUE),
               write_qs = T, inp_resdir = path_resdir)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  selected_learner = target("oversample.classif.ranger"),

  tasks_featsel = target(
    select_features(
      in_bm = rfbm_classif,
      in_lrnid =  selected_learner,
      in_task = tasks$classif,
      pcutoff = 0.05,
      inp_resdir = path_resdir
    )
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
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
    # ,
    # trigger  = trigger(mode = "condition", condition =TRUE)
  ),

  rfresampled_featsel = target(
    dynamic_resamplebm(in_task = in_taskfeatsel,
                       in_bm = rfbm_classif,
                       in_lrnid =  selected_learner,
                       in_resampling = in_resampling,
                       store_models = store_models,
                       inp_resdir = path_resdir,
                       type = 'classif'),
    transform= map(in_taskfeatsel = c(tasks_featsel[[1]], tasks_featsel[[2]],
                                      tasks_featsel[[1]], tasks_featsel[[2]]),
                   in_resampling = c(featsel_cv, featsel_cv,
                                     featsel_spcv, featsel_spcv),
                   store_models = c(FALSE, TRUE, FALSE, FALSE),
                   .names = c('res_all_cv', 'res_featsel_cv',
                              'res_all_spcv', 'res_featsel_spcv'))
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  rfeval_featall = target(c(res_all_cv))
  ,

  rfeval_featsel = target(c(res_featsel_cv, res_featsel_spcv), #Cannot use combine as lead to BenchmarkResult directly in the branching
  )
  ,

  rfbm_featall= target(
    analyze_benchmark(in_bm = rfeval_featall,
                      in_measure = measures)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  rfbm_featsel = target(
    analyze_benchmark(in_bm = rfeval_featsel,
                      in_measure = measures)
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),

  interthresh = target(0.5), #target(rfbm_featsel$interthresh_dt[learner == selected_learner, max(thresh)]),

  # Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
  rftuned = target(
    selecttrain_rf(in_rf = res_featsel_cv,
                   in_learnerid = selected_learner,
                   in_task = "inter_class_featsel")
    # ,
    # trigger  = trigger(mode = "condition", condition =FALSE)
  ),


  vimp_plot = ggvimp(in_rftuned = rftuned, in_predvars = predvars,
                     varnum=20, spatial_rsp = FALSE),

  pd_plot = ggpartialdep(in_rftuned=rftuned,
                         in_predvars=predvars,
                         colnums=1:27,
                         nvariate=1,  nodupli = FALSE, ngrid = 20, parallel = F,
                         spatial_rsp = FALSE),

  table_allbm = target(
    tabulate_benchmarks(in_bm, in_bmid, inp_resdir = path_resdir, interthresh=0.5),
    transform = map(
      in_bm = list(
        rfbm_classif,
        c(rfresampled_regr_res, rfresampled_regover_res),
        c(rfeval_featall, rfeval_featsel)),
      in_bmid = list('classif1', 'regr1', 'classif2'),
      .names = c('tablebm_classif1', 'tablebm_regr1', 'tablebm_classif2'))
  ),

  bin_rftunedmisclass = bin_misclass(in_resampleresult = res_featsel_spcv,
                                     binvar = 'dis_m3_pyr',
                                     binfunc = 'manual',
                                     binarg = c(0.1, 1, 10, 100, 1000, 10000, 1000000),
                                     interthresh=interthresh,
                                     spatial_rsp=TRUE
  ),

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
                                  insamp_nfolds =  4, insamp_nevals = 25),

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

  rfpreds_network = target(
    write_netpreds(
      in_network = rivernetwork,
      in_rftuned = list(rftuned_u10, rftuned_o1),
      predcol = 'predbasic800',
      discharge_interval = list(c(0, 10), c(1, Inf)),
      interthresh = list(interthresh_u10, interthresh_o1),
      in_predvars = predvars,
      outp_riveratlaspred = outpath_riveratlaspred
    )
   #  ,
   # trigger = trigger(mode = "condition", condition =TRUE)
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

  rfpreds_network_mdur30 = target(
    write_netpreds(
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
    # ,
    # trigger = trigger(mode = "condition", condition =FALSE)
  )
)
########################### plan_getoutputs #####################################
plan_getoutputs_1 <- drake_plan(
  gauges_plot = target(
    gggauges(in_gaugepred =rfpreds_gauges,
             in_basemaps = basemaps,
             binarg = c(30, 60, 100),
             binvar = 'totalYears_kept_o1800'
    )
    # ,
    # trigger = trigger(mode = "condition", condition =FALSE)
  )
  ,

  envhist = target(
    layout_ggenvhist(in_rivernetwork = rivernetwork,
                     in_gaugepred = rfpreds_gauges,
                     in_predvars = predvars)
    # ,
    # trigger = trigger(mode = "condition", condition =FALSE)
  ),

  bin_finalmisclass = target(
    bin_misclass(in_predictions = gpredsdt,
                 binvar = 'dis_m3_pyr',
                 binfunc = 'manual',
                 binarg = c(0.1, 1, 10, 100, 1000, 10000, 1000000),
                 interthresh=0.5,
                 rspcol = rspcol
    ),
    map(rspcol = c('IRpredcat_CVnosp', 'IRpredcat_CVsp')
    )
  )
  ,

  threshold_sensitivity = test_thresholdsensitivity(
    in_gpredsdt = gpredsdt,
    in_rivpred = rivpred,
    threshrange_gauges = seq(0.25, 0.75, 0.01),
    threshrange_network = seq(0.45, 0.55, 0.01),
    mincutoff = 0.1,
    gaugescol = 'IRpredprob_CVnosp',
    netcol = 'predbasic800'
  ),

  gaugeIPR_plot = gggaugeIPR(in_gpredsdt = gpredsdt,
                             in_predvars = predvars,
                             spatial_rsp = TRUE,
                             yearthresh = 1800)
  ,

  
  rivpred = target(
    netpredformat(in_rivernetwork = rivernetwork,
                  outp_riveratlaspred = rfpreds_network)
    # ,
    # trigger = trigger(mode = "condition", condition =FALSE)
  ),
  
  basinBACC = target(
    map_basinBACC(in_gaugepred = gpredsdt,
                  inp_basin = path_bas03,
                  in_rivernetwork = rivpred,
                  outp_basinerror = outpath_bas03error,
                  spatial_rsp = TRUE)
    # ,
    # trigger = trigger(mode = "condition", condition =FALSE)
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
)


plan_getoutputs_2 <- drake_plan(
  IRESextra = extrapolate_IRES(in_rivpred = rivpred,
                               in_extranet = netlength_extra,
                               min_cutoff = 0.1,
                               valuevar = 'predbasic800cat',
                               interactive = F),
  
  globaltables = target(
    tabulate_globalsummary(outp_riveratlaspred = rfpreds_network,
                           inp_riveratlas = path_riveratlas,
                           inp_riveratlas_legends = path_riveratlas_legends,
                           interthresh = 0.5,
                           idvars = in_idvars,
                           castvar = 'dis_m3_pyr',
                           castvar_num = FALSE,
                           weightvar = 'LENGTH_KM',
                           valuevar = 'predbasic800',
                           valuevarsub = 1,
                           binfunc = 'manual',
                           binarg = c(0.1, 1, 10, 100, 1000, 10000, 1000000),
                           na.rm=T,
                           tidy = FALSE,
                           nolake = TRUE,
                           mincutoff = 0.1),
    transform = map(in_idvars = c('gad_id_cmj', 'fmh_cl_cmj',
                                  'tbi_cl_cmj', 'clz_cl_cmj'))
    # ,
    # trigger = trigger(mode = "condition", condition =TRUE)
  ),

  globaltable_clzextend = extend_globalsummary_clz(
    in_IRESextra = IRESextra,
    in_globaltable = globaltables_clz_cl_cmj,
    inp_riveratlas_legends = path_riveratlas_legends
  ),


  IRpop = compute_IRpop(in_rivpred = rivpred,
                        inp_linkpop = path_linkpop,
                        valuevar = 'predbasic800cat'
  )
)

########################### plan_compareresults ################################
plan_compareresults <- drake_plan(
  fr_plot = compare_fr(inp_frdir = path_frresdir,
                       in_rivpred = rivpred,
                       predcol = 'predbasic800cat',
                       mincutoff = 0.1,
                       binarg = c(0.1, 1, 10, 100, 1000, 10000, 100000)
  )
  ,

  us_plot = compare_us(inp_usresdir = path_usresdir,
                       inp_usdatdir = path_usdatdir,
                       in_rivpred = rivpred,
                       predcol = 'predbasic800cat',
                       binarg = c(0.1, 1, 10, 100, 1000, 10000, 100000),
                       mincutoff = 0.1
  ),

  au_plot = compare_au(inp_resdir = path_resdir,
                       in_rivpred = rivpred,
                       predcol = 'predbasic800cat',
                       binarg = c(100, 10^3, 10^4, 10^5, 10^6, 10^7)
                       ),

  pnw_plot = qc_pnw(inp_pnwresdir = path_pnwresdir,
                    in_rivpred = rivpred,
                    predcol = 'predbasic800cat',
                    interthresh = 0.5,
                    mincutoff = 0.1),

  onde_plot = qc_onde(inp_ondedatdir = path_ondedatdir,
                      inp_onderesdir = path_onderesdir,
                      inp_riveratlas = path_riveratlas,
                      in_rivpred = rivpred,
                      predcol = 'predbasic800cat',
                      interthresh=0.5,
                      mincutoff = 0.1
  )
)


########################### BIND AND BRANCH PLANS ##############################
plan_runmodels_branches_default <- lapply(
  c('_all', '_u10', '_o1'), function(suffix) {
    return(
      branch_plan(
        plan = plan_runmodels,
        branch_suffix = suffix,
        external_arguments_to_modify = c('tasks', 'gaugestats_format'),
        verbose = FALSE)
    )
  }
) %>%
  do.call(bind_plans, .) %>%
  .[!(.$target %in% c('res_all_spcv_all', 'res_featsel_spcv_all',
                      'rfeval_featall_all', 'rfeval_featsel_all',
                      'rfbm_featall_all', 'rfbm_featsel_all',
                      'bin_rftunedmisclass_all', 'res_all_spcv_o1',
                      'res_all_spcv_u10', 'rftuned_all', 'gpredsdt_all',
                      'vimp_plot_all', 'pd_plot_all',
                      'tablebm_classif1_all', 'tablebm_classif2_all', 'tablebm_regr1_all')
      ),]

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
  plan = plan_getoutputs_2,
  branch_suffix = '_mdur30',
  external_arguments_to_modify = c('gpredsdt',
                                   'rfpreds_gauges',
                                   'rfpreds_network'),
  verbose = FALSE) %>%
  as.data.table(.) %>%
  .[substr(target, 1, 12) == 'globaltables',
    command :=
      lapply(command, function(call) {
        rlang::call_modify(.call=call,
                           valuevar = "predbasic800_mdur30")
      })] %>%
  .[target %in% c('IRESextra_mdur30', 'IRpop_mdur30'),
    command := lapply(command, function(call) {
      rlang::call_modify(.call=call,
                         valuevar = "predbasic800_mdur30cat")
    })] %>%
  as_tibble

plan_compareresults_30d <- branch_plan(
  plan = plan_compareresults,
  branch_suffix = '_mdur30',
  external_arguments_to_modify = 'rivpred',
  verbose = FALSE
) %>%
  as.data.table(.) %>%
  .[,
    command :=
      lapply(command, function(call) {
        rlang::call_modify(.call=call,
                           predcol = "predbasic800_mdur30cat")
      })] %>%
  as_tibble


plan <- bind_plans(plan_preprocess,
                   plan_setupdata,
                   plan_runmodels_branches_default,
                   plan_runsimplemodels_branches_30d,
                   plan_getpreds,
                   plan_getpreds_30d,
                   plan_getoutputs_1,
                   plan_getoutputs_2,
                   plan_getoutputs_30d,
                   plan_compareresults,
                   plan_compareresults_30d
)
