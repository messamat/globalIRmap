
plan <- drake_plan(
  ##################### DEFINE INPUT AND OUTPUT FILES ##########################
  #Note: for dependency tracking, drake only reads string literals inside file_in and file_out; hence the use of absolute paths
  #Tidy evaluation (!!) didn't work for some reason -- can look into it later
  path_resdir = "C:\\globalIRmap\\results",
  path_GRDCgaugep = "C:\\globalIRmap\\results\\spatialoutputs.gdb\\grdcstations_cleanjoin",
  path_GRDCgaugedir =  file_in("C:\\globalIRmap\\data\\GRDCdat_day"),
  path_GSIMgaugep = "C:\\globalIRmap\\results\\GSIM\\GSIM.gdb\\GSIMstations_cleanriverjoin",
  path_GSIMindicesdir =  file_in("C:\\globalIRmap\\data\\GSIM\\GSIM_indices"),
  path_riveratlas_meta = file_in('C:\\globalIRmap\\data\\HydroATLAS\\HydroATLAS_metadata_MLMv11.xlsx'),
  path_riveratlas_legends = file_in('C:\\globalIRmap\\data\\HydroATLAS\\HydroATLAS_v10_Legends.xlsx'),
  path_monthlynetdischarge = 'C:\\globalIRmap\\data\\HydroSHEDS\\HS_discharge_monthly.gdb\\Hydrosheds_discharge_monthly',
  path_riveratlas = file_in('C:\\globalIRmap\\results\\RiverATLAS_v10tab.csv'),
  path_riveratlas2 = file_in('C:\\globalIRmap\\results\\RiverATLAS_v11tab.csv'),
  path_compresdir = file.path('C:\\globalIRmap\\results\\Comparison_databases'),
  path_frresdir = file.path(path_compresdir, 'france.gdb'),
  path_usdatdir = file.path('C:\\globalIRmap\\data\\Comparison_databases', 'US'),
  path_usresdir = file.path(path_compresdir, 'us.gdb'),
  path_insitudatdir = file.path('C:\\globalIRmap\\data\\Insitu_databases'),
  path_insituresdir = file.path('C:\\globalIRmap\\results\\Insitu_databases'),
  path_pnwresdir =  file.path(path_insituresdir, 'pnw.gdb'),
  path_ondedatdir = file.path(path_insitudatdir, 'OndeEau'),
  path_onderesdir = file.path(path_insituresdir, 'ondeeau.gdb'),

  outpath_rfclassif = 'C:\\globalIRmap\\results\\rfbm_classif.qs',
  outpath_gaugep = file_out('C:\\globalIRmap\\results\\GRDCstations_predbasic800.gpkg'),
  outpath_riveratlaspred = file_out('C:\\globalIRmap\\results\\RiverATLAS_predbasic800.csv'),
  path_bufrasdir = file.path('C:\\globalIRmap\\results\\bufrasdir'),
  outpath_krigingtif = file_out("C:\\globalIRmap\\results\\prederror_krigingtest.tif"),


  ##################### ANALYSIS ############################################
  monthlydischarge = read_monthlydis(in_path = path_monthlynetdischarge),


  gaugep = read_gaugep(inp_GRDCgaugep = path_GRDCgaugep,
                       inp_GSIMgaugep = path_GSIMgaugep,
                       inp_riveratlas2 = path_riveratlas2,
                       in_monthlydischarge = monthlydischarge),

  GRDCgauged_filenames = read_GRDCgauged_paths(
    inp_GRDCgaugedir = path_GRDCgaugedir,
    in_gaugep = gaugep),

  GSIMgaugedmo_filenames = read_GSIMgauged_paths(
    inp_GSIMindicesdir = path_GSIMindicesdir,
    in_gaugep = gaugep,
    timestep = 'month'),

  GSIMgaugedse_filenames = read_GSIMgauged_paths(
    inp_GSIMindicesdir = path_GSIMindicesdir,
    in_gaugep = gaugep,
    timestep = 'season'),

  GRDCgaugestats = future_map(GRDCgauged_filenames, #To map
                              comp_GRDCdurfreq, #Function to run on each file name
                              maxgap = 20,
                              in_gaugep = gaugep,
                              monthsel = NULL, #Other arguments in function
                              mdurthresh = 1,
                              .progress = TRUE)
  #
  # gaugestats_format = format_gaugestats(in_gaugestats = gaugestats,
  #                                       in_gaugep = gaugep),
  #
  # predvars = selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta,
  #                                  in_gaugestats = gaugestats_format),
  #
  # tasks = create_tasks(in_gaugestats = gaugestats_format,
  #                      in_predvars = predvars),
  #
  # baselearners = target(
  #   create_baselearners(tasks),
  #   dynamic = map(tasks)),
  #
  # #Subdivide the baselearners dynamic targets from 3 to 8
  # seplearners = target(readd(baselearners, subtarget_list = FALSE)),
  #
  # measures = target(list(classif = msr("classif.bacc"),
  #                        regr = msr("regr.mae"))
  # ),
  #
  # autotuningset = target(
  #   set_tuning(in_learner = seplearners,
  #              in_measures = measures,
  #              nfeatures = length(tasks$classif$feature_names),
  #              insamp_nfolds = 4, insamp_neval = 100,
  #              insamp_nbatch = parallel::detectCores(logical=FALSE)-1
  #   ),
  #   dynamic = map(seplearners)
  # ),
  #
  # resamplingset = set_cvresampling(rsmp_id = 'repeated_cv',
  #                                  in_task = tasks$classif,
  #                                  outsamp_nrep = 2,
  #                                  outsamp_nfolds = 3),
  #
  # rfresampled_classif = target(
  #   dynamic_resample(in_task = tasks$classif,
  #                    in_learner = autotuningset,
  #                    in_resampling = resamplingset,
  #                    store_models = TRUE,
  #                    type = 'classif'),
  #   dynamic = map(autotuningset)
  # ),
  #
  # rfresampled_regr = target(
  #   dynamic_resample(in_task = in_tsk,
  #                    in_learner = in_lrn,
  #                    in_resampling = resamplingset,
  #                    store_models = TRUE,
  #                    type = 'regr'),
  #   transform = map(in_tsk = c(tasks$regr, tasks$regover),
  #                   in_lrn = c(autotuningset[[7]], autotuningset[[8]]),
  #                   .names = c('rfresampled_regr_res', 'rfresampled_regover_res'))
  #
  # ),
  #
  # rfbm_classif = target(
  #   combine_bm(readd(rfresampled_classif, subtarget_list = TRUE),
  #              out_qs = outpath_rfclassif)
  # ),
  #
  # # bm_checked = target(
  # #   analyze_benchmark(in_bm, in_measure = measures),
  # #   transform = map(in_bm = list(file_in(!!file.path(resdir,'rfbm_classif.qs')),
  # #                                c(rfresampled_regr_res, rfresampled_regover_res)))#,
  # #                 #in_measure = list(measures$classif, measures$regr))
  # # ),
  #
  # selected_learner = target("oversample.classif.ranger"),
  #
  # resamplingset_featsel = target(
  #   set_cvresampling(rsmp_id = in_strategy,
  #                    in_task = tasks$classif,
  #                    outsamp_nrep = in_outrep,
  #                    outsamp_nfolds = in_outfolds),
  #   transform = map(in_strategy = c('repeated_cv', "repeated-spcv-coords"),
  #                   in_outrep = c(2, 2),
  #                   in_outfolds = c(3, 3),
  #                   .names = c('featsel_cv', 'featsel_spcv'))
  # ),
  #
  # tasks_featsel = select_features(
  #   in_bm = rfbm_classif,
  #   in_lrnid =  selected_learner,
  #   in_task = tasks$classif,
  #   pcutoff = 0.05
  # ),
  #
  # rfresampled_featsel = target(
  #   dynamic_resamplebm(in_task = in_taskfeatsel,
  #                      in_bm = rfbm_classif,
  #                      in_lrnid =  selected_learner,
  #                      in_resampling = in_resampling,
  #                      store_models = TRUE,
  #                      type = 'classif'),
  #   transform= cross(in_taskfeatsel = c(tasks_featsel[[1]], tasks_featsel[[2]]),
  #                    in_resampling = c(featsel_cv, featsel_spcv),
  #                    .names = c('res_all_cv', 'res_all_spcv',
  #                               'res_featsel_cv', 'res_featsel_spcv'))
  # ),
  #
  # rfeval_featall = target(c(res_all_cv, res_all_spcv)),
  #
  # rfeval_featsel = target(c(res_featsel_cv, res_featsel_spcv) #Cannot use combine as lead to BenchmarkResult directly in the branching
  # ),
  #
  # rfbm_featall= analyze_benchmark(in_bm = rfeval_featall,
  #                                 in_measure = measures$classif),
  #
  # rfbm_featsel = analyze_benchmark(in_bm = rfeval_featsel,
  #                                  in_measure = measures$classif),
  #
  # interthresh = target(rfbm_featsel$interthresh_dt[learner == selected_learner,
  #                                                  max(thresh)]),
  #
  # #  Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
  # rftuned = target(
  #   selecttrain_rf(in_rf = res_featsel_spcv,
  #                  in_learnerid = selected_learner,
  #                  in_taskid = "inter_basicsp_featsel",
  #                  insamp_nfolds = 4,
  #                  insamp_nevals = 100)),
  #
  # rivernetwork = rformat_network(in_predvars = predvars,
  #                                in_monthlydischarge = monthlydischarge,
  #                                inp_riveratlasmeta = path_riveratlas_meta,
  #                                inp_riveratlas = path_riveratlas,
  #                                inp_riveratlas2 = path_riveratlas2),
  #
  # gaugeIPR_plot = gggaugeIPR(in_rftuned = rftuned$rf_outer,
  #                            in_gaugestats = gaugestats_format,
  #                            in_predvars = predvars,
  #                            spatial_rsp = FALSE,
  #                            interthresh = interthresh),
  #
  # rfpreds = write_preds(in_gaugep = gaugep,
  #                       in_gaugestats = gaugestats_format,
  #                       in_network = rivernetwork,
  #                       in_rftuned = rftuned,
  #                       in_predvars = predvars,
  #                       in_gaugeIPR = gaugeIPR_plot[['out_gaugeIPR']],
  #                       interthresh = interthresh,
  #                       outp_gaugep = outpath_gaugep,
  #                       outp_riveratlaspred = outpath_riveratlaspred),
  #
  # rivpred = netpredformat(in_rivernetwork = rivernetwork,
  #                         outp_riveratlaspred = rfpreds[["rivpredpath"]]),
  #
  # vimp_plot = ggvimp(in_rftuned = rftuned, in_predvars = predvars,
  #                    varnum=20, spatial_rsp = FALSE),
  #
  # pd_plot = ggpd_bivariate(in_rftuned=rftuned, in_predvars=predvars, colnums=1:20,
  #                          nvariate=1,  nodupli = TRUE, ngrid = 20, parallel = T,
  #                          spatial_rsp = FALSE),
  #
  # basemaps = get_basemapswintri(),
  #
  # gauges_plot = gggauges(in_gaugepred = rfpreds[["out_gaugep"]],
  #                        in_basemaps = basemaps,
  #                        binarg <- c(30, 60, 100),
  #                        binvar <- 'totalYears_kept'),
  #
  # envhist = layout_ggenvhist(in_rivernetwork = rivernetwork,
  #                            in_gaugepred =  rfpreds[["out_gaugep"]],
  #                            in_predvars = predvars),
  #
  # table_allbm = target(
  #   tabulate_benchmarks(in_bm, in_bmid),
  #   transform = map(
  #     in_bm = list(
  #       rfbm_classif,
  #       c(rfresampled_regr_res, rfresampled_regover_res),
  #       c(rfeval_featall, rfeval_featsel)),
  #     in_bmid = list('classif1', 'regr1', 'classif2'),
  #     .names = c('tablebm_classif1', 'tablebm_regr1', 'tablebm_classif2'))
  # ),
  #
  # misclass_format = target(
  #   formatmisclass_bm(in_bm = in_bm, in_bmid = in_bmid),
  #   transform = map(
  #     in_bm = list(
  #       rfbm_classif,
  #       c(rfresampled_regr_res, rfresampled_regover_res),
  #       c(rfeval_featall, rfeval_featsel)),
  #     in_bmid = list('classif1', 'regr1', 'classif2'),
  #     .names = c('misclass_classif1', 'misclass_regr1', 'misclass_classif2'))
  # ),
  #
  # misclass_plot = ggmisclass_bm(list(misclass_classif1,
  #                                    misclass_regr1,
  #                                    misclass_classif2)),
  #
  # # krigepreds = krige_spgaugeIPR(in_rftuned = rftuned,
  # #                               in_gaugep = gaugep,
  # #                               in_gaugestats = gaugestats_format,
  # #                               inp_bufrasdir = path_bufrasdir,
  # #                               kcutoff=50000, overwrite=T),
  # #
  # # krigepreds_mosaic = mosaic_kriging(in_filestructure = filestructure,
  # #                                    in_kpathlist = krigepreds,
  # #                                    outp_krigingtif = outpath_krigingtif,
  # #                                    overwrite = TRUE),
  #
  # globaltables = target(
  #   tabulate_globalsummary(outp_riveratlaspred = rfpreds[["rivpredpath"]],
  #                          inp_riveratlas = path_riveratlas,
  #                          inp_riveratlas_legends = path_riveratlas_legends,
  #                          idvars = in_idvars,
  #                          castvar = 'dis_m3_pyr',
  #                          castvar_num = FALSE,
  #                          weightvar = 'LENGTH_KM',
  #                          valuevar = 'predbasic800cat',
  #                          valuevarsub = 1,
  #                          binfunc = 'manual',
  #                          binarg = c(0.1, 1, 10, 100, 100, 10000, 100000),
  #                          na.rm=T,
  #                          tidy = FALSE),
  #   transform = map(in_idvars = c('gad_id_cmj', 'fmh_cl_cmj',
  #                                 'tbi_cl_cmj', 'clz_cl_cmj'))
  # )
  # ,
  #
  # fr_plot = compare_fr(inp_frdir = path_frresdir,
  #                      in_rivpred = rivpred,
  #                      binarg = c(10,20,50,100,200,500,1000,
  #                                 2000,5000,10000,50000,100000,150000)),
  #
  # us_plot = compare_us(inp_usresdir = path_usresdir,
  #                      inp_usdatdir = path_usdatdir,
  #                      in_rivpred = rivpred,
  #                      binarg = c(10,20,50,100,200,500,1000,
  #                                 2000,5000,10000,50000,100000,2000000, 3200000))
  # ,
  #
  # pnw_plot = qc_pnw(inp_pnwresdir = path_pnwresdir,
  #                   in_rivpred = rivpred,
  #                   interthresh = interthresh)
################## END OF PLAN  #################################################################
)

