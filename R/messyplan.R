################################## PREPROCESS_PLAN #############################
preprocess_plan <- drake_plan(
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

  gauges_div = target(
    gaugestats_format[(dis_m3_pyr >= discharge_interval[1]) &
                        (dis_m3_pyr < discharge_interval[2]),]
    ,
    transform = map(discharge_interval = list(c(0, 10), c(10, Inf)),
                    .names = c('gauges_u10', 'gauges_o10')
    )
  ),

  tasks = target(
    create_tasks(in_gaugestats = in_gauges,
                 in_predvars = predvars,
                 id_suffix = in_id,
                 include_discharge = map_includedis),
    transform = map(
      in_gauges = c(gauges_u10, gauges_o10),
      in_id = '_u10', '_o10',
      map_includedis = c(FALSE, TRUE),
      .names = c('task_classifu10', 'task_classifo10'),
      tag_in = task,
      tag_out = size
    )
  ),

  measures = target(list(classif = msr("classif.bacc"),
                         regr = msr("regr.mae"))
  ),
)


################################## RUNMODEL_PLAN_BASIC #########################

runmodel_plan_basic <- drake_plan(
  baselearners = target(
    create_baselearners(tasks),
    dynamic = map(tasks)),

  #Subdivide the baselearners dynamic targets from 3 to 8
  seplearners = target(readd(baselearners, subtarget_list = FALSE)),

  autotuningset = target(
    set_tuning(in_learner = seplearners,
               in_measures = measures,
               nfeatures = length(tasks$classif$feature_names),
               insamp_nfolds = 4, insamp_neval = 100,
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
               out_qs = outpath_rfclassif)
  ),

  # bm_checked = target(
  #   analyze_benchmark(in_bm, in_measure = measures),
  #   transform = map(in_bm = list(file_in(!!file.path(resdir,'rfbm_classif.qs')),
  #                                c(rfresampled_regr_res, rfresampled_regover_res)))#,
  #                 #in_measure = list(measures$classif, measures$regr))
  # ),

  selected_learner = target("oversample.classif.ranger"),

  resamplingset_featsel = target(
    set_cvresampling(rsmp_id = in_strategy,
                     in_task = tasks$classif,
                     outsamp_nrep = in_outrep,
                     outsamp_nfolds = in_outfolds),
    transform = map(in_strategy = c('repeated_cv', "repeated-spcv-coords"),
                    in_outrep = c(2, 2),
                    in_outfolds = c(3, 3),
                    .names = c('featsel_cv', 'featsel_spcv'))
  ),

  tasks_featsel = select_features(
    in_bm = rfbm_classif,
    in_lrnid =  selected_learner,
    in_task = tasks$classif,
    pcutoff = 0.05
  ),

  rfresampled_featsel = target(
    dynamic_resamplebm(in_task = in_taskfeatsel,
                       in_bm = rfbm_classif,
                       in_lrnid =  selected_learner,
                       in_resampling = in_resampling,
                       store_models = TRUE,
                       type = 'classif'),
    transform= cross(in_taskfeatsel = c(tasks_featsel[[1]], tasks_featsel[[2]]),
                     in_resampling = c(featsel_cv, featsel_spcv),
                     .names = c('res_all_cv', 'res_featsel_cv',
                                'res_all_spcv', 'res_featsel_spcv'))
  ),

  rfeval_featall = target(c(res_all_cv, res_all_spcv)),

  rfeval_featsel = target(c(res_featsel_cv, res_featsel_spcv) #Cannot use combine as lead to BenchmarkResult directly in the branching
  ),

  rfbm_featall= analyze_benchmark(in_bm = rfeval_featall,
                                  in_measure = measures$classif),

  rfbm_featsel = analyze_benchmark(in_bm = rfeval_featsel,
                                   in_measure = measures$classif),

  interthresh = target(rfbm_featsel$interthresh_dt[learner == selected_learner,
                                                   max(thresh)]),

  #  Assertion on 'uhash' failed: Must be element of set {'f00f1b58-0316-4828-814f-f30310b47761','1b8bb7dc-69a0-49a2-af2e-f377fb162a5a'}, but is not atomic scalar.
  rftuned = target(
    selecttrain_rf(in_rf = res_featsel_spcv,
                   in_learnerid = selected_learner,
                   in_taskid = "inter_basicsp_featsel",
                   insamp_nfolds = 4,
                   insamp_nevals = 100))
)
################################## AMEND PLANS #########################
library(rlang)
runmodel_plan_u10 <- runmodel_plan_basic
runmodel_plan_u10$target <- paste0(runmodel_plan_u10$target, '_u10')


target_suffix <- '_u10'
target_regex <- paste(paste0('(', runmodel_plan_basic$target, '(?=([$].+)*))'),
                      collapse='|')

rename_arg <- function(arg, in_regex, in_suffix) {
  if (inherits(arg, 'name')) {
    argc <- as.character(arg)
    argc_match <- grep(in_regex, argc, value = T, perl = T)
    if (length(argc_match)==1) {
      argrep <- gsub(in_regex, paste0(argc_match, in_suffix), argc, perl=T)
      print(paste('Changing argument', argc, 'to', argrep))
      return(as.name(argrep))
    }
  } else if (inherits(arg, 'call')) {
    rename_call_args(call=arg,
                     in_regex = target_regex,
                     in_suffix = target_suffix)
  }
}

########## To continue troubleshot
rename_call_args <- function(call, in_regex, in_suffix) {
  args <- rlang::call_args(call)

  args_renamed <- lapply(args, rename_arg,
                         in_regex=in_regex, in_suffix=in_suffix) %>%
    plyr::compact()

  arg_names <- call_args_names(call)
  if (any(duplicated(arg_names)) & length(arg_names) == 2) { #If there are duplicate argument names (as in a primitive function like foo$bar)
    call_modified <- call2(call_fn(call),
                           args_renamed[[1]],
                           args_renamed[[2]])

  } else { #if all argument names are unique
    call_modified <- rlang::call_modify(call, !!!arg_renamed)
  }

  return(call_modified)
}

###TRY TO NEST THEM
bla <- runmodel_plan_u10$command[[3]]
rename_call_args(call=bla,
                 in_regex = target_regex,
                 in_suffix = target_suffix)


#Deal with primitives
hash <-
do.call(call_fn(call), list(quote(runmodel_plan_u10), quote(command)))


meh <- list(cha=1)
call2(call_fn(call)(as.name('meh'), as.name('cha')))

call_fn(call)(as.name('meh'), as.name('cha'))

gsub(paste0('seplearners', '(?=([)"$,]|\\[|\\s)))'),
     '1',
     runmodel_plan_u10$command[[1]],
     perl = T
)


#list

# Modify an existing argument
call_modify(call, na.rm = FALSE)
#> mean(x, na.rm = FALSE)
call_modify(call, x = quote(y))
#> mean(x, na.rm = TRUE, x = y)



runmodel_plan_o10 <- runmodel_plan_basic
runmodel_plan_o10$target <- paste0(runmodel_plan_u10$target, '_o10')

final_plan <- bind_plans(preprocess_plan, runmodel_plan_u10, runmodel_plan_o10)

drake::vis_drake_graph(final_plan, targets_only = T)

