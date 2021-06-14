library(drake)
#r_outdated()
r_make()


drake::drake_cache("C:\\globalIRmap\\src/globalIRmap/.drake")$unlock()
drake::drake_cache("E:\\Mathis/src/globalIRmap/.drake")$unlock()

#################### Interactive manipulations of plan #########################
source("_drake.R")
drake::vis_drake_graph(plan, targets_only = T)
drake::sankey_drake_graph(plan, targets_only = T)

cached()

history_last <- drake_history(analyze = FALSE) %>%
  setDT %>%
  .[, built := as_datetime(built)] %>%
  setorder(built) %>%
  .[, .SD[.N,], by=target]

cache <- drake_cache()

rfbm_regr<- cache$get_value(history_last[target=='rfbm_regr',hash])
rfbm_classif<- cache$get_value(history_last[target=='rfbm_classif',hash])
bmcheck_classif <- cache$get_value(
  history_last[target=='bm_checked_rfbm_classif.bm_classif_rfbm_classif.measure_classif', hash])
bmcheck_regr<- cache$get_value(
  history_last[target=='bm_checked_rfbm_regr.bm_regr_rfbm_regr.meassure_regr', hash])
rfeval_featsel <- cache$get_value(history_last[target=='rfeval_featsel',hash])
rftuned <- cache$get_value(history_last[target=='rftuned',hash])

library(vctrs)

loadd(baselearners)
# same as readd(model)
s <- subtargets(baselearners)[1]
check <- vec_c(
  readd(s[1], character_only = TRUE),
  readd(s[2], character_only = TRUE),
  readd(s[3], character_only = TRUE)
)

#drake::clean(destroy = F, garbage_collection = T)

# "recover=T
# drake’s data recovery feature is another way to avoid rerunning commands. It is useful if:
#
#   You want to revert to your old code, maybe with git reset.
# You accidentally clean()ed a target and you want to get it back.
# You want to rename an expensive target."

# drake_config(plan,
#   verbose = 2,
#   targets = c("pathogen_maps_debugging", "prediction_pathogens"),
#   lazy_load = "promise",
#   console_log_file = "log/drake.log",
#   caching = "worker",
#   template = list(log_file = "log/worker%a.log", n_cpus = 16, memory = 60000,
#     job_name = "paper1"),
#   prework = list(quote(set.seed(1, "L'Ecuyer-CMRG")),
#     quote(future::plan(future.callr::callr, workers = 10)),
#     quote(parallelStart(
#       mode = "multicore", cpus = ignore(16), level = "mlr.resample"
#     ))
#   ),
#   garbage_collection = TRUE, jobs = 3, parallelism = "clustermq",
#   lock_envir = FALSE, keep_going = TRUE
# )


