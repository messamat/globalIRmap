library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

memory.limit(size=50000)

drake_config(plan,
             verbose=1L,
             cache_log_file = TRUE,
             memory_strategy="preclean",
             garbage_collection = TRUE,
             prework = quote(future::plan(future.callr::callr)),
             log_make = "log/drake.log")
#,workers=availableCores()-1)))

#drake_history(plan)
#loadd()



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