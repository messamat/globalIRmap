#Intermittent river analysis master script 
#Author: Mathis L. Messager
#Contact info: mathis.messager@mail.mcgill.ca
#Affiliation: 
#Global HydroLAB, Department of Geography, McGill University
#EcoFlows Lab, RiverLy Research Unit, INRAE Lyon

source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmap_plan.R')

future::plan(future::multiprocess) 
drake_config(plan, verbose = 2, parallelism = "future", jobs = )
make(plan)

vis_drake_graph(plan)

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