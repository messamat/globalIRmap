library(drake)
source('R/IRmapping_packages.R')

#### Choose whether to run full plan or simplified plan to get main results ######
plan_choice <-'full' #'trimmed'


if (plan_choice == 'trimmed') {
  source('R/IRmapping_plan_trimmed.R')
} else {
  source('R/IRmapping_plan.R')
}

drake_config(plan,
             verbose=1L,
             cache_log_file = TRUE,
             memory_strategy="preclean",
             garbage_collection = TRUE,
             seed = 0,
             prework = quote(future::plan(future.callr::callr)),
             log_make = "log/drake.log"
)

