library(drake)
#r_outdated()
r_make()

vis_drake_graph(plan, targets_only = T)

cached()
#drake::drake_cache("C:/globalIRmap/src/globalIRmap/.drake")$unlock()
# drake::clean(list = cached_unplanned(plan), garbage_collection = TRUE)

history_last <- drake_history(analyze = TRUE) %>%
  setDT %>%
  .[, built := as_datetime(built)] %>%
  setorder(built) %>%
  .[, .SD[.N,], by=target]

cache <- drake_cache()
rfbm_regr<- cache$get_value(history_last[target=='rfbm_regr',hash])
rfbm_classif<- cache$get_value(history_last[target=='rfbm_classif',hash])
bmcheck_classif <- cache$get_value(history_last[target=='bm_checked_rfbm_classif.bm_classif_rfbm_classif.measure_classif', hash])
bmcheck_regr<- cache$get_value(history_last[target=='bm_checked_rfbm_regr.bm_regr_rfbm_regr.meassure_regr', hash])
rfeval_featsel <- cache$get_value(history_last[target=='rfeval_featsel',hash])
