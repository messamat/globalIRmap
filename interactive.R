library(drake)
#r_outdated()
r_make()

vis_drake_graph(plan, targets_only = T)
#drake::drake_cache("D:/Mathis/PhD/globalIRmap/src/globalIRmap/.drake")$unlock()
# drake::clean(list = cached_unplanned(plan), garbage_collection = TRUE)
