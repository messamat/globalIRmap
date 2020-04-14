library(drake)
#r_outdated()
r_make()

vis_drake_graph(plan, targets_only = T)
# drake::drake_cache("D:/Mathis/PhD/globalIRmap/globalIRmap/.drake")$unlock()
