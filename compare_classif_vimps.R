library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

loadd(rfbm)
loadd(predvars)

inbm_cforest <- rfbm$bm_classif$clone()$filter(learner_ids = "oversample.classif.cforest")
lrn_cforest <-  inbm_cforest$resample_result(uhash=unique(as.data.table(inbm_cforest)$uhash))

inbm_tunedranger <- rfbm$bm_classif$clone()$filter(learner_ids = "oversample.classif.ranger.tuned")
lrn_tunedranger <- inbm_tunedranger$resample_result(uhash=unique(as.data.table(inbm_tunedranger)$uhash))
  
vimp_cforest <- ggvimp(in_rftuned=lrn_cforest, in_predvars=predvars)
vimp_ranger <- ggvimp(in_rftuned=lrn_tunedranger, in_predvars=predvars)

grid.arrange(vimp_cforest, vimp_ranger)
