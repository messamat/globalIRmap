library(drake)
source('R/IRmapping_packages.R')
source('R/IRmapping_functions.R')
source('R/IRmapping_plan.R')

loadd(rftuned)

outer_lrns <- rftuned$rf_outer$learners
#Get performance measure and train+test time for each resampling iteration
timeacc_summary <- lapply(seq_along(outer_lrns), function(outerf) {
  farchive <- outer_lrns[[outerf]]$tuning_instance$archive() 
  
  outdt <- lapply(seq_along(farchive$resample_result), function(innerf) {
    return(
      data.table(outf = outerf,
                 inf = innerf,
                 meantime=mean(farchive$resample_result[[innerf]]$score(msr('time_both'))$time_both),
                 meanclassifbac=mean(farchive$resample_result[[innerf]]$score(msr('classif.bacc'))$classif.bacc)) %>%
        cbind(as.data.table(farchive$params[[innerf]]), .)
    )
  }) %>%
    do.call(rbind, .)
  
  return(outdt)
}) %>%
  do.call(rbind, .) 

#Compute time and performance relative to mean within each inner fold
timeacc_summary[, `:=`(meantime_rel = meantime/mean(meantime),
                       meanbac_rel = meanclassifbac/mean(meanclassifbac)),
                by=.(outf)]

timeacc_summary_melt <- melt(timeacc_summary,
                             id.vars=colnames(timeacc_summary[, -c('classif.ranger.mtry',
                                                                   'classif.ranger.min.node.size',
                                                                   'classif.ranger.sample.fraction')]))

ggplot(timeacc_summary_melt, aes(x=value, y=meantime_rel, 
                                 group=as.factor(outf))) + 
  geom_point(aes(size = meanbac_rel, color=meanbac_rel)) + 
  geom_smooth(method='lm', se = FALSE) + 
  scale_color_distiller(palette='Spectral') +
  facet_wrap(~variable, scales='free') + 
  theme_bw()

ggplot(timeacc_summary_melt, aes(x=value, y=meanbac_rel, 
                                 group=as.factor(outf))) + 
  geom_point(aes(color=meantime_rel, size = meantime_rel)) + 
  geom_smooth(method='lm', se = FALSE) + 
  scale_color_distiller(palette='Spectral') +
  facet_wrap(~variable, scales='free')


ggplot(timeacc_summary_melt, aes(x=value, y=meantime, 
                                 group=as.factor(outf), color=as.factor(outf))) + 
  geom_point() + 
  facet_wrap(~variable, scales='free') +
  geom_smooth(method='lm')


ggplot(timeacc_summary, aes(x=classif.ranger.mtry, 
                            y=classif.ranger.sample.fraction)) + 
  geom_point(aes(color=meantime_rel)) +
  scale_color_distiller(palette='Spectral')

ggplot(timeacc_summary, aes(x=classif.ranger.mtry, 
                            y=classif.ranger.min.node.size)) + 
  geom_point(aes(color=meantime_rel)) +
  scale_color_distiller(palette='Spectral')

ggplot(timeacc_summary, aes(x=classif.ranger.sample.fraction, 
                            y=classif.ranger.min.node.size)) + 
  geom_point(aes(color=meantime_rel)) +
  scale_color_distiller(palette='Spectral')