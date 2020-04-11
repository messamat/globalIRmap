write_preds <- function(filestructure, gaugep, gaugestats_format, rftuned, predvars) {
  
  
}


# ---- Output GRDC predictions as points ----
cols_toditch<- colnames(gaugep)[colnames(gaugep) != 'GRDC_NO']
st_write(obj=merge(gaugep, 
                   gaugestats_format[, !cols_toditch, with=F], by='GRDC_NO'),
         dsn=filestructure['out_gauge'], 
         driver = 'gpkg', 
         delete_dsn=T)

# ---- Generate predictions to map on river network ----
cols_tokeep <-  c("HYRIV_ID", predvars[!is.na(ID),varcode],
                'ele_mt_cav','ele_mt_uav',
                paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
                paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
                paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0)))

riveratlas <- fread_cols(filestructure['in_riveratlas'], 
                         cols_tokeep = cols_tokeep) %>%
  comp_derivedvar


predict(rftuned, type='response', riveratlas[1:10,])

#Predict model (should take ~10-15 minutes for 9M rows x 40 cols — chunk it up by climate zone to avoid memory errors)
for (clz in unique(riveratlas$clz_cl_cmj)) {
  print(clz)
  tic()
  riveratlas[!is.na(cly_pc_cav) & !is.na(cly_pc_uav) & clz_cl_cmj == clz, 
             predbasic800 := predict(ranger_basictrained, type='response', data=.SD)$predictions[,'1'], 
             .SDcols = c("HYRIV_ID", rfpredcols)]
  toc()
}

riveratlas[, predbasic800cat := ifelse(predbasic800>=0.3, 1, 0)]
fwrite(riveratlas[, c('HYRIV_ID', 'predbasic800cat'), with=F], file.path(resdir, 'RiverATLAS_predbasic800.csv'))









##################################################################################################################################################
############################ MORE ALREADY SPED UP/TRANSCRIBED ####################################################################################
#Set default ranger learner with explicit parameter set
# rangerlrn <- mlr3::lrn('classif.ranger', id = 'ranger', 
#                        num.trees = 500,
#                        #mtry = , #Default is square root of number of variables
#                        min.node.size = 1, #Default is 1 for classification, 5 for regression and 10 for probability
#                        replace = FALSE, #in ranger package, default replace is True but Boulesteix et al. 2012 show that less biased
#                        sample.fraction = 0.632, #Default for sampling without replacement
#                        split.select.weights = ,
#                        always.split.variables= ,
#                        respect.unordered.factors = 'ignore',
#                        importance='permutation',
#                        write.forest=TRUE, 
#                        scale.permutation.importance = FALSE,
#                        num.threads = bigstatsr::nb_cores(),
#                        save.memory = FALSE,
#                        verbose = TRUE,
#                        splitrule = "gini",
#                        num.random.splits = 1L,
#                        keep.inbag = FALSE,
#                        predict_type = 'prob',
#                        max.depth = NULL #Unlimited depth default
# )



#---- Get mDur and mFreq in winter and non-winter periods to assess whether intermittency is due to freezing ----

#Compute 3-month period of minimum temperature
window <- 3
mintemplist<- dcast(
  melt(
    reaches_with_points[, c('GRDC_NO', grep('^tmp.*_c[0-9]{2}', colnames(reaches_with_points), value = T)), with=F], 
    id.var = 'GRDC_NO'),
  variable~GRDC_NO, value.var='value')[,-1] %>% #Pre-format temperature
  frollapply(n=window, FUN=mean, align="center") %>% #Compute 3-month rolling mean
  lapply(function(x) {
    mincenter <- which(x==min(x,na.rm=T)) #Get center month index of 3-month period with minimum temp
    return(seq(mincenter-floor(window/3), mincenter+ceiling(window/3)))}) #Recreate window of months
names(mintemplist) <- reaches_with_points$GRDC_NO

#Compute mdur and mfreq for winter months
gaugestats_winter <- durfreq_parallel(pathlist=fileNames, maxgap=20, monthsel_list=mintemplist, reverse=FALSE)

#Compute mdur and mfreq for non-winter months
gaugestats_nowinter <- durfreq_parallel(pathlist=fileNames, maxgap=20, monthsel_list=mintemplist, reverse=TRUE)

#---- Check month of intermittency ----
monthinter_melt <- melt(gaugestats[, c('GRDC_NO', "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), with=F],
                        id.var='GRDC_NO')
ggplot(monthinter_melt, aes(x=variable, y=value, fill=GRDC_NO)) +
  geom_area(stat='identity')

#Base ranger training
#To adapt
OOBfulldat <- OOBCurve(mod=mod_basic, measures = list(mmce, auc, brier), task = intermittent.task, data = gaugestats_format[!is.na(cly_pc_cav),]) %>%
  setDT %>%
  .[, numtrees := .I]

#Plot it
ggplot(melt(OOBfulldat, id.vars = 'numtrees'), aes(x=numtrees, y=value, color=variable)) + 
  geom_point(size=1, alpha=0.5) + 
  scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.1)) + 
  scale_x_continuous(expand=c(0,0), breaks = seq(0,1000,100)) + 
  theme_bw()
