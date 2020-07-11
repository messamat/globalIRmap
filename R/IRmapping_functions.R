#Functions for intermittent river analysis

##### -------------------- Utility functions ---------------------------------

#------ diny -----------------
#' Number of days in the year
#'
#' Computes the number of days in a given year, taking in account leap years
#'
#' @param year Numeric or integer vector of the year (format: 4-digit year, %Y).
#'
#' @return Number of days in that year.
#'
#' @examples
#' diny(1999)
#' diny(2000)
#' diny(2004)
#' diny(2100)
#' diny(1600)
#'
#' @export
diny <- function(year) {
  365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
}

#------ zero_lomf -----------------
#' Last \[non-zero\] Observation Moved Forward (lomf)
#'
#' Finds the index, for each row, of the previous row with a non-zero value
#'
#' @param x Numeric vector.
#' @param first (logical) Whether to consider first value as a non-zero value
#'   whose index is moved forward even if it is zero. This prevents having NAs
#'   in the results and somewhat assumes that, for a time series, the day prior
#'   to the first value is non-zero.
#'
#' @return Numeric vector of the indices of the previous non-zero for each
#'   element of the input vector.
#'
#' @examples
#' test1 <- c(1,1,1,0,0,0,0,1,1)
#' zero_lomf(test)
#' test2 <- c(0,0,0,0,0,1,1,0,1)
#' zero_lomf(test2, first=FALSE)
#' zero_lomf(test2, first=TRUE)
#'
#' @export
zero_lomf <- function(x, first=TRUE) {
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if(first==T & x[1]==0)
      non.zero.idx=c(1,non.zero.idx)
    #Repeat index of previous row with non-zero as many times gap until next non-zero values
    rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
  }
}

#------ comp_ymean -------------
#Compute weighted mean across 12 columns corresponding to 12 months
#mstart is the start of the two-number month
comp_ymean<- function(in_dt, fieldex, mstart) {
  refd <- data.table(month=c('01', '02', '03', '04', '05', '06',
                             '07', '08', '09', '10', '11', '12'),
                     dimonth = c(31, 28.25, 31, 30, 31, 30,
                                 31, 31, 30, 31, 30, 31))
  refd[, mcols := paste0(substr(fieldex, 1, mstart-1),
                         month,
                         substr(fieldex, mstart+2, nchar(fieldex)))]
  in_dt2 <- copy(in_dt[, refd$mcols, with=F])
  in_dt2[, comp_ymeanID := .I]
  in_dt[, outcol] <- melt(in_dt2,
                          id.vars='comp_ymeanID',
                          variable.name = 'mcols') %>%
    .[refd, on='mcols'] %>%
    .[, weighted.mean(value, dimonth), by=comp_ymeanID] %>%
    .[,'V1',with=F]
}

#------ comp_derivedvar -----------------
#' Compute derived variables
#'
#' Format and compute derived a set of environmental variables from input
#' dataset that contains
#' \href{https://www.hydrosheds.org/page/hydroatlas}{HydroATLAS} variables for
#' the purpose of predicting intermittency.
#'
#' @param dt data.table which contains rows as records and HydroATLAS
#'   environmental variables as columns.
#' @param copy (logical) whether to make a copy of \code{dt} (if TRUE) or to
#'   modify it in place (if FALSE)
#'
#' @details This function only formats and compute derived variables deemed
#'   relevant for the intermittency analysis. \cr
#'   \cr
#'   Steps include: \cr
#'   1. Convert some -9999 that should be 0 after checking them
#'   2. Count the number of -9999 values per column
#'   3. Convert the rest of the -9999 values to NA
#'   4. Compute a set of derived variables based on existing variables in
#'   HydroATLAS
#'   (e.g. \code{pre_mm_cvar= fifelse(pre_mm_cmx==0, 0, pre_mm_cmn/pre_mm_cmx)})
#'   5. Correct scaling from RiverATLAS
#'   (e.g. Degree of regulation is out of 10,000 in HydroATLAS)
#'
#' ```
#' @source Linke, S., Lehner, B., Dallaire, C. O., Ariwi, J., Grill, G., Anand,
#'   M., ... & Tan, F. (2019). Global hydro-environmental sub-basin and river
#'   reach characteristics at high spatial resolution. Scientific Data, 6(1),
#'   1-15.
#'
#' @export

comp_derivedvar <- function(in_dt, copy=FALSE) {
  if (copy) {
    in_dt2 <- copy(in_dt)
  } else {
    in_dt2 <- in_dt
  }


  #---- Inspect and correct -9999 values ----
  print('Inspect and correct -9999 values')
  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  in_dt2[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  in_dt2[snw_pc_cmx == -9999, snw_pc_cmx:=0]

  print('Number of NA values per column')
  colNAs<- in_dt2[, lapply(.SD, function(x) sum(is.na(x) | x==-9999))]
  print(colNAs)

  #Convert -9999 values to NA
  gladcols <- c('wloss_pc_cav', 'wdryp_pc_cav', 'wwetp_pc_cav',
                   'whfrq_pc_cav', 'wseas_pc_cav', 'wperm_pc_cav',
                   'wfresh_pc_cav',
                   'wloss_pc_uav', 'wdryp_pc_uav', 'wwetp_pc_uav',
                   'whfrq_pc_uav', 'wseas_pc_uav', 'wperm_pc_uav',
                   'wfresh_pc_uav')
  sgcols <- c('cly_pc_cav','slt_pc_cav', 'snd_pc_cav',
              'cly_pc_uav','slt_pc_uav', 'snd_pc_uav')
  biocols <- unlist(lapply(c('cav', 'uav'),
                           function(s) paste0('bio', 1:11, '_dc_', s)
                           ))
  biocolsnegative <- grep('bio[15689][01]*_.*', biocols, value=T)

  for (j in which(sapply(in_dt2,is.numeric))) { #Iterate through numeric column indices
    if (!(j %in% c(gladcols, sgcols))) {
      set(in_dt2,which(in_dt2[[j]]==-9999),j, NA) #Set those to 0 if -9999
      }
    }

  cmi_cmcols <- paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))
  cmi_umcols <- paste0('cmi_ix_u', str_pad(1:12, width=2, side='left', pad=0))

  #Scale variables based on HydroATLAS v1.0 documentation and v1.1 processing
  in_dt2[, `:=`(
    ari_ix_cav = ari_ix_cav/100,
    ari_ix_uav = ari_ix_uav/100,
    dor_pc_pva = dor_pc_pva/100,
    lka_pc_cse = lka_pc_cse/10,
    lka_pc_use = lka_pc_use/10,
    gwt_m_cav = gwt_cm_cav/100
  )]

  in_dt2[, (gladcols) := lapply(.SD, function(x) x/100), .SDcols = gladcols]
  in_dt2[, (biocols) := lapply(.SD, function(x) (x/100)), .SDcols = biocols]
  in_dt2[, (biocolsnegative) := lapply(.SD, function(x) x-100), .SDcols=biocolsnegative]
  in_dt2[, (cmi_cmcols) := lapply(.SD, function(x) x-100), .SDcols=cmi_cmcols]
  in_dt2[, (cmi_umcols) := lapply(.SD, function(x) x-100), .SDcols=cmi_umcols]

  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  comp_ymean(in_dt=in_dt2, fieldex = 'cmi_ix_c01', mstart=8, outcol='cmi_ix_cyr')
  comp_ymean(in_dt=in_dt2, fieldex = 'cmi_ix_u01', mstart=8, outcol='cmi_ix_uyr')
  comp_ymean(in_dt=in_dt2, fieldex = 'pet_mm_c01', mstart=8, outcol='pet_mm_cyr')
  comp_ymean(in_dt=in_dt2, fieldex = 'pet_mm_u01', mstart=8, outcol='pet_mm_uyr')

  in_dt2[, `:=`(cmi_ix_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
                cmi_ix_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
         .SDcols= cmi_cmcols] %>%
    .[, `:=`(cmi_ix_umn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
             cmi_ix_umx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
      .SDcols= cmi_umcols] %>%
    .[, swc_pc_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
      .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly swc
    .[, `:=`(#min/max monthly CMI
             cmi_ix_cvar= fifelse(cmi_ix_cmx==0, 0, cmi_ix_cmn/cmi_ix_cmx),
             cmi_ix_uvar= fifelse(cmi_ix_umx==0, 0, cmi_ix_umn/cmi_ix_umx),
             #min/max monthly watershed discharge
             dis_m3_pvar=fifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx),
             #min monthly/average yearly watershed discharge
             dis_m3_pvaryr=fifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),
             #catchment average elv - watershec average elev
             ele_pc_rel = fifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav),
             #runoff coefficient (runoff/precipitation)
             runc_ix_cyr = run_mm_cyr/pre_mm_cyr,
             #Specific discharge
             sdis_ms_uyr = dis_m3_pyr/UPLAND_SKM,
             sdis_ms_umn = dis_m3_pmn/UPLAND_SKM
             )]
  return(in_dt2)
}

#------ threshold_misclass -----------------
#' Threshold misclassification
#'
#' Compute the misclassification rate, sensitivity, and specificity of probability
#' predictions to a binary classification problem from a random forest learner
#' based on a given threshold.
#'
#' @param i (numeric) Threshold to be used for converting probabilities
#'   predictions to binary predictions, limits=\[0,1\].
#' @param in_preds Either:
#' 1. a \link[mlr3]{PredictionClassif} or
#' 2. a data.table of predictions for a set of CV repetitions as formatted by
#' \code{\link{analyze_benchmark}}. The column of the predicted probability of
#' getting a 'positive' must be called \code{prob.1}.
#'
#' @return A single-row data.table with four columns:
#' \itemize{
#'   \item i - the threshold used to compute the classification statistics
#'   \item misclass - misclassification rate \[0-1\]
#'   \item sens - sensitivity \[0-1\]
#'   \item spec - specificity \[0-1\]
#' }
#'
#' @details
#' \itemize{
#' \item Misclassification rate is the proportion of misclassified
#' records
#' \item Sensitivity is the proportion of actual positives (1) that are
#' correctly identified as such (the proportion of intermittent rivers that
#' are identified as intermittent)
#' \item Specificity is the proportion of actual negatives that are correctly
#' identified as such (the proportion of perennial rivers that are identified as
#' perennial)
#' } \cr
#' See \url{https://en.wikipedia.org/wiki/Sensitivity_and_specificity}
#'
#' @section Warning:
#' This function was only tested for the outputs of a probability random forest
#' (using a classification framework) and a regression forest (using 0 or 1 as
#' dummy continuous variables)
#'
#' @export

threshold_misclass <- function(i=0.5, in_preds) {
  #---- Get confusion matrix ----
  if (inherits(in_preds, 'PredictionClassif')) {
    confu <- as.data.table(in_preds$set_threshold(1-i)$confusion) #Get confusion matrix directly
  }

  if (is.data.table(in_preds)) {
    #If task associated with predictions is a classification and has records
    if (in_preds[task_type == 'classif',.N] > 0) {
      #For each CV repetition:
      #   1. set the probability threshold to compute a confusion matrix to 1-i
      #     (i being the threshold to classify something as 1, set_threshold
      #     being based on prob.0, not prob.1)
      #   2. Compute confusion matrix
      confu <- in_preds[, as.data.table(pred[[1]]$set_threshold(1-i)$confusion),
                        by=outf] %>%
        #Aggregate confusion matrices across repetitions
        .[, .(N=sum(N)), by=.(response, truth)]
    }

    #If task associated with predictions is a regression and has records
    if (in_preds[task_type == 'regr', .N] > 0) {
      #Reclassify continuous predictions into binary response across all records
      confu <- in_preds[, response := fifelse(prob.1>=i, '1', '0')] %>%
        .[, truth := as.character(truth)] %>%
        #Create aggregate confusion matrix
        .[, .N, by=.(response, truth)]
    }
  }

  #---- Compute statistics based on confusion matrix and format into data.table----
  outvec <- data.table(
    i,
    misclas = confu[truth != response, sum(N)] / confu[, sum(N)],
    sens = confu[truth == '1' & response == '1', N] / confu[truth=='1', sum(N)],
    spec  = confu[truth=='0' & response==0, N]/confu[truth=='0', sum(N)])
  return(outvec)
}

#------ threshold_dat ----------------
threshold_dat <- function(bmres) {
  bmres_dt <- as.data.table(bmres) %>%
    .[, learner_id := unlist(lapply(learner, function(x) x$id))] %>%
    .[, task_id := unlist(lapply(task, function(x) x$id))] %>%
    .[, resampling_id := unlist(lapply(resampling, function(x) x$id))] %>%
    unique(by='uhash')

  if (bmres$task_type == 'regr') {
    preds <- lapply(seq_len(bmres_dt[,.N]), function(rsmp_i) {
      preds <- bmres_dt$prediction[[rsmp_i]]$test %>%
        as.data.table %>%
        .[, `:=`(outf = bmres_dt$iteration[[rsmp_i]],
                 task = bmres_dt$task[[rsmp_i]]$id,
                 task_type = bmres$task_type,
                 learner = bmres_dt$learner[[rsmp_i]]$id)]
      return(preds)
    }) %>%
      do.call(rbind, .)

    if (!('prob.1' %in% names(preds)) & 'response' %in% names(preds)) {
      preds[, prob.1 := response]
    }
  }

  outtab <- lapply(1:bmres$n_resample_results, function(rsmp_i) {
    print(rsmp_i)
    rsmp_preds <- bmres$resample_result(rsmp_i)$prediction()

    if (bmres$task_type == 'regr') {
      rsmp_preds <- preds
    }

    baccthresh <- ldply(seq(0,1,0.01), function(threshold_class) {
      print(threshold_class)
      cbind(
        rsmp_bacc(bmres, rsmp_i=rsmp_i, threshold_class=threshold_class),
        threshold_misclass(i=threshold_class, in_preds=rsmp_preds)
      ) %>%
        .[, i:=NULL]
    }) %>%
      setDT %>%
      cbind(bmres_dt[rsmp_i, .(learner_id, task_id, resampling_id)])
    return(baccthresh)
  }) %>%
    rbindlist


  outtab_melt <- melt(outtab, id.vars=c('threshold_class',
                                        'learner_id', 'task_id', 'resampling_id'))
  return(outtab_melt)
}
#------ get_outerrsmp -----------------
#' Get outer ResampleResult
#'
#' Extract the outer resampling \link[mlr3]{ResampleResult} from the output
#' from \code{\link{selecttrain_rf}} with the option to grab the output from
#' the spatial resampling strategy or not (when there are multiple resampling
#' strategies).
#'
#' @param in_rftuned (list or \link[mlr3]{ResampleResult}) Output from
#' \code{\link{selecttrain_rf}}.
#' @param spatial_rsp (logical) Whether to return the results from the
#' spatial resampling strategy if there is one.
#'
#' @details
#' If \code{in_rftuned} is already a \link[mlr3]{ResampleResult}, then simply
#' return it.
#'
#' @return a \link[mlr3]{ResampleResult}.
#'
#' @section Warning:
#' This will only return the \link[mlr3]{ResampleResult} from a single resampling
#' strategy.
#'
#' @export

get_outerrsmp <- function(in_rftuned, spatial_rsp=FALSE) {
  #Adapt whether return resample result or output from selecttrain_rf
  if (inherits(in_rftuned, 'list')) {
    #If there is more than one type of outer resampling
    if (length(in_rftuned$rf_outer) > 0) {
      #Check which resampling is spatial — only works if one is spatial
      sp_i <- which(unlist(lapply(in_rftuned$rf_outer, function(x) {
        grepl('.*Resampling.*Sp.*', x$resampling$format())
      })))

      #If user request that spatial resampling be used
      if (spatial_rsp==TRUE) {

        if (length(sp_i)>0) {
          rsmp_res <- in_rftuned$rf_outer[[min(sp_i)]]
        } else { #But if there is no spatial resampling provided
          stop("spatial_rsp==TRUE but the in_rftuned does not include
               any Spatial Resampling")
        }
        #If user didn't request spatial resampling to be used, grab the first
        #resampling that is not spatial
      } else {
        rsmp_res <- in_rftuned$rf_outer[[
          min((1:length(in_rftuned$rf_outer))[-sp_i])]]
      }

      #If there is only one type of outer resampling
    } else {
      rsmp_res <- in_rftuned$rf_outer
    }
    #If in_rftuned is alreadyh a ResampleResult, simply return it
  } else if (inherits(in_rftuned, "ResampleResult")) {
    rsmp_res <- in_rftuned
  }
  return(rsmp_res)
}

#------ weighted_sd -----------------
#' Weighted standard deviation
#'
#' Computes a weighted standard deviation
#'
#' @param x An object containing the values whose weighted mean is to be
#' computed.
#' @param w An integer vector of weights the same length as x giving the
#' weights to use for elements of x (weights are > 1, not decimals).
#' @param na.rm	(logical) Whether NA values in x should be stripped before the
#' computation proceeds.
#'
#' @return A length-one numeric vector.
#'
#' @examples
#' ```{r}
#' wt <- c(5, 4, 3 ,2 , 1)
#' x <- c(1, 2, 3, 4, 5)
#' xm <- weighted.mean(x, wt)
#' xsd <- weighted_sd(x, wt)
#' ```
#'
#' @export

weighted_sd <- function(x, w=NULL, na.rm=FALSE) {
  if (na.rm) {
    x <-  na.omit(x)
    if (length(w) > length(x)) {
      w <- w[-which(is.na(x))]
    }
  }

  if (length(w)==0) {
    w <- rep(1, length(x))
  }

  #Compute weighted standard deviation
  return(sqrt(sum((w) * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)))
}

#------ extract_impperf_nestedrf -----------------
#' Extract variable importance and performance from a trained RF learner
#'
#' Computes the performance, variable importance and associated p-value from
#' either a trained \link[mlr3]{AutoTuner} or a trained \link[mlr3]{Learner}.
#'
#' @param in_rflearner An \link[mlr3]{AutoTuner} or a trained
#' \link[mlr3]{Learner}.
#' @param imp (logical) Whether to compute variable importance.
#' @param perf (logical) Whether to compute the performance measure.
#' @param pvalue (logical) Whether to compute p-values.
#' @param pvalue_permutn (integer) number of permutations to use in p-value
#' calculations
#'
#' @return A data.table with, when all logicals are TRUE, the following columns.
#' \describe{
#'   \item{varnames} - predictor variable name
#'   \item{imp} - variable importance value
#'   \item{pvalue} - variable importance p-value
#'   \item{\[perf\]} - performance value (same for all rows), name changes
#' }
#'
#' @details
#' Accepts both a trained \link[mlr3]{AutoTuner} or a trained \
#' link[mlr3]{Learner} of class \link[mlr3learners]{mlr_learners_classif.ranger}
#' or \link[mlr3learners.partykit]{mlr_learners_classif.cforest}. \cr
#'
#' Also accept learners of class \link[mlr3pipelines]{GraphLearner} \cr
#'
#' If \code{p-value==TRUE} and \link[mlr3learners]{mlr_learners_classif.ranger},
#' then compute p-values with \link[ranger]{importance_pvalues} with Altmann
#' permutation method using \code{pvalue_permutn} permutations. \cr
#'
#' If \link[mlr3learners.partykit]{mlr_learners_classif.cforest} is provided,
#' \code{p-value} is ignored for now.
#'
#' @seealso \link{weighted_vimportance_nestedrf}
#' @section documentation to-do:
#' Can add an example down the line, add source.
#' @export

extract_impperf_nestedrf <- function(in_rflearner,
                                     imp = TRUE, perf = TRUE,
                                     pvalue = TRUE, pvalue_permutn = 100) {

  if (inherits(in_rflearner, "AutoTuner")) {
    sublrn <- in_rflearner$model$learner
  } else {
    sublrn <- in_rflearner
  }

  print(paste0("Computing variable importance for resampling instance hash #",
               sublrn$hash))

  return(cbind(if (imp) {
    ####################### IF GraphLearner ####################################
    if (inherits(sublrn, "GraphLearner")) {

      if ('classif.ranger' %in% names(sublrn$model)) {

        if (pvalue == TRUE) {
          in_task <- in_rflearner$model$tuning_instance$task
          in_formula <- as.formula(paste0(in_task$target_names, '~.'))

          importance_pvalues(
            sublrn$model$classif.ranger$model,
            method = "altmann",
            num.permutations = pvalue_permutn,
            data = in_task$data(),
            formula= in_formula
          )

        } else {
          data.table(importance=sublrn$model$classif.ranger$model$variable.importance)
        }
      }

      else if ('classif.cforest' %in% names(sublrn$model)) {
        if (pvalue == TRUE) {
          warning("p_value calculation is only available for ranger classification rf, ignoring p_value.
                  In addition, default parameters were used in partykit::varimp, adjust as needed.")
        }
        data.table(importance=
                     partykit::varimp(sublrn$model$classif.cforest$model,
                                      nperm = 1,
                                      OOB = TRUE,
                                      risk = "misclassification",
                                      conditional = FALSE,
                                      threshold = .2))
      }
    } else { ####################### IF direct model ####################################
      if (pvalue == TRUE) { #If want pvalue associated with predictor variables
        in_task <- in_rflearner$model$task
        in_formula <- as.formula(paste0(in_task$target_names, '~.'))

        importance_pvalues(
          in_rflearner$model,
          method = "altmann",
          num.permutations = pvalue_permutn,
          data = in_task$data(),
          formula= in_formula
        )
      } else { #If pvalue == FALSE
        data.table(importance= in_rflearner$model$learner$importance())
      }

    }
  },
  if (perf) {
    outperf <- in_rflearner$tuning_result$perf
    data.table(outperf) %>% setnames(names(outperf))
  }
  ))
}
#------ weighted_vimportance_nestedrf -----------------
#' Weighted mean of variable importance for resampled RF learner
#'
#' Compute mean and standard deviation of variable importance and mean of
#' p-value across resampling instances (e.g. folds and repetitions) weighted by
#' resampling prediction accuracy.
#'
#' @param rfresamp the \link[mlr3]{ResampleResult} from a classification RF
#' of type \link[mlr3learners]{mlr_learners_classif.ranger}
#' or \link[mlr3learners.partykit]{mlr_learners_classif.cforest}
#' @inheritParams extract_impperf_nestedrf
#'
#' @return A data.table with, when all logicals are TRUE, the following columns.
#' \describe{
#'   \item{varnames} - predictor variable name
#'   \item{imp_wmean} - weighted mean of variable importance across resampling
#'   instances
#'   \item{imp_wsd} - weighted standard deviation of variable importance across
#'   resampling instances
#'   \item{pvalue_wmean} -weighted mean of variable p-value across resampling
#'   instances
#' }
#'
#' @details
#' See \link{extract_impperf_nestedrf} for more details on computations.
#'
#' @seealso \link{extract_impperf_nestedrf} and \link{ggvimp} and
#' \link{benchmark_featsel}
#'
#' @section Warning:
#' Does not accept error measure for weighting
#' - could be amended to be based on 1-error
#'
#' @section documentation to-do:
#' Can add an example down the line, add source.
#' @export

weighted_vimportance_nestedrf <- function(rfresamp,
                                          pvalue = TRUE, pvalue_permutn = 100) {
  varnames <- rfresamp$task$feature_names

  vimportance_all <- lapply(rfresamp$learners, #Extract vimp and perf for each resampling instance
                            extract_impperf_nestedrf,
                            imp=T, perf=T, pvalue=pvalue, pvalue_permutn) %>%
    do.call(rbind, .) %>%
    cbind(., varnames)

  ####!!!!!!!!!!Adapt to allow for other measure than classif.bacc!!!!!!!!######
  out_vimportance <- vimportance_all[
    , list(imp_wmean = weighted.mean(importance, classif.bacc), #Compute weighted mean
           imp_wsd =  weighted_sd(importance, classif.bacc)), #Compute weighted sd
    by=varnames]

  if (pvalue) {
    out_vimportance <- cbind(
      out_vimportance,
      vimportance_all[,
                      list(imp_pvalue = weighted.mean(pvalue, classif.bacc)), #Compute weighted mean of pvalue
                      by=varnames][, !'varnames']
    )
  }

  return(out_vimportance)
}

#------ extract_pd_nestedrf -----------------
#' Extract partial dependence (and performance) from a trained RF learner.
#'
#' Computes the marginal relationship between a subset of the predictors
#' (here, two variables at a time) and the model’s predictions by averaging
#' over the marginal distribution of the compliment of this subset of
#' the predictors, taking in account the interaction between the chosen
#' predictors.
#'
#' @param learner_id (integer) Index of the outer resampling instance to be
#' analyzed.
#' @param in_rftuned \link[mlr3]{ResampleResult} from a classification RF.
#' @param datdf Data from the task that was used to train RF.
#' @param selcols Character vector of the predictor variables to analyze.
#' @param ngrid (integer) Number of values of the  predictor variables over
#' which to compute the marginal relationship.
#'
#' @return A data.table with the following columns.
#' \describe{
#'   \item{value1} - value of the first predictor variable in the pair
#'   \item{value2} - value of the second predictor variable in the pair
#'   \item{0} - predicted probability of 0 (e.g. probability that the river
#'   is perennial) at value1 and value2
#'   \item{1} - predicted probability of 1 (e.g. probability that the river
#'   is intermittent) at value1 and value2
#'   \item{\[perf\]} - performance value (same for all rows), name changes
#'   \item{var1} - name of the first predictor variable
#'   \item{var2} - name of the second predictor variable
#' }
#'
#'
#' @details
#' Also accept learners of class \link[mlr3pipelines]{GraphLearner}. \cr
#' Uses \link[edarf]{partial_dependence} for computing.
#'
#' @seealso \link{weighted_vimportance_nestedrf},
#'  \link{ggpd_bivariate}
#'
#' @section Warning:
#' Has only been tested on \link[mlr3learners]{mlr_learners_classif.ranger}
#'
#' @section documentation to-do:
#' Can add an example down the line, add source.
#'
#' @export

extract_pd_nestedrf <- function(learner_id=1, in_rftuned, datdf,
                                selcols, nvariate, ngrid) {
  in_mod <- in_rftuned$learners[[learner_id]]
  #in_mod <- nestedresamp_ranger$learners[[1]]

  if (inherits(in_mod$learner, "GraphLearner")) {
    in_fit <- in_mod$learner$model$classif.ranger$model
  } else {
    in_fit <- in_mod$learner$model
  }

  #Get fold-specific performance measure
  foldperf <- extract_impperf_nestedrf(in_mod, imp=F, perf=T, pvalue=F)

  # selcols <- in_vimp_plot$data %>% #Can use that if extracting from tunredrf is expensive
  #   setorder(-imp_wmean) %>%
  #   .[colnums, variable]

  #Make dataset of all combinations of selected column names, two at a time
  if (nvariate == 1) {





  } else if (nvariate == 2) {
    vargrid <- combn(selcols, 2, simplify=F) %>%
      do.call(rbind, .)

    ngridvec <- c(ngrid, ngrid)

    #Get marginal distribution of the effect of two columns at a time
    pdcomb <- mapply(function(i, j) {
      pdout <- edarf::partial_dependence(in_fit, vars = c(i, j), n = ngridvec,
                                         interaction = TRUE, data = datdf) %>% #Warning: does not work with data_table
        setDT %>%
        .[,(names(foldperf)) := foldperf] %>%
        .[, `:=`(var1=i, var2=j)] %>%
        setnames(c(i,j), c('value1', 'value2'))


      return(pdout)
    }, vargrid[,1], vargrid[,2], SIMPLIFY = FALSE) %>%
      do.call(rbind, .)
  } else {
    print('Warning: function cannot yet work with more than two variables at a time')
  }

  return(pdcomb)
}

#------ fread_cols -----------------
#' fread columns
#'
#' Fast data.table-based reading of a subset of columns from a table
#'
#' @param file_name path of the table to be read
#' @param cols_tokeep character vector, names of columns to read
#'
#' @return a data.table with the \code{cols_tokeep} that are found in the table
#'
#' @examples
#' fread_cols(iris, c('Sepal.Length', 'Sepal.Width')
#'
#' @export

fread_cols <- function(file_name, cols_tokeep) {
  #Only read the first row from the file
  header <- fread(file_name, nrows = 1, header = FALSE)
  #Check which columns are in the table
  keptcols <- cols_tokeep[cols_tokeep %chin% unlist(header)]
  missingcols <- cols_tokeep[!(cols_tokeep %chin% unlist(header))]
  paste('Importing', file_name, 'with ', length(keptcols),
        'columns out of ', length(cols_tokeep), 'supplied column names')
  #Reading in table
  paste(missingcols, 'columns are not in the file...')
  return(fread(input=file_name, header=TRUE, select=keptcols, verbose=TRUE))
}

#------ get_oversamp_ratio -----------------
#' Get oversample ratio
#'
#' Identify minority class and compute ratio between the number of observations
#' in the majority and the minority classes of a binary
#' \link[mlr3]{TaskClassif }.
#'
#' @param in_task binary \link[mlr3]{TaskClassif }
#'
#' @return named list with following items:
#' \describe{
#' \item{minoclass} - Value of minority class (e.g. '1' for intermittent rivers)
#' \item{ratio} - numeric ratio between number of items in majority and minority
#' class
#' }
#'
#' @examples
#' ```{r}
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "in_dt", backend = in_dt, target ='intermittent')
#' get_oversamp_ratio(task)
#' ```
#'
#' @export

#When given an mlr3 task with binary classification target, gets which class is minority and the ratio


get_oversamp_ratio <- function(in_task) {
  return(
    in_task$data()[, .N, by=get(in_task$target_names)] %>%
      setorder(N) %>%
      .[, list(minoclass=get[1], ratio=N[2]/N[1])]
  )
}

#------ convert_clastoregrtask -----------------
#' Convert classification task to regression task
#'
#' Convert a binary \link[mlr3]{TaskClassif} or
#' \link[mlr3spatiotempcv]{TaskClassifST}
#' to a \link[mlr3]{TaskRegr}, optionally oversampling the minority class first.
#'
#' @param in_task \link[mlr3]{TaskClassif} or \link[mlr3spatiotempcv]{TaskClassifST}
#' @param in_id the id of the output task
#' @param oversample (logical) whether to oversample minority class beforehand
#' (this is not exactly equivalent to a PipeOp) but oversample PipeOp is not yet
#' available for regression tasks. If TRUE, will oversample so that classes are
#' equally represented in task.
#'
#' @return a \link[mlr3]{TaskRegr}
#'
#' @examples
#' ```{r}
#' in_dt <- data.table(intermittent=c(rep(0, 300), rep(1, 300)))
#' task = mlr3::TaskClassif$new(id = "task_classif", backend = in_dt, target ='intermittent')
#' convert_clastoregrtask(task, id="task_regr")
#' ```
#'
#' @export


convert_clastoregrtask <- function(in_task, in_id, oversample=FALSE) {
  if (oversample) {
    oversamp_ratio <- get_oversamp_ratio(in_task)
    mino_subdat <- in_task$data()[
      eval(in_task$target_names) == oversamp_ratio$minoclass,] #Get part of the underlying dataset of the minority class
    nmino <- mino_subdat[,.N]
    oversamp_index <- sample(nmino, nmino*(oversamp_ratio$ratio-1), replace=T)

    newdat <- rbind(in_task$data(),
                    mino_subdat[oversamp_index,]) %>%
      .[, eval(in_task$target_names):= as.numeric(as.character(get(in_task$target_names)))]
  } else {
    newdat <- in_task$data()[, eval(in_task$target_names) :=
                               as.numeric(as.character(get(in_task$target_names)))]
  }

  return(mlr3::TaskRegr$new(id =in_id,
                            backend = newdat,
                            target = in_task$target_names))
}

#Not used
#------ durfreq_parallel -----------------
#' Parallel wrapper for \link{comp_durfreq}
#'
#' Wrapper to run \link{comp_durfreq} (which computes the annual mean number of zero
#' flow days and frequency of zero-flow periods for a streamflow time series)
#' in parallel across a list of time series (paths).
#'
#' @param pathlist vector or list of paths to time series formatted GRDC-style
#' @inheritParams comp_durfreq
#' @param monthsel_list  named list with names being the gauge ID (e.g. GRDC_NO)
#' and values the selected months to compute the statistics over.
#' @param reverse (logical) whether to use monthsel_list as excluding months
#' rather than subsetting
#'
#' @return data.table with each row a gauging station and columns inherited
#' from \link{comp_durfreq}
#'
#' @export
durfreq_parallel <- function(pathlist, maxgap, monthsel_list=NULL,
                             reverse=FALSE) {

  cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)

  return(as.data.table(rbindlist(
    foreach(j=pathlist,
            .packages = c("data.table", "magrittr"),
            .export=c('comp_durfreq', 'zero_lomf', 'diny'))
    %dopar% {
      if (is.null(monthsel_list)) {
        return(comp_durfreq(j, maxgap = 20))
      } else {
        monthsublist = monthsel_list[[strsplit(basename(j), '[.]')[[1]][1]]]

        #Run duration and frequency computation for subset of months
        return(comp_durfreq(j, maxgap = maxgap,
                            monthsel= if (reverse) setdiff(1:12, monthsublist) else monthsublist)
        )
      }
    }
  )))
}

#------ rsmp_bbrier ----------------
rsmp_bbrier <- function(bmres, rsmp_i) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  if (inherits(rsmp_pred, 'PredictionClassif')) {
    bbrier <- sum((as.numeric(as.character(rsmp_pred$truth)) -
                     rsmp_pred$prob[,'1'])^2)/length(rsmp_pred$row_ids)
  }

  if (inherits(rsmp_pred, 'PredictionRegr')) {
    bbrier <- sum((as.numeric(as.character(rsmp_pred$truth)) -
                     rsmp_pred$response)^2)/length(rsmp_pred$row_ids)
  }
  return(bbrier)
}

#------ rsmp_auc ----------------
rsmp_auc <- function(bmres, rsmp_i) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred<- rsmp$prediction()
  auc <- Metrics::auc(rsmp_pred$truth, rsmp_pred$response)
  return(auc)
}

#------ rsmp_bacc ----------------
rsmp_bacc <- function(bmres, rsmp_i, threshold_class) {
  rsmp <- bmres$resample_result(rsmp_i)
  rsmp_pred <- rsmp$prediction()

  if (inherits(rsmp_pred, 'PredictionClassif')) {
    prob_resp <- rsmp_pred$prob[,'1']
  }

  if (inherits(rsmp_pred, 'PredictionRegr')) {
    prob_resp <- rsmp_pred$response
  }

  response <- as.data.table(prob_resp)[
    , c(fifelse(prob_resp>=threshold_class, '1', '0'))]

  bacc <- mlr3measures::bacc(factor(as.character(rsmp_pred$truth),
                                    levels=c('0', '1')),
                             factor(response, levels=c('0','1')))
  return(data.table(bacc=bacc, threshold_class=threshold_class))
}

#------ bm_paramstime ----------------
bm_paramstime <- function(bmres) {
  lrns_dt <- as.data.table(bmres$data)

  ##############################################################################
  params_format <- lapply(lrns_dt$learner, function(lrn) {
    print(lrn$id)
    #Get general parameters-----------------------------------------------------
    if (inherits(lrn ,'AutoTuner')) {
      lrn_in <- lrn$learner
    } else {
      lrn_in <- lrn
    }

    #If GraphLearner, get parameters for each Pipe Operator
    if (inherits(lrn_in, 'GraphLearner')) {
      paramtab <- lapply(lrn_in$graph$pipeops, function(po) {
        #Get default params
        def <- as.data.table(po$param_set$default) %>%
          .[, po_id := po$id] %>%
          melt(id.var='po_id')

        #Get set param
        value <- as.data.table(po$param_set$values) %>%
          .[, po_id := po$id]%>%
          melt(id.var='po_id')

        #Merge default and set params, keeping set ones when present
        if (ncol(def)==1) {
          tabmerge <- value
        } else if (ncol(value)==1) {
          tabmerge <- def
        } else {
          tabmerge <- merge(def, value[, .(variable, value, po_id)],
                            by=c('po_id', 'variable'),
                            all.x=T, all.y=T) %>%
            .[is.na(value.y), value.y := value.x] %>%
            setnames('value.y', 'value') %>%
            .[, value.x := NULL]
        }

        graphparams <- tabmerge %>%
          .[, variable := gsub('[.]', '', variable)] %>%
          unique(by='variable')

        return(graphparams)
      }) %>%
        rbindlist(use.names=TRUE)

      #If direct learner
    } else {
      #Get default params
      def <- as.data.table(lrn_in$param_set$default) %>%
        melt(id.vars=integer())

      #Get set param
      value <- as.data.table(lrn_in$param_set$values) %>%
        melt(id.vars=integer())

      if (def[,.N]>0 & value[,.N]>0) {
        tabmerge <- merge(def, value, by='variable', all.x=T, all.y=T)%>%
          .[is.na(value.y), value.y := value.x] %>%
          setnames('value.y', 'value') %>%
          .[, value.x := NULL]
      } else {
        tabmerge <- rbind(def, value)
      }

      #Merge default and set params, keeping set ones when present
      paramtab <- tabmerge %>%
        .[, variable := gsub('[.]', '', variable)] %>%
        unique(by='variable') #Models are not consistent in their use of dots
    }

    #Get inner sampling parameters----------------------------------------------
    if ('instance_args' %in% names(lrn)) {
      param_ranges <- as.list(
        paste0(lrn$instance_args$param_set$lower,
               '-',
               lrn$instance_args$param_set$upper))
      names(param_ranges) <- lrn$instance_args$param_set$ids()


      inner_smp  <- cbind(
        as.data.table(
          lrn$instance_args$resampling$param_set$get_values()),
        as.data.table(
          lrn$instance_args$terminator$param_set$get_values())
      ) %>%
        .[, nruns := prod(.SD)] %>%
        setnames(names(.), paste0('inner_', names(.))) %>%
        cbind(as.data.table(param_ranges))

      #If no inner sampling
    } else {
      inner_smp <- data.table(inner_nruns = 1)
    }


    #Remove pipe operator names and dot from names (inconsistent use across learners)
    setnames(inner_smp,
             old=names(inner_smp),
             new=gsub('[.]', '',
                      gsub(paste(paste0('(', unique(paramtab$po_id), '[.]*)'),
                                 collapse = '|'),
                           '',
                           names(inner_smp))))

    #---------------------------------------------------------------------------
    #Merge inner sampling parameters and general parameters, then format
    return(
      paramtab %>%
        .[!(variable %in% names(inner_smp)),] %>%
        dcast(.~variable, value.var='value') %>%
        cbind(inner_smp)
    )
  }) %>%
    rbindlist(., fill=T)
  ##############################################################################

  #Get task, learner and resampling ids
  lrns_time <- cbind(lrns_dt, params_format) %>%
    .[, learner_id := unlist(lapply(learner, function(x) x$id))] %>%
    .[, task_id := unlist(lapply(task, function(x) x$id))] %>%
    .[, resampling_id := unlist(lapply(resampling, function(x) x$id))]

  #Get number of features
  lrns_time[, npredictors :=
              lapply(lrns_dt$task,function(x) length(x$feature_names)) %>%
              unlist]

  #Compute train and predict time by learner, task and resampling
  bm_ttot <- as.data.table(bmres$aggregate(list(msr('time_train'),
                                                msr('time_predict')))) %>%
    .[lrns_time, on=c('learner_id', 'task_id', 'resampling_id')] %>%
    unique(by='nr') %>%
    .[, time_train := time_train/inner_nruns]

  return(bm_ttot)
}

#------ bm_msrtab ----------------
bm_msrtab <- function(bmres) {
  bbrier_vec <- lapply(1:bmres$n_resample_results, function(i) {
    rsmp_bbrier(bmres=bmres, rsmp_i=i)
  }) %>%
    unlist

  auc_vec <- lapply(1:bmres$n_resample_results, function(i) {
    rsmp_auc(bmres=bmres, rsmp_i=i)
  }) %>%
    unlist

  outer_smp <- lapply(bmres$resamplings$resampling, function(rsmp_design) {
    as.data.table(rsmp_design$param_set$get_values()) %>%
      setnames(names(.), paste0('outer_', names(.))) %>%
      .[, resampling_id := rsmp_design$id]
  }) %>%
    rbindlist(use.names=T)

  bacc_vec <- lapply(1:bmres$n_resample_results, function(i) {
    baccthresh <- ldply(seq(0,1,0.01), function(threshold_class) {
      rsmp_bacc(bmres, rsmp_i=i, threshold_class=threshold_class)}) %>%
      setorder(bacc) %>%
      setDT %>%
      .[, .SD[.N, .(bacc, threshold_class)]]
  }) %>%
    rbindlist

  #Get inner sampling params, general params, and time
  params_time_dt <- bm_paramstime(bmres)

  #Get time and add to bbrier
  moddt <- bm_paramstime(bmres) %>%
    cbind(., outer_smp) %>%
    cbind(., bacc_vec) %>%
    .[,`:=`(bbrier = bbrier_vec,
            auc = auc_vec)]

  return(moddt)
}

#------ format_modelcompdat --------------
format_modelcompdat <- function(bmres, typecomp=c('classif1', 'regr1', 'classif2')) {
  if (typecomp == 'classif1') {
    bmres[, `:=`(selection = 'Algorithm',
                 type = 'Classif.')]  %>%
      .[, learner_format := dplyr::case_when(
        learner_id == 'classif.ranger'~'default RF',
        learner_id == 'oversample.classif.ranger'~'default RF-oversampled',
        learner_id == 'classweights.classif.ranger'~'default RF-weighted classes',
        learner_id == 'classif.cforest'~'CIF',
        learner_id == 'oversample.classif.cforest'~'CIF-oversampled',
        learner_id == 'classweights.classif.cforest'~'CIF-weighted classes',
      )]

  } else if (typecomp == 'regr1') {
    bmres[, `:=`(selection = 'Algorithm',
                 type = 'Regr.')] %>%
      .[, learner_format := dplyr::case_when(
        task_id == 'inter_regr' ~'MAXSTAT',
        task_id == 'inter_regrover' ~'MAXSTAT-oversampled'
      )]

  } else if (typecomp == 'classif2') {
    bmres[, `:=`(selection = 'Predictors',
                 type = 'Classif.')] %>%
      .[, learner_format := dplyr::case_when(
        task_id == 'inter_basicsp' ~ 'default RF-oversampled-all variables',
        task_id == 'inter_basicsp_featsel' ~ 'default RF-oversampled-selected variables'
      )]
  } else {
    stop('typecomp is not recognized')
  }
}

#------ sfformat_wintri ----------------------
#Convert st to sf and transform to Winkel Trippel projection
sfformat_wintri <- function(in_sp) {
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"

  return(st_as_sf(in_sp) %>%
           st_transform(crs_wintri, use_gdal = FALSE)
  )
}



#------ bin_dt -----------------
#Bin a data_table on a given variable - see tabulate global summary
bin_dt <- function(in_dt, binvar, valuevar, binfunc, binarg,
                   bintrans=NULL, na.rm=FALSE) {
  #Inspired from rbin, adapted to dt and simplified
  in_dt <- copy(in_dt)

  el_freq <- function(byd, bins) {

    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    append(min(byd, na.rm = TRUE), min(byd, na.rm = TRUE) + (bin_length * seq_len(bins)))[1:bins]

  }

  eu_freq <- function(byd, bins) {

    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    ufreq      <- min(byd, na.rm = TRUE) + (bin_length * seq_len(bins))
    n          <- length(ufreq)
    ufreq[n]   <- max(byd, na.rm = TRUE) + 1
    return(ufreq)

  }

  binvar_orig <- copy(binvar)

  #Remove NAs
  if (na.rm) {
    in_dt <- in_dt[!is.na(get(eval(valuevar))),]
  }

  #Transform data if trans
  if (!is.null(bintrans)) {
    transvar <- paste0(binvar, '_bintrans')
    if (bintrans == 'log') {
      nneg <- in_dt[get(eval(binvar)) <= 0, .N]
      warning(paste0('There are ', nneg, ' records with', binvar, ' <= 0...',
                     'removing them for log transformation'))
      in_dt[, eval(transvar) :=
              log(get(eval(binvar)))]

    } else if (is.numeric(bintrans)) {
      in_dt[, eval(transvar) :=
              get(eval(binvar))^eval(bintrans)]

    }
    binvar = transvar
  }

  byd <- in_dt[, get(eval(binvar))]

  if (binfunc == 'manual') {
    l_freq    <- append(min(byd), binarg)
    u_freq    <- c(binarg, (max(byd, na.rm = TRUE) + 1))
    bins      <- length(binarg) + 1
  }

  if (binfunc == 'equal_length') {
    bins = round(binarg)
    l_freq    <- el_freq(byd, bins)
    u_freq    <- eu_freq(byd, bins)
  }

  if (binfunc == 'equal_freq') {
    bins = round(binarg)
    bin_prop     <- 1 / bins
    bin_length   <- in_dt[, round(.N/bins)]
    first_bins   <- (bins - 1) * bin_length
    residual     <- in_dt[, .N - first_bins]
    bin_rep      <- c(rep(seq_len((bins - 1)), each = bin_length),
                      rep(residual, residual))
    l_freq        <- c(1, (bin_length * seq_len((bins - 1)) + 1))
    u_freq       <- c(bin_length * seq_len((bins - 1)), in_dt[,.N])
    setorderv(in_dt, cols= eval(binvar))
    in_dt[, binid := .I]
    binvar = 'binid'
  }

  for (i in seq_len(bins)) {
    in_dt[get(eval(binvar)) >= l_freq[i] & get(eval(binvar)) < u_freq[i],
          bin := i]
    in_dt[bin == i, `:=`(bin_lmin = min(get(eval(binvar_orig)), na.rm=T),
                         bin_lmax = max(get(eval(binvar_orig)), na.rm=T))] %>%
      .[bin == i, bin_lformat := paste(bin_lmin, bin_lmax, sep='-')]

    if (i == bins) {
      in_dt[get(eval(binvar)) == u_freq[i],  bin := i]
    }
  }

  if (binfunc == 'equal_freq') {in_dt[, binid := NULL]}

  return(in_dt)
}

#------ label_manualbins ------------
label_manualbins <- function(binarg, minval) {
  minlabel <- paste(minval, binarg[1], sep=" - ")
  otherlabels <- mapply(function(x, y) {paste(x, y-1, sep=" - ")},
                        binarg[1:(length(binarg)-1)], binarg[2:length(binarg)])
  return(c(minlabel, otherlabels))
}

#------ formathistab -----------------
#Format summary table
formathistab <- function(in_dt, castvar, valuevar, valuevarsub,
                         weightvar, binfunc, binarg, binlabels,
                         datname) {
  rivbin <- bin_dt(in_dt = as.data.table(in_dt),
                   binvar = castvar,
                   valuevar = valuevar,
                   binfunc = binfunc,
                   binarg = binarg)

  netstat <- as.data.table(rivbin)[,sum(get(weightvar)),
                                   by=c(eval(valuevar), "bin")]
  tidyperc_riv <- netstat[, list(perc = 100*.SD[get(valuevar)==valuevarsub,
                                                sum(V1)]/sum(V1),
                                 binsumlength = sum(V1)),
                          by=c('bin')] %>%
    setorder(bin) %>%
    .[, `:=`(dat=datname,
             binformat = binlabels[bin])]
  return(tidyperc_riv)
}

#------ ggcompare -------------------
#Plot comparing intermittency prevalence and network length by drainage area for two datasets
ggcompare <- function(datmerge, binarg) {
  x_tick <- c(0, unique(datmerge$bin)) + 0.5
  binarg_tick <- c(0, binarg)
  len <- length(x_tick)

  plot_size <- ggplot(datmerge, aes(x=bin, y=binsumlength, fill=dat)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(name = 'Dataset', values=c('#a6cee3', '#1f78b4')) +
    scale_x_continuous(breaks = c(sort(unique(datmerge$bin)), x_tick)[seq(1, 2*len, 2)],
                       labels = c(rep(c(""), len), binarg)[seq(1, 2*len, 2)]) +
    labs(x= '', #bquote('Drainage area'~(km^2)),
         y='Total river length (km)') +
    coord_cartesian(expand=FALSE, clip="off") +
    theme_classic() +
    theme(legend.position = 'none',
          text = element_text(size=9),
          plot.background = element_blank(),
          axis.text.x = element_text(),
          axis.ticks.x = element_line(color = c(rep(NA, len/2 - 1),
                                                rep("black", len/2))))

  plot_inter <- ggplot(datmerge, aes(x=bin, y=perc, fill=dat)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(name = 'Dataset', values=c('#a6cee3', '#1f78b4')) +
    scale_x_continuous(breaks = c(sort(unique(datmerge$bin)), x_tick),
                       labels = c(rep(c(""), len), binarg)) +
    coord_cartesian(expand=FALSE, ylim=c(0, 1.4*max(datmerge$perc)), clip="off") +
    labs(x= bquote('Drainage area'~(km^2)),
         y='Prevalence of intermittency (% river length)') +
    theme_classic() +
    theme(legend.position = c(0.17, 0.90),
          legend.background = element_blank(),
          axis.ticks.x = element_line(color = c(rep(NA, len - 1), rep("black", len))))


  xmin_inset <- datmerge[, 0.4*(max(bin)-min(bin))]
  xmax_inset <- datmerge[, 1*max(bin)]
  ymin_inset <- datmerge[,round((0.8 * max(perc, na.rm=T) - min(perc, na.rm=T)) +
                                  min(perc, na.rm=T))]
  plot_join <- plot_inter +
    annotation_custom(grob = ggplotGrob(plot_size),
                      xmin = xmin_inset,
                      xmax = xmax_inset,
                      ymin = ymin_inset,
                      ymax = Inf)
  return(plot_join)
}

##### -------------------- Workflow functions ---------------------------------
#------ def_filestructure -----------------
#' Define file structure
#'
#' Defines file structure for running RF prediction of global intermittency
#'
#' @return Named list with paths to directories and files as well as paths for
#' analysis output
#'
#' @details all paths are relative to file in which project/package is contained. \cr
#' Project must be contained in ~/project_name/src \cr
#' Data must be contained in ~/project_name/data \cr
#' A directory called ~/project_name/results must exist to write outputs
#'
#' @export
#'
#'

def_filestructure <- function() {
  # Get main directory for project
  rootdir <- find_root(has_dir("src"))
  #Directory where raw data are located
  datdir <-  file.path(rootdir, 'data')
  # Directory where results and figures are written
  resdir <- file.path(rootdir, 'results')
  # File geodatabase to write outputs
  outgdb <- file.path(resdir, 'spatialoutputs.gdb')


  # hydrometric stations that have been joined to RiverATLAS
  in_gaugep <- file.path(outgdb, 'GRDCstations_riverjoin')
  # Directory containing hydrometric data
  in_gaugedir <-  file.path(datdir, 'GRDCdat_day')
  # River atlas formatted variables
  in_riveratlas_meta <- file.path(datdir, 'HydroATLAS', 'HydroATLAS_metadata_MLMv11.xlsx')
  #Average monthly discharge for HydroSHEDS network
  in_monthlynetdischarge <- file.path(datdir, 'HydroSHEDS',
                                      'HS_discharge_monthly.gdb',
                                      'Hydrosheds_discharge_monthly')
  #Rasters of dissolved buffers around gauge stations
  in_bufrasdir <- file.path(resdir, 'bufrasdir')
  # River atlas attribute data1
  in_riveratlas <- file.path(resdir, 'RiverATLAS_v10tab.csv')

  # River atlas attribute data2
  in_riveratlas2 <- file.path(resdir, 'RiverATLAS_v11tab.csv')

  # French river network for comparison
  compresdir <- file.path(resdir, 'Comparison_databases')
  in_netfr <- file.path(compresdir, 'france.gdb', 'network')
  in_basfr <- file.path(compresdir, 'france.gdb', 'hydrobasins12')

  # Output geopackage of hydrometric stations with appended predicted intermittency class
  out_gauge <- file.path(resdir, 'GRDCstations_predbasic800.gpkg')
  # River atlas predictions table
  out_riveratlas <- file.path(resdir, 'RiverATLAS_predbasic800.csv')

  return(c(rootdir=rootdir, datdir=datdir, resdir=resdir, outgdb=outgdb,
           in_gaugep=in_gaugep, in_gaugedir=in_gaugedir, in_riveratlas_meta=in_riveratlas_meta,
           in_monthlynetdischarge = in_monthlynetdischarge,
           in_bufrasdir = in_bufrasdir,
           in_netfr = in_netfr, in_basfr = in_basfr,
           out_gauge=out_gauge,
           in_riveratlas=in_riveratlas, in_riveratlas2=in_riveratlas2,
           out_riveratlas = out_riveratlas))
}

#------ read_monthlydis -------------
read_monthlydis <- function(in_filestructure) {
  monthlydischarge <- as.data.frame(st_read(
    dsn = dirname(in_filestructure['in_monthlynetdischarge']),
    layer = basename(in_filestructure['in_monthlynetdischarge'])
  )) %>%
    .[, c(paste0('DIS_', str_pad(seq(1, 12), 2, pad=0), '_CMS'),
          'REACH_ID')]
  return(monthlydischarge)
}

#------ read_gaugep -----------------
#' Read gauge points
#'
#' Import streamgauging station spatial data and attributes, only keeping those
#' whose point representation is within a given distance from the river network.
#'
#' @param in_filestructure named list containing path to point data for gauging
#' stations, named \code{in_gaugep}
#' @param dist maximum distance from the river network beyond which gauging
#' stations are excluded.
#'
#' @return object of class sf
#'
#' @details The distance from the river network must have been determined
#' beforehand and be an attribute of the gauge stations point data called
#' \code{station_river_distance}
#'
#' @export

read_gaugep <- function(in_filestructure, dist, in_monthlydischarge) {
  #Import gauge stations and only keep those < dist m from a HydroSHEDS reach

  gaugep <- st_read(dsn = dirname(in_filestructure['in_gaugep']),
                    layer = basename(in_filestructure['in_gaugep'])) %>%
    .[.$station_river_distance<dist,]

  riveratlas2 <- fread(in_filestructure['in_riveratlas2'])
  setnames(riveratlas2,
           names(riveratlas2),
           gsub('_11$', '', names(riveratlas2)))

  keepcols <- c('REACH_ID', 'Shape',
                names(gaugep)[!(names(gaugep) %in% names(riveratlas2))])
  gaugep_attriall <- merge(gaugep[,keepcols, with=F], riveratlas2,
                           by = 'REACH_ID')

  gaugep_monthlydischarge <- merge(gaugep_attriall, in_monthlydischarge,
                                   by.x = 'HYRIV_ID', by.y = 'REACH_ID',
                                   all.x=TRUE, all.y=FALSE)

  return(gaugep_monthlydischarge)
}

#------ read_gauged_paths -----------------
#' Read file paths to gauge flow data
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @inheritParams read_gaugep
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables
#'
#' @export

read_gauged_paths <- function(in_filestructure, in_gaugep) { #, gaugeid = 'GRDC_NO' down the line
  #Get data paths of daily records for gauge stations
  fileNames <- file.path(in_filestructure['in_gaugedir'], paste(in_gaugep$GRDC_NO,  ".txt", sep=""))
  #Check whether any GRDC record does not exist
  print(paste(length(which(do.call(rbind, lapply(fileNames, file.exists)))),
              'GRDC records do not exist...'))
  return(fileNames)
}
#------ comp_durfreq -------------------------------
#' Compute intermittency statistics for a gauging station
#'
#' Determine general characteristics of the whole time series and of the subset
#' of years that have less than a given threshold of missing data as well as
#' intermittency statistics. The intermittency statistics can be computed for a
#' subset of months of the year (e.g. only winter months)
#'
#' @param path file path to a GRDC-formatted streamflow time series table
#' @param maxgap maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param monthsel selected months to compute the statistics over
#'
#' @return One row data.table with the following columns: \cr
#' \describe{
#' \item{GRDC_NO} - (char) unique identifier for the gauge
#' \item{firstYear} - (num) first year on full record
#' \item{lastYear} - (num) last year on full record
#' \item{totalYears} - (int) total number of years on full record
#' \item{firstYear_kept} - (num) first year on record with < maxgap missing days
#' \item{lastYear_kept} - (num) first year on record with < maxgap missing days
#' \item{totalYears_kept} - (int) total number of years with < maxgap missing days
#' \item{totaldays} - (num) total number of days with discharge data
#' \item{sumDur} - (int) total number of days with discharge = 0
#' \item{mDur} - (num) mean number of days/year with discharge = 0
#' \item{mFreq} - (num) mean number of periods with discharge = 0
#' (with at least one day of flow between periods)
#' }
#'
#' @export

comp_durfreq <- function(path, maxgap, monthsel=NULL) {


  gaugetab <- cbind(fread(path, header=T, skip = 40, sep=";",
                          colClasses=c('character', 'character', 'numeric', 'numeric', 'integer')),
                    GRDC_NO = strsplit(basename(path), '[.]')[[1]][1])%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)

  gaugetab[, year:=as.numeric(substr(dates, 1, 4))] %>% #Create year column and total number of years on record
    .[, prevflowdate := gaugetab[zero_lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA]

  #Compute number of missing days per year
  gaugetab[!(Original == -999 | is.na(Original)), `:=`(missingdays = diny(year)-.N,
                                                       datadays = .N), by= 'year']


  gaugetab[,month:=as.numeric(substr(dates, 6, 7))] #create a month column

  if (gaugetab[(missingdays < maxgap), .N>0]) { #Make sure that there are years with sufficient data
    monthlyfreq <- gaugetab[
      (missingdays < maxgap),
      .(GRDC_NO = unique(GRDC_NO),
        monthrelfreq = length(unique(na.omit(prevflowdate)))/length(unique(year))), by='month'] %>% #Count proportion of years with zero flow occurrence per month
      dcast(GRDC_NO~month, value.var='monthrelfreq') %>%
      setnames(as.character(seq(1,12)),
               c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  } else {
    monthlyfreq <- data.table(GRDC_NO = gaugetab[, unique(GRDC_NO)], month=1:12, monthrelfreq=rep(NA,12)) %>%
      dcast(GRDC_NO~month, value.var='monthrelfreq') %>%
      setnames(as.character(seq(1,12)),
               c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  }

  #If analysis is only performed on a subset of months
  if (!is.null(monthsel)) {
    gaugetab <- gaugetab[month %in% monthsel, ]
  }
  #Compute number of days of zero flow for years with number of gap days under threshold
  gaugetab_yearly <- merge(gaugetab[, .(missingdays=max(missingdays, na.rm=T),
                                        datadays=max(datadays, na.rm=T)), by='year'],
                           gaugetab[Original == 0, .(dur=.N,
                                                     freq=length(unique(prevflowdate))), by='year'],
                           by = 'year', all.x = T) %>%
    .[!is.finite(missingdays), missingdays := diny(year)] %>%
    .[!is.finite(datadays), datadays := diny(year)] %>%
    .[dur==diny(year), freq:=1] %>%
    .[is.na(dur), `:=`(dur=0, freq=0)]


  gaugetab_all <- cbind(monthlyfreq,
                        gaugetab_yearly[, .(firstYear=min(year),
                                            lastYear=max(year),
                                            totalYears=length(unique(year)))],
                        gaugetab_yearly[missingdays <= maxgap,
                                        .(firstYear_kept=min(year),
                                          lastYear_kept=max(year),
                                          totalYears_kept=length(unique(year)),
                                          totaldays = sum(datadays),
                                          sumDur = sum(dur),
                                          mDur = mean(dur),
                                          mFreq = mean(freq),
                                          intermittent = factor(fifelse(mean(dur)>=1, 1, 0), levels=c('0','1')))]
  )
  return(gaugetab_all)
}
#------ format_gaugestats --------------------------------------------------------
#' Format gauge statistics
#'
#' Format gauge attributes, s
#'
#' @param path file path to a GRDC-formatted streamflow time series table
#' @param maxgap maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param monthsel selected months to compute the statistics over
#'
#' @return One row data.table with the following columns: \cr
#' \describe{
#' \item{GRDC_NO} - (char) unique identifier for the gauge
#' \item{firstYear} - (num) first year on full record
#' }
#'
#' @export

format_gaugestats <- function(in_gaugestats, in_gaugep) {
  #Format gaugestats into data.table

  #Join intermittency statistics to predictor variables and subset to only include those gauges with at least
  gaugestats_join <- do.call(rbind, in_gaugestats) %>%
    setDT %>%
    .[as.data.table(in_gaugep), on='GRDC_NO'] %>%
    .[!is.na(totalYears_kept) & totalYears_kept>=10,] %>% # Only keep stations with at least 10 years of data
    .[, c('X', 'Y') := as.data.table(sf::st_coordinates(geometry))] %>%
    comp_derivedvar #Compute derived variables, rescale some variables, remove -9999

  return(gaugestats_join)
}

#------ selectformat_predvars -----------------
selectformat_predvars <- function(in_filestructure, in_gaugestats) {
  #---- List predictor variables ----
  monthlydischarge_preds <- paste0('DIS_',
                                   str_pad(seq(1, 12), 2, pad=0),
                                   '_CMS')

  predcols<- c(
    monthlydischarge_preds,
    'UPLAND_SKM',
    'dis_m3_pyr',
    'dis_m3_pmn',
    'dis_m3_pmx',
    'dis_m3_pvar',
    'dis_m3_pvaryr',
    'run_mm_cyr',
    'runc_ix_cyr', #runoff coefficient (runoff/precipitation)
    'sdis_ms_uyr', #specific discharge
    'inu_pc_umn',
    'inu_pc_umx',
    'inu_pc_cmn',
    'lka_pc_cse',
    'lka_pc_use',
    'dor_pc_pva',
    'gwt_m_cav',
    'ele_pc_rel',
    'slp_dg_cav',
    'slp_dg_uav',
    'clz_cl_cmj',
    'pet_mm_uyr',
    'snw_pc_uyr',
    'snw_pc_cyr',
    'snw_pc_cmx',
    'glc_cl_cmj',
    'pnv_cl_cmj',
    'wet_pc_cg1',
    'wet_pc_cg2',
    'wet_pc_ug1',
    'wet_pc_ug2',
    'for_pc_use',
    'for_pc_cse',
    'ire_pc_use',
    'ire_pc_cse',
    'gla_pc_use',
    'gla_pc_cse',
    'prm_pc_use',
    'prm_pc_cse',
    'soc_th_uav',
    'soc_th_cav',
    'swc_pc_uyr',
    'swc_pc_cyr',
    'swc_pc_cmn',
    'lit_cl_cmj',
    'kar_pc_use',
    'kar_pc_cse',
    'ppd_pk_cav',
    'ppd_pk_uav',
    'urb_pc_cse',
    'urb_pc_use',
    'hft_ix_c93',
    'hft_ix_u93',
    'hft_ix_c09',
    'hft_ix_u09',
    'hdi_ix_cav',

    'cly_pc_cav',
    'cly_pc_uav',
    'slt_pc_cav',
    'slt_pc_uav',
    'snd_pc_cav',
    'snd_pc_uav',
    'pre_mm_c01',
    'pre_mm_c02',
    'pre_mm_c03',
    'pre_mm_c04',
    'pre_mm_c05',
    'pre_mm_c06',
    'pre_mm_c07',
    'pre_mm_c08',
    'pre_mm_c09',
    'pre_mm_c10',
    'pre_mm_c11',
    'pre_mm_c12',
    'pre_mm_u01',
    'pre_mm_u02',
    'pre_mm_u03',
    'pre_mm_u04',
    'pre_mm_u05',
    'pre_mm_u06',
    'pre_mm_u07',
    'pre_mm_u08',
    'pre_mm_u09',
    'pre_mm_u10',
    'pre_mm_u11',
    'pre_mm_u12',
    'pet_mm_cyr',
    'pet_mm_uyr',
    'cmi_ix_cyr',
    'cmi_ix_uyr',
    'cmi_ix_cmn',
    'cmi_ix_umn',
    'cmi_ix_cvar',
    'cmi_ix_uvar',
    'cmi_ix_c01',
    'cmi_ix_c02',
    'cmi_ix_c03',
    'cmi_ix_c04',
    'cmi_ix_c05',
    'cmi_ix_c06',
    'cmi_ix_c07',
    'cmi_ix_c08',
    'cmi_ix_c09',
    'cmi_ix_c10',
    'cmi_ix_c11',
    'cmi_ix_c12',
    'cmi_ix_u01',
    'cmi_ix_u02',
    'cmi_ix_u03',
    'cmi_ix_u04',
    'cmi_ix_u05',
    'cmi_ix_u06',
    'cmi_ix_u07',
    'cmi_ix_u08',
    'cmi_ix_u09',
    'cmi_ix_u10',
    'cmi_ix_u11',
    'cmi_ix_u12',
    'ari_ix_cav',
    'ari_ix_uav',
    'wloss_pc_cav',
    'wdryp_pc_cav',
    'wwetp_pc_cav',
    'whfrq_pc_cav',
    'wseas_pc_cav',
    'wperm_pc_cav',
    'wfresh_pc_cav',
    'wloss_pc_uav',
    'wdryp_pc_uav',
    'wwetp_pc_uav',
    'whfrq_pc_uav',
    'wseas_pc_uav',
    'wperm_pc_uav',
    'wfresh_pc_uav',
    'bio1_dc_cav',
    'bio2_dc_cav',
    'bio3_dc_cav',
    'bio4_dc_cav',
    'bio5_dc_cav',
    'bio6_dc_cav',
    'bio7_dc_cav',
    'bio8_dc_cav',
    'bio9_dc_cav',
    'bio10_dc_cav',
    'bio11_dc_cav',
    'bio12_mm_cav',
    'bio13_mm_cav',
    'bio14_mm_cav',
    'bio15_mm_cav',
    'bio16_mm_cav',
    'bio17_mm_cav',
    'bio18_mm_cav',
    'bio19_mm_cav',
    'bio1_dc_uav',
    'bio2_dc_uav',
    'bio3_dc_uav',
    'bio4_dc_uav',
    'bio5_dc_uav',
    'bio6_dc_uav',
    'bio7_dc_uav',
    'bio8_dc_uav',
    'bio9_dc_uav',
    'bio10_dc_uav',
    'bio11_dc_uav',
    'bio12_mm_uav',
    'bio13_mm_uav',
    'bio14_mm_uav',
    'bio15_mm_uav',
    'bio16_mm_uav',
    'bio17_mm_uav',
    'bio18_mm_uav',
    'bio19_mm_uav')

  #Check that all columns are in dt
  message(paste(length(predcols[!(predcols %in% names(in_gaugestats))]),
                'variables are missing from formatted gauge dataset'))

  #---- Associate HydroATLAS column names with variables names ----

  #Get predictor variable names
  metaall <- readxl::read_xlsx(in_filestructure['in_riveratlas_meta'],
                               sheet='Overall') %>%
    setDT

  metascale <- readxl::read_xlsx(in_filestructure['in_riveratlas_meta'],
                                 sheet='scale') %>%
    setDT %>%
    setnames(c('Key','Spatial representation'),
             c('Keyscale', 'Spatial.representation'))

  metastat <- readxl::read_xlsx(in_filestructure['in_riveratlas_meta'],
                                sheet='stat') %>%
    setDT %>%
    setnames(c('Key','Temporal or statistical aggregation or other association'),
             c('Keystat', 'Temporal.or.statistical.aggregation.or.other.association'))

  meta_format <- as.data.table(expand.grid(`Column(s)`=metaall$`Column(s)`,
                                           Keyscale=metascale$Keyscale,
                                           Keystat=metastat$Keystat)) %>%
    .[metaall, on='Column(s)'] %>%
    .[metascale, on = 'Keyscale'] %>%
    .[metastat, on = 'Keystat',
      allow.cartesian=TRUE]

  meta_format[, `:=`(
    unit = substr(`Column(s)`, 5, 6),
    varcode = paste0(gsub('[-]{3}', '', `Column(s)`),
                     Keyscale,
                     Keystat),
    varname = paste(Attribute,
                    Spatial.representation,
                    Temporal.or.statistical.aggregation.or.other.association))]

  #Add newly generated variables to meta_format (variable labels)
  addedvars <- data.table(varname=c('Precipitation catchment Annual min/max',
                                    'Discharge watershed Annual min/max',
                                    'Discharge watershed Annual min/average',
                                    'Elevation catchment average - watershed average',
                                    'Runoff coefficient catchment Annual average',
                                    'Specific discharge watershed Annual average',
                                    'Specific discharge watershed Annual min',
                                    paste0('Discharge watershed ', month.name)),
                          varcode=c('pre_mm_cvar',
                                    'dis_m3_pvar', 'dis_m3_pvaryr',
                                    'ele_pc_rel',
                                    'runc_ix_cyr',
                                    'sdis_ms_uyr', 'sdis_ms_umn',
                                    monthlydischarge_preds
                          )
  )

  oldcolnames <- c('Spatial.representation',
                   'Temporal.or.statistical.aggregation.or.other.association',
                   'Source Data')
  newcolnames <- c('Spatial representation',
                   'Temporal/Statistical aggreg.',
                   'Source')

  predcols_dt <- merge(data.table(varcode=predcols),
                       rbind(meta_format, addedvars, fill=T),
                       by='varcode', all.x=T, all.y=F)   %>%
    setnames(oldcolnames, newcolnames) %>%
    setorder(Category, Attribute,
             `Spatial representation`, `Temporal/Statistical aggreg.`)

  #Format table
  predcols_dt[varcode=='dis_m3_pvar', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='dis_m3_pvaryr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='runc_ix_cyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Runoff coefficient',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2, WorldClim v2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='sdis_ms_uyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='sdis_ms_umn', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[grepl('DIS_[0-9]{2}.*', varcode), `:=`(
    Category = 'Hydrology',
    Attribute= 'Nature Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`= gsub('[A-Z_]', '', varcode),
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]

  predcols_dt[varcode=='ele_pc_rel', `:=`(
    Category = 'Physiography',
    Attribute= 'Elevation',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='(cav-uav)/uav',
    Source = 'EarthEnv-DEM90',
    Citation = 'Robinson et al. 2014'
  )]

  return(predcols_dt)
}
#------ create_tasks  -----------------
create_tasks <- function(in_gaugestats, in_predvars) {
  #Create subset of gauge data for analysis (in this case, remove records with missing soil data)
  datsel <- in_gaugestats[!is.na(cly_pc_cav),
                          c('intermittent',in_predvars$varcode, 'X', 'Y'),
                          with=F]

  #Basic task for classification
  task_classif <- mlr3spatiotempcv::TaskClassifST$new(
    id="inter_basicsp",
    backend = datsel,
    target = "intermittent",
    coordinate_names = c("X", "Y"))

  #Basic task for regression without oversampling
  task_regr <- convert_clastoregrtask(in_task = task_classif,
                                      in_id = 'inter_regr',
                                      oversample=FALSE)

  #Basic task for regression with oversampling to have the same number of minority and majority class
  task_regrover <- convert_clastoregrtask(in_task = task_classif,
                                          in_id = 'inter_regrover',
                                          oversample=TRUE)
  return(list(classif=task_classif, regr=task_regr, regover=task_regrover))
}

#------ create_baselearners -----------------
create_baselearners <- function(in_task) {
  #---------- Create learners --------------------------------------------------
  lrns <- list()

  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }

  if (inherits(in_task, 'TaskClassif')) {
    #Compute ratio of intermittent to perennial observations
    imbalance_ratio <- get_oversamp_ratio(in_task)$ratio

    #Create basic learner
    lrns[['lrn_ranger']] <- mlr3::lrn('classif.ranger',
                                      num.trees = 500,
                                      sample.fraction = 0.632,
                                      replace = FALSE,
                                      splitrule = 'gini',
                                      predict_type = 'prob',
                                      importance = 'impurity_corrected',
                                      respect.unordered.factors = 'order')

    #print(lrn_ranger$param_set)

    #Create a conditional inference forest learner with default parameters
    # mtry = sqrt(nvar), fraction = 0.632
    lrns[['lrn_cforest']] <- mlr3::lrn('classif.cforest',
                                       ntree = 500,
                                       fraction = 0.632,
                                       replace = FALSE,
                                       alpha = 0.05,
                                       mtry = round(sqrt(length(in_task$feature_names))),
                                       predict_type = "prob")

    #Create mlr3 pipe operator to oversample minority class based on major/minor ratio
    #https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
    #https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
    #Sampling happens only during training phase.
    po_over <- mlr3pipelines::po("classbalancing", id = "oversample", adjust = "minor",
                                 reference = "minor", shuffle = TRUE,
                                 ratio = imbalance_ratio)
    #table(po_over$train(list(in_task))$output$truth()) #Make sure that oversampling worked

    #Create mlr3 pipe operator to put a higher class weight on minority class
    po_classweights <- mlr3pipelines::po("classweights", minor_weight = imbalance_ratio)

    #Create graph learners so that oversampling happens systematically upstream of all training
    lrns[['lrn_ranger_overp']] <- mlr3pipelines::GraphLearner$new(
      po_over %>>% lrns[['lrn_ranger']])
    lrns[['lrn_cforest_overp']] <- mlr3pipelines::GraphLearner$new(
      po_over %>>% lrns[['lrn_cforest']])

    #Create graph learners so that class weighin happens systematically upstream of all training
    lrns[['lrn_ranger_weight']] <- mlr3pipelines::GraphLearner$new(
      po_classweights %>>% lrns[['lrn_ranger']])
    lrns[['lrn_cforest_weight']] <- mlr3pipelines::GraphLearner$new(
      po_classweights  %>>% lrns[['lrn_cforest']])
  }


  if (inherits(in_task, 'TaskRegr')) {
    #Create regression learner with maxstat
    lrns[['lrn_ranger_maxstat']] <- mlr3::lrn('regr.ranger',
                                              num.trees=500,
                                              sample.fraction = 0.632,
                                              min.node.size = 10,
                                              replace=FALSE,
                                              splitrule = 'maxstat',
                                              importance = 'impurity_corrected',
                                              respect.unordered.factors = 'order')
  }

  return(lrns)
}

#------ set_tuning -----------------
set_tuning <- function(in_learner, in_measures, nfeatures,
                       insamp_nfolds, insamp_neval, insamp_nbatch) {

  if (is.list(in_learner)) {
    in_learner <- in_learner[[1]]
  }

  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)

    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T),
                   lower = 5,
                   upper = floor(nfeatures/2)), #Half number of features
      ParamDbl$new(grep(".*fraction", prmset, value=T),
                   lower = 0.2,
                   upper = 0.8)
    ))

    in_split =in_lrn$param_set$get_values()[
      grep(".*split(rule|stat)", prmset, value=T)]

    if (in_split == 'maxstat') {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)],
                     lower = 0.01, upper = 0.1)
      )

    } else if (any(grepl(".*min.node.size", prmset))) {
      tune_rf$add(
        ParamInt$new(prmset[grep(".*min.node.size", prmset)],
                     lower = 1, upper = 10)
      )
    } else if (any(grepl(".*splitstat", prmset))) {
      tune_rf$add(
        ParamDbl$new(prmset[grep(".*alpha", prmset)],
                     lower = 0.01, upper = 0.1)
      )
    }
  }

  #Define inner resampling strategy
  rcv_rf = rsmp("cv", folds=insamp_nfolds) #aspatial CV repeated 10 times

  #Define termination rule
  evalsn = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds

  if (in_learner$task_type == 'classif') {
    if (grepl('classif[.]cforest$', in_learner$id)) {
      learnertune <- in_learner
    } else if (grepl('classif[.]ranger$', in_learner$id)) {
      learnertune <- AutoTuner$new(learner= in_learner,
                                   resampling = rcv_rf,
                                   measures = in_measures$classif,
                                   tune_ps = regex_tuneset(in_learner),
                                   terminator = evalsn,
                                   tuner =  tnr("random_search",
                                                batch_size = insamp_nbatch)) #batch_size determines level of parallelism
    } else{
      stop('The classification learner provided is not configurable with this workflow yet...')
    }
  } else if (in_learner$task_type == 'regr') {
    learnertune <- AutoTuner$new(learner= in_learner,
                                 resampling = rcv_rf,
                                 measures = in_measures$regr,
                                 tune_ps = regex_tuneset(in_learner),
                                 terminator = evalsn,
                                 tuner =  tnr("random_search",
                                              batch_size = insamp_nbatch))
  }

  learnertune$id <- in_learner$id

  return(learnertune)
}

#------ instantiate resampling ---------------
set_cvresampling <- function(rsmp_id, in_task, outsamp_nrep, outsamp_nfolds) {
  #repeated_cv or repeated-spcv-coords
  outer_resampling = rsmp(rsmp_id,
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)
  outer_resampling$instantiate(in_task)

  return(outer_resampling)
}

#------ dynamic_resample ------------------
#Run resample on in_task, selected learner (in_lrnid) from in_bm, in_resampling
dynamic_resample <- function(in_task, in_learner, in_resampling, type,
                             store_models = TRUE) {
  if (is.list(in_learner)) {
    in_learner <- in_learner[[1]]
  }

  if (is.list(in_task)) {
    in_task <- in_task[[1]]
  }

  if (inherits(in_learner, 'BenchmarkResult')) {
    print(('BenchmarkResults was provided, getting the learner...'))
    in_learner <- in_learner$learners$learner[[1]]
  }

  if ((in_learner$task_type == 'classif' & type=='classif') |
      (in_learner$task_type == 'regr' & type=='regr')) {
    resmp_rs <- mlr3::resample(learner = in_learner, task = in_task,
                               resampling = in_resampling, store_models = store_models
    )
    return(resmp_rs)
  }
}

#------ dynamic_resamplebm ------------------
#Run resample on in_task, selected learner (in_lrnid) from in_bm, in_resampling
dynamic_resamplebm <- function(in_task, in_bm, in_lrnid, in_resampling, type,
                               store_models = TRUE) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  #get desired resampled_results/learner
  in_rf <- in_bm$filter(learner_ids = in_lrnid)

  return(
    dynamic_resample(in_task = in_task,
                     in_learner = in_rf,
                     in_resampling = in_resampling,
                     type = type,
                     store_models = store_models)
  )
}

#------ combined resample results into benchmark results -------------
combine_bm <- function(in_resampleresults, out_qs) {
  #When tried as_benchmark_result.ResampleResult, got "Error in setcolorder(data, slots) :
  # x has some duplicated column name(s): uhash. Please remove or rename the
  # duplicate(s) and try again.". SO use this instead
  print('Converting to benchmark results...')
  if (length(in_resampleresults) > 1) {
    bmres_list <- lapply(
      in_resampleresults[!sapply(in_resampleresults, is.null)],
      function(rsmpres) {
        print(rsmpres)
        if (!is.null(rsmpres)) {
          as_benchmark_result(rsmpres)
        }
      })
    #BenchmarkResult$new(rsmpres$data)})

    print('Combining...')
    bmrbase = bmres_list[[1]]
    for (i in 2:length(bmres_list)) {
      if (in_resampleresults[[i]]$task$task_type ==
          in_resampleresults[[1]]$task$task_type) {
        print(i)
        bmrbase$combine(bmres_list[[i]])
      } else {
        warning('ResampleResult #', i,
                'is not of the same task type as the first ResampleResult you provided, skipping...')
      }
    }
  } else {
    warning('You provided only one resample result to combine_bm,
            simply returning output from as_benchmark_result...')
    bmrbase = BenchmarkResult$new(in_resampleresults[[1]])
  }
  print('Done combining, now writing to qs...')
  qs::qsave(bmrbase, out_qs)
}

#------ analyze_benchmark -----------------
analyze_benchmark <- function(in_bm, in_measure) {

  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  print(paste('It took',
              in_bm$aggregate(msr('time_both'))$time_both,
              'seconds to train and predict with the',
              in_bm$aggregate(msr('time_both'))$learner_id,
              'model...'))

  bmdt <- as.data.table(in_bm)

  if (in_bm$task_type == 'regr') {
    print(in_bm$aggregate(in_measure$regr))
    boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure$regr)

    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- bmdt$prediction[[rsmp_i]]$test %>%
        as.data.table %>%
        .[, `:=`(outf = bmdt$iteration[[rsmp_i]],
                 task = bmdt$task[[rsmp_i]]$id,
                 task_type = in_bm$task_type,
                 learner = bmdt$learner[[rsmp_i]]$id)]
      return(preds)
    }) %>%
      do.call(rbind, .)

    if (!('prob.1' %in% names(preds)) & 'response' %in% names(preds)) {
      preds[, prob.1 := response]
    }
  }


  if (in_bm$task_type == 'classif') {
    print(in_bm$aggregate(in_measure$classif))
    boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure$classif)

    preds <- lapply(seq_len(bmdt[,.N]), function(rsmp_i) {
      preds <- data.table(outf = bmdt$iteration[[rsmp_i]],
                          task = bmdt$task[[rsmp_i]]$id,
                          learner = bmdt$learner[[rsmp_i]]$id,
                          task_type = in_bm$task_type,
                          pred = list(bmdt$prediction[[rsmp_i]]$test))
      return(preds)
    }) %>%
      do.call(rbind, .)
  }

  tasklearner_unique <- preds[, expand.grid(unique(task), unique(learner))] %>%
    setnames(c('task', 'learner')) %>%
    setDT

  tasklearner_unique[, learner_format := dplyr::case_when(
    learner == 'classif.ranger.tuned'~'default RF',
    learner == 'oversample.classif.ranger.tuned'~'default RF - oversampled',
    learner == 'classweights.classif.ranger.tuned'~'default RF - weighted classes',
    learner == 'classif.cforest'~'CIF',
    learner == 'oversample.classif.cforest'~'CIF - oversampled',
    learner == 'classweights.classif.cforest'~'CIF - weighted classes',
  )]

  glist <- lapply(1:nrow(tasklearner_unique), function(tsklrn) {
    print(tasklearner_unique[tsklrn,])
    subpred <- preds[task ==tasklearner_unique$task[tsklrn] &
                       learner == tasklearner_unique$learner[tsklrn],]

    gout <- ggmisclass_single(in_predictions = subpred) +
      ggtitle(paste(tasklearner_unique$task[tsklrn],
                    tasklearner_unique$learner_format[tsklrn])) +
      labs(x='Threshold', y='Value')

    if (tsklrn < nrow(tasklearner_unique)) {
      gout <- gout +
        theme(legend.position = 'none')
    }

    return(ggplotGrob(gout))
  })

  return(list(bm_misclasscomp=do.call("grid.arrange", list(grobs=glist)),
              bm_boxcomp = boxcomp))
}
#------ select_features ------------------
select_features <- function(in_bm, in_lrnid, in_task, pcutoff) {

  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  #get desired resampled_results/learner
  in_rf <- in_bm$filter(learner_ids = in_lrnid)

  #Apply feature/variable selection
  vimp <- weighted_vimportance_nestedrf(
    rfresamp = in_rf$resample_result(uhash=unique(as.data.table(in_rf)$uhash)),
    pvalue = TRUE) %>%
    .[,imp_wmeanper := imp_wmean/sum(imp_wmean)]

  task_featsel <- in_task$clone()$select(
    vimp[imp_pvalue <= pcutoff, as.character(varnames)])
  task_featsel$id <- paste0(in_task$id, '_featsel')

  return(list(in_task, task_featsel))
}

#------ selecttrain_rf -----------------
selecttrain_rf <- function(in_rf, in_learnerid, in_taskid,
                           insamp_nfolds =  NULL, insamp_nevals = NULL) {
  #Prepare autotuner for full training
  in_bmsel <- in_rf$clone()$filter(learner_ids = in_learnerid,
                                   task_id = in_taskid)

  lrn_autotuner <- in_bmsel$clone()$learners$learner[[1]]
  in_task <-in_bmsel$tasks$task[[1]]

  if (!is.null(insamp_nfolds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- insamp_nfolds
  }

  if (!is.null(insamp_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- insamp_nevals
  }

  #Return outer sampling object for selected model (or list of outer sampling objects)
  uhashes <- unique(as.data.table(in_bmsel)$uhash)
  if (length(uhashes) == 1) {
    outer_resampling_output <- in_bmsel$resample_result(uhash=uhashes)
  } else {
    outer_resampling_output <- lapply(uhashes, function(x) {
      in_bmsel$resample_result(uhash=x)
    })
  }

  #Train learners
  lrn_autotuner$param_set$values = mlr3misc::insert_named(
    lrn_autotuner$param_set$values,
    list(classif.ranger.importance = 'permutation')
  )
  lrn_autotuner$train(in_task)

  return(list(rf_outer = outer_resampling_output, #Resampling results
              rf_inner = lrn_autotuner, #Core learner (with hyperparameter tuning)
              task = in_task)) #Task
}


#------ rformat_network ------------------
rformat_network <- function(in_filestructure, in_predvars, in_monthlydischarge) {
  cols_tokeep <-  c("HYRIV_ID", "HYBAS_L12", "UPLAND_SKM", "LENGTH_KM",
                    in_predvars[!is.na(ID), varcode],
                    'ele_mt_cav','ele_mt_uav', 'gwt_cm_cav', 'ORD_STRA',
                    'hft_ix_c09', 'hft_ix_u09',
                    paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
                    paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
                    paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0)))

  riveratlas <- fread_cols(file_name=in_filestructure['in_riveratlas'],
                           cols_tokeep = cols_tokeep) %>%
    merge(as.data.table(in_monthlydischarge),
          by.x='HYRIV_ID', by.y='REACH_ID') %>%
    comp_derivedvar
}

#------ write_preds -----------------
write_preds <- function(in_filestructure, in_gaugep, in_gaugestats,
                        in_network, in_rftuned, in_predvars, in_uncertainty) {
  # ---- Output GRDC predictions as points ----
  in_gaugestatsformat <- in_gaugestats[
    !is.na(cly_pc_cav),  #SHould change that to not have to adjust subselection in multiple spots
    IRpredprob := in_rftuned$rf_inner$predict(in_rftuned$task)$prob[,2]] %>%
    .[in_uncertainty, on='GRDC_NO']

  cols_toditch<- colnames(in_gaugep)[colnames(in_gaugep) != 'GRDC_NO']


  out_gaugep <- merge(in_gaugep,
                      in_gaugestatsformat[, !cols_toditch, with=F],
                      by='GRDC_NO')

  st_write(obj=out_gaugep,
           dsn=in_filestructure['out_gauge'],
           driver = 'gpkg',
           delete_dsn=T)

  #Predict model (should take ~10-15 minutes for 9M rows x 40 cols — chunk it up by climate zone to avoid memory errors)
  for (clz in unique(in_network$clz_cl_cmj)) {
    print(clz)
    tic()
    in_network[!is.na(cly_pc_cav) & !is.na(cly_pc_uav) & clz_cl_cmj == clz,
               predbasic800 := in_rftuned$rf_inner$predict_newdata(as.data.frame(.SD))$prob[,2],]
    toc()
  }

  in_rftuned$rf_inner$predict_newdata(as.data.frame(
    in_network[!is.na(cly_pc_cav) & !is.na(cly_pc_uav) & clz_cl_cmj == clz,]
  ))

  in_network[, predbasic800cat := ifelse(predbasic800>=0.5, 1, 0)]
  fwrite(in_network[, c('HYRIV_ID', 'HYBAS_L12', 'predbasic800', 'predbasic800cat'), with=F],
         in_filestructure['out_riveratlas'])

  # --------- Return data for plotting ------------------------
  return(out_gaugep)
}


##### -------------------- Diagnostics functions -------------------------------

#------ ggvimp -----------------
ggvimp <- function(in_rftuned, in_predvars, varnum = 10, spatial_rsp=FALSE) {
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)

  #Get variable importance and format them
  varimp_basic <- weighted_vimportance_nestedrf(rfresamp = rsmp_res,
                                                pvalue = FALSE) %>%
    merge(., in_predvars, by.x='varnames', by.y='varcode') %>%
    .[, varname := factor(varname, levels=varname[order(-imp_wmean)])]

  #Plot 'em
  ggplot(varimp_basic[1:varnum,],aes(x=varname)) +
    geom_bar(aes(y=imp_wmean), stat = 'identity') +
    geom_errorbar(aes(ymin=imp_wmean-2*imp_wsd, ymax=imp_wmean+2*imp_wsd)) +
    scale_x_discrete(name='Predictor',
                     labels = function(x) stringr::str_wrap(x, width = 10)) +
    scale_y_continuous(name='Variable importance', expand=c(0,0)) +
    coord_cartesian(ylim=c(0,100)) +
    theme_classic() +
    theme(axis.text.x = element_text(size=8))
}

#------ ggmisclass_single -----------------
ggmisclass_single <-  function(in_predictions=NULL, in_rftuned=NULL, spatial_rsp=FALSE) {
  #Get predicted probabilities of intermittency for each gauge
  # in_gaugestats[!is.na(cly_pc_cav), intermittent_predprob :=
  #                 as.data.table(in_predictions)[order(row_id), mean(prob.1), by=row_id]$V1]
  #Get misclassification error, sensitivity, and specificity for different classification thresholds
  #i.e. binary predictive assignment of gauges to either perennial or intermittent class
  if (!is.null(in_rftuned)) {
    rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)
    in_predictions <- rsmp_res$prediction()
  }

  threshold_confu_dt <- ldply(seq(0,1,0.01), threshold_misclass, in_predictions) %>%
    setDT

  #Get classification threshold at which sensitivity and specificity are the most similar
  balanced_thresh <- threshold_confu_dt[which.min(abs(spec-sens)),]
  print(paste('Sensitivity =', round(balanced_thresh$sens,2),
              'and Specificity =', round(balanced_thresh$spec,2),
              'at a classification threshold of', balanced_thresh$i))

  gout <- ggplot(melt(threshold_confu_dt, id.vars='i'),
                 aes(x=i, y=value, color=variable, linetype=variable)) +
    geom_line(size=1.2) +
    geom_vline(xintercept=balanced_thresh$i, alpha=1/2) +
    geom_hline(yintercept=balanced_thresh$spec, alpha=1/2) +
    annotate('text', x=(balanced_thresh$i), y=0.4,
             label=balanced_thresh$i, angle=-90) +
    annotate('text', x=0.9, y=(balanced_thresh$spec),
             label=round(balanced_thresh$sens,2)) +
    scale_x_continuous(expand=c(0,0), name='Threshold') +
    scale_y_continuous(expand=c(0,0), name='Value') +
    scale_color_brewer(palette='Dark2',  #colorblind friendly
                       labels=c('Misclassification rate',
                                'Sensitivity (true positives)',
                                'Specificity (true negatives)')) +
    theme_bw()

  #Plot it
  return(gout)
}

#------ ggpd_bivariate -----------------
ggpd_bivariate <- function (in_rftuned, in_predvars, colnums, ngrid, nodupli=T,
                            parallel=T, spatial_rsp=FALSE) {

  #Get outer resampling of interest
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)

  #Get partial dependence across all folds
  nlearners <- rsmp_res$resampling$param_set$values$folds
  datdf <- as.data.frame(rsmp_res$task$data()) #This may be shortened
  varimp <- weighted_vimportance_nestedrf(rsmp_res, pvalue=FALSE)

  if (nodupli) {
    selcols <- as.character(
      varimp$varnames[!duplicated(substr(varimp$varnames, 1,3))][colnums])
  } else {
    selcols <- as.character(
      varimp$varnames[colnums])
  }

  if (parallel) {
    print(paste("Computing partial dependence with future.apply across", nlearners,
                "CV folds"))
    pd <- future.apply::future_lapply(seq_len(nlearners),
                                      extract_pd_nestedrf,
                                      in_rftuned = rsmp_res,
                                      datdf = datdf,
                                      selcols = selcols,
                                      nvariate=2,
                                      ngrid = ngrid,
                                      future.scheduling = structure(TRUE,ordering = "random"),
                                      future.packages = c("data.table","edarf","ranger"))

  } else {
    print(paste("Computing partial dependence iteratively across", nlearners,
                "CV folds"))
    pd <- lapply(seq_len(nlearners),
                 extract_pd_nestedrf,
                 in_rftuned = rsmp_res,
                 datdf = datdf,
                 selcols = selcols,
                 ngrid = ngrid)
  }

  #Get weighted mean
  pdformat <- do.call(rbind, pd) %>%
    setDT %>%
    .[, list(mean1 = weighted.mean(`1`, classif.bacc)),
      by= c('var1', 'var2', 'value1', 'value2')] %>%
    .[, variables := paste(var1, var2)]

  vargrid <- t(combn(1:length(selcols), 2))
  #leglims <- pdformat[, c(min(mean1), max(mean1))]

  #Iterate over every pair of variables

  tileplots_l <- pdformat[,list(list(ggplotGrob(
    ggplot(.SD, aes(x=value1, y=value2)) +
      geom_tile(aes(fill = mean1)) +
      scale_fill_distiller(palette='Grey') +
      geom_jitter(data=datdf,
                  aes_string(color='intermittent', x=eval(var1),y=eval(var2)),
                  alpha=1/3) +
      scale_color_manual(values=c('#0F9FD6','#ff9b52')) +
      labs(x=stringr::str_wrap(in_predvars[varcode==eval(var1), varname],
                               width = 20),
           y=stringr::str_wrap(in_predvars[varcode==eval(var2), varname],
                               width = 20)) +
      theme_bw() +
      theme(text = element_text(size=12))
  )))
  , by=.(var1, var2)]



  pagelayout <-   lapply(1:(nrow(tileplots_l) %/% 9), function(p_i) {
    (p_i-1)*9+(1:9)
  })
  if (nrow(tileplots_l) %% 9 > 0) {
    pagelayout[[nrow(tileplots_l) %/% 9 + 1]] <- (p_i-1)*9+(1:(nrow(tileplots_l) %% 9))
  }

  tileplots_multipl <- lapply(pagelayout, function(page) {
    print(page)
    return(do.call("grid.arrange", list(grobs=(tileplots_l[page,V1]))))
  })
  return(tileplots_multipl)
}

#------ gguncertainty -----------------
gguncertainty <- function(in_rftuned, in_gaugestats, in_predvars, spatial_rsp) {
  #Get outer resampling of interest
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=spatial_rsp)

  #Get average predictions for oversampled rows
  gaugepred <-  rsmp_res$prediction() %>%
    as.data.table %>%
    .[, list(truth=first(truth), prob.1=mean(prob.1)), by=row_id] %>%
    setorder(row_id)

  predattri <- cbind(in_gaugestats[!is.na(cly_pc_cav),], gaugepred) %>%
    .[, `:=`(preduncert = prob.1-as.numeric(as.character(intermittent)),
             yearskeptratio = totalYears_kept/totalYears)]

  #Plot numeric variables
  predmelt_num <- predattri[, which(as.vector(unlist(lapply(predattri, is.numeric)))), with=F] %>%
    cbind(predattri[, c('GRDC_NO', 'intermittent'), with=F]) %>%
    melt(id.vars=c('GRDC_NO', 'intermittent', 'prob.1', 'preduncert'))


  #Set variable labels
  varlabels <- copy(in_predvars) %>%
    .[, .(varcode, varname)] %>%
    setkey(varcode) %>%
    .[levels(predmelt_num$variable)] %>%
    .[is.na(varname), varname:=varcode] %>%
    .[varcode=='totalYears_kept', varname := 'Years of record kept'] %>%
    .[varcode=='yearskeptratio', varname :='Years of record kept/All years'] %>%
    .[varcode=='mDur', varname :='Mean annual # of dry days'] %>%
    .[varcode=='mFreq', varname :='Mean annual # of dry periods'] %>%
    .[varcode=='station_river_distance', varname :='Station distance to network (m)'] %>%
    .[varcode=='ORD_STRA', varname :='Strahler river order'] %>%
    .[varcode=='UPLAND_SKM', varname :='Drainage area (km2)']


  levels(predmelt_num$variable) <-  varlabels$varname
  varstoplot <- varlabels[varcode %in% c('totalYears_kept', 'yearskeptratio',
                                     'mDur', 'mFreq', 'station_river_distance',
                                     'UPLAND_SKM', 'ORD_STRA',
                                     'dis_m3_pyr', 'dor_pc_pva',
                                     'cmi_ix_uyr','ari_ix_uav'),]

  plotdt <- predmelt_num[variable %in% varstoplot$varcode |
                           variable %in% varstoplot$varname,]

  colorpal <- c('#1f78b4', '#ff7f00')
  rectdf <- data.table(
    xmin=rep(-Inf, 4),
    xmax=rep(Inf, 4),
    ymin=c(-1, -0.5, 0, 0.5),
    ymax=c(-0.5, 0, 0.5, 1),
    fillpal = rep(colorpal, 2)
  )
  uncertainty_numplot <-
    ggplot(plotdt) +
      geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
                            ymin=ymin, ymax=ymax, fill=fillpal),
           alpha=1/4) +
      scale_fill_manual(values=colorpal,
                        name='Predicted regime',
                        labels = c('Perennial', 'Intermittent')) +
      geom_point(aes(x=value, y=preduncert, color=intermittent), alpha = 1/4) +
      geom_hline(yintercept=0, alpha=1/2) +
      new_scale_fill() +
      geom_smooth(aes(x=value, y=preduncert, color=intermittent),
                  method='gam', formula = y ~ s(x, k=3)) +
      # annotate("text", x = Inf-5, y = 0.5, angle = 90,
      #          label = "Pred:Int, Obs:Per",
      #          color = colorpal[[1]]) +
      # annotate("text", x = Inf-5, y = -0.5, angle = 90,
      #          label = "Pred:Per, Obs:Int",
      #          color = colorpal[[2]]) +
      scale_color_manual(values=colorpal,
                         name='Observed regime',
                         labels = c('Perennial', 'Intermittent')) +
      labs(x='Value', y='Intermittency Prediction Residuals (IPR)') +
      #scale_x_sqrt(expand=c(0,0)) +
      coord_cartesian(expand=FALSE, clip='off') +
      facet_wrap(~variable, scales='free') +
      theme_classic() +
      theme(legend.position = c(0.85, 0.1))


    #Plot categorical variables
    predmelt_cat <- predattri[, c('GRDC_NO', 'intermittent', 'preduncert',
                                  'ENDORHEIC', 'clz_cl_cmj'), with=F] %>%
      melt(id.vars=c('GRDC_NO', 'intermittent', 'preduncert'))
    levels(predmelt_cat$variable) <- c('Endorheic',
                                       'Climate Zone (catchment majority)')

    uncertainty_catplot <-
      ggplot(predmelt_cat) +
      geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
                                 ymin=ymin, ymax=ymax, fill=fillpal),
                alpha=1/4) +
      scale_fill_manual(values=colorpal,
                        name='Predicted regime',
                        labels = c('Perennial', 'Intermittent')) +
      #geom_boxplot(alpha = 0.75) +
      new_scale_fill() +
      geom_violin(aes(x=as.factor(value), y=preduncert,
                      fill=intermittent, color=intermittent),
                  alpha=0.75, color=NA) +
      geom_hline(yintercept=0, alpha=1/2) +
      labs(x='Value', y='Intermittency Prediction Residuals (IPR)') +
      coord_cartesian(expand=FALSE, clip='off') +
      scale_fill_manual(values=colorpal,
                        name='Observed regime',
                        labels = c('Perennial', 'Intermittent')) +
      scale_color_manual(values=c('#175885', '#9e3f00'),
                         name='Observed regime',
                         labels = c('Perennial', 'Intermittent')) +
      facet_wrap(~variable, scales='free', labeller=label_value) +
      theme_classic()

    return(list(uncertainty_numplot=uncertainty_numplot,
              uncertainty_catplot=uncertainty_catplot,
              out_uncertainty = predattri[, .(GRDC_NO, preduncert)]))
}

#------ krige_spuncertainty----
krige_spuncertainty <- function(in_filestructure, in_rftuned,
                                in_gaugep, in_gaugestats, kcutoff=50000,
                                overwrite = FALSE) {
  rsmp_res <- get_outerrsmp(in_rftuned, spatial_rsp=TRUE)
  predsp <- rsmp_res$prediction() %>%
    as.data.table %>%
    .[, list(IRpredprob_spcv = 100*mean(prob.1),
             IRpredcoefcv_spcv = sd(100*prob.1, na.rm=T)/mean(100*prob.1, na.rm=T)),
      by=.(row_id, truth)] %>%
    .[, prederror := IRpredprob_spcv-100*as.numeric(as.character(truth))] %>%
    .[, prederror_abs := abs(prederror)] %>%
    .[, IRpedcat_spcv := as.character(fifelse(IRpredprob_spcv > 40, 1, 0))] %>%
    setorder(row_id) %>%
    cbind(in_gaugestats[!is.na(cly_pc_cav), list(GRDC_NO=GRDC_NO)])

  predsp_gaugep <- merge(in_gaugep, predsp, by='GRDC_NO') %>%
    .[order(.[,'prederror_abs']$prederror_abs),]

  bufras_vec <- file.path(in_filestructure[['in_bufrasdir']],
                          grep('bufras_T.*proj[.]tif$',
                               list.files(in_filestructure[['in_bufrasdir']]),
                               value=T)
  )

  #Convert to SpatialPointDataFrame to use in gstats and project to Goode Homolosine
  crs_aeqd <- "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  predsp_gaugep_df <- st_transform(predsp_gaugep, crs=crs_aeqd)  %>%
    as_Spatial()

  # Compute the sample variogram; note that the f.1 trend model is one of the
  # parameters passed to variogram(). This tells the function to create the
  # variogram on the de-trended data.
  var_smpl <- gstat::variogram(prederror ~ 1, predsp_gaugep_df,
                               cutoff=kcutoff, width=500, cressie=TRUE) #cloud=T)

  ggplot(var_smpl, aes(x=dist, y=gamma)) +
    geom_point(alpha=1/3) +
    coord_cartesian(expand=F) +
    scale_y_log10() +
    theme_classic()

  # scale_x_sqrt(breaks=c(0,100,1000,10000,25000,50000, 100000),
  #              labels=c(0,100,1000,10000,25000,50000, 100000)) +

  # Compute the variogram model by passing the nugget, sill and range values
  # to fit.variogram() via the vgm() function.
  varmods <- as.character(vgm()[,'short'])
  dat_fit  <- fit.variogram(var_smpl, fit.ranges = TRUE, fit.sills = FALSE,
                            vgm(varmods[!(varmods %in% c('Pow', 'Int'))]),
                            fit.kappa = TRUE)
  # The following plot allows us to assess the fit
  plot(var_smpl, dat_fit)

  # Import gauge buffer mask


  #Predict error
  kmod <- gstat(formula=prederror~1, locations=predsp_gaugep_df, model=dat_fit,
                maxdist = kcutoff)
  kpl <-  lapply(seq_along(bufras_vec), function(i) {
    bufmask <- raster(bufras_vec[[i]]) %>%
      as('SpatialGrid')
    kp <- round(raster(predict(kmod, bufmask)))

    outras = file.path(dirname(bufras_vec[[i]]),
                       gsub('bufras', 'krigpred', basename(bufras_vec[[i]])))


    print(paste0('Writing ', outras, '...'))
    writeRaster(kp, outras, datatype = "INT2S", overwrite=overwrite)
    return(outras)
  }) %>%
    do.call(rbind, .)

  return(kpl)
}
#------ mosaic_kriging -------------
mosaic_kriging <- function(in_kpathlist, in_filestructure, overwrite) {
  kpl <- in_kpathlist

  kplextents <- lapply(kpl, function(ras_path) {
    ras <- raster(ras_path)
    data.table(x=c(xmin(ras),xmax(ras)), y= c(ymin(ras),ymax(ras)))
  }) %>%
    rbindlist %>%
    .[, lapply(.SD, function(i) c(min(i), max(i)))]

  template <- raster(extent(kplextents[, x], kplextents[, y]))
  out_krigingtif = file.path(in_filestructure[['resdir']],
                             "prederror_krigingtest.tif")
  writeRaster(template, file=out_krigingtif, datatype='INT2S',
              format="GTiff", overwrite=overwrite)
  gdalUtils::mosaic_rasters(gdalfile=as.vector(kpl),
                            dst_dataset=out_krigingtif,
                            force_ot = "Int16",
                            of="GTiff")
  gdalinfo(out_krigingtif)

  return(out_krigingtif)
}

##### -------------------- Report functions -----------------------------------
#------ get_basemaps ------------
get_basemapswintri <- function(in_filestructure) {
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"

  wcountries <- rnaturalearth::ne_countries(
    scale = "medium", returnclass = "sf") %>%
    sfformat_wintri
  wland <- rnaturalearth::ne_download(
    scale = 110, type = 'land', category = 'physical') %>%
    sfformat_wintri
  wlakes <- rnaturalearth::ne_download(
    scale = 110, type = 'lakes', category = 'physical') %>%
    sfformat_wintri

  grat_wintri <-
    st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
    st_transform(crs = crs_wintri)

  # riv_query <- "SELECT HYRIV_ID, ORD_STRA, dis_m3_pyr FROM \"RiverATLAS_v10\" WHERE ORD_STRA > 4"
  # riv_lines <- st_read(dsn = dirname(in_filestructure['in_riveratlas']),
  #                      layer = basename(in_filestructure['in_riveratlas']),
  #                      query = riv_query)
  #
  # riv_wintri <- sfformat_wintri(riv_lines)
  # riv_simple <-  sf::st_simplify(riv_wintri, preserveTopology = TRUE, dTolerance = 1000)

  return(list(wcountries = wcountries,
              wland = wland,
              wlakes = wlakes,
              grat = grat_wintri))
  #,riv_simple = riv_simple
}

#------ ggrivers -----------------------
ggrivers <- function(in_basemaps) {
  p <- ggplot() +
    geom_sf(data = in_basemaps[['grat']],
            color = alpha("black", 1/5), size = 0.25/.pt) +
    geom_sf(data = in_basemaps[['wland']],
            color='#dee0e0') +
    geom_sf(data = in_basemaps[['wlakes']],
            alpha=1/3, fill='black', color=NA) +
    # geom_sf(data=in_basemaps['riv_simple'][riv_simple$ORD_STRA < 6,][1:10000,],
    #         aes(size=ORD_STRA), color=alpha('black', 1/20), show.legend = F) +
    # geom_sf(data=in_basemaps['riv_simple'][riv_simple$ORD_STRA >= 6,][1:10000,],
    #         aes(size=ORD_STRA), color=alpha('black', 1/10), show.legend = F) +
    scale_size_continuous(range=c(0.2, 1.0)) +
    coord_sf(datum = NA, expand=F,
             xlim=c(-14500000, 18000000), ylim=c(-7000000,8700000)) +
    theme_map() +
    theme(legend.position=c(0.1,0),
          legend.direction="horizontal",
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.key.height = unit(0.1,"in"))

  return(p)
}

#------ gggauges --------------------------
binarg <- c(30, 50, 100)
binvar <- 'totalYears_kept'

gggauges <- function(in_gaugepred, in_basemaps,
                     binarg) {
  gaugepred <- in_gaugepred %>%
    sfformat_wintri

  #Bin data
  dt <- as.data.table(gaugepred)
  byd <- dt[, binvar, with=F]
  l_freq    <- append(min(byd), binarg)
  u_freq    <- c(binarg, (max(byd, na.rm = TRUE) + 1))
  bins      <- length(binarg) + 1

  for (i in seq_len(bins)) {
    gaugepred[dt[,get(binvar)] >= l_freq[[i]] &
                dt[,get(binvar)] < u_freq[[i]], 'bin'] <- i

    if (i == bins) {
      gaugepred[dt[,get(binvar)]  == u_freq[i], 'bin'] <- i
    }
  }

  gaugepred$bin <- factor(gaugepred$bin, levels=seq_len(bins))


  #Create color scales
  cs_per <- scales::seq_gradient_pal('#45c6f7', '#152540', "Lab")(
    seq(0,1,length.out=bins)) #'#0F9FD6'
  cs_ir <- scales::seq_gradient_pal('#ffb88f','#940404', "Lab")( #'#ff9b52'  '#a32d18'
    seq(0,1,length.out=bins)) #


  #Subset data into perennial and intermittent rivers
  perennial_gauges <- gaugepred[gaugepred$intermittent=='0',]
  ir_gauges <-  gaugepred[gaugepred$intermittent=='1',]


  #Make histograms
  gggaugehist <- function(in_gdf, cs) {
    x_tick <- c(0, unique(as.numeric(levels(in_gdf$bin)))) + 0.5
    binarg_tick <- c(0, seq_len(bins))
    len <- length(x_tick)

    permean <- as.data.table(in_gdf)[,  mean(get(binvar))]
    permeanbindiff <- c(permean - l_freq)
    meanbin <- which(permeanbindiff > 0)[which.min(permeanbindiff[permeanbindiff  > 0])]
    meanpos <- meanbin + (permean-l_freq[meanbin])/(u_freq[meanbin]-l_freq[meanbin])

    gaugehist <- ggplot(in_gdf, aes(x=as.numeric(bin), fill=bin)) +
      geom_histogram(stat="count", width=1) +
      geom_vline(xintercept = meanpos - 1) +
      annotate(geom='text', angle=90,
               x=meanpos-0.6,
               y=0.5*max(as.data.table(in_gdf)[, .N, by=bin]$N),
               label=paste(round(permean), 'y')) +
      scale_fill_manual(values=cs) +
      scale_x_continuous(name = 'Years of data',
                         breaks = c(binarg_tick, x_tick),
                         labels = c(rep(c(""), len), c(l_freq, u_freq[bins]))) +
      # scale_y_continuous(name = paste0('Number of gauging stations (total: ',
      #                                  nrow(in_gdf),
      #                                  ')'))+
      coord_cartesian(clip='off', expand=c(0,0)) +
      theme_classic() +
      theme(legend.position = 'none',
            text = element_text(size=12),
            plot.background = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_text(vjust=3),
            axis.ticks.x = element_line(color = c(rep(NA, len - 1),
                                                  rep("black", len))))
    return(gaugehist)
  }

  histper <- gggaugehist(perennial_gauges, cs=cs_per)
  histir <- gggaugehist(ir_gauges, cs=cs_ir)


  #Plot it
  p_pr <- ggrivers(in_basemaps) +
    geom_sf(data=perennial_gauges, aes(color=bin), size=1.2, alpha=0.8) +
    scale_color_manual(values=cs_per) +
    coord_sf(datum = NA, expand=F,
             xlim=c(-16000000, 18000000), ylim=c(-7000000,8700000)) +
    theme(legend.position = 'none') +
    annotation_custom(grob = ggplotGrob(histper),
                      xmin = -16000000,
                      xmax = -8000000,
                      ymin = -9000000,
                      ymax = -1000000)


  p_ir <- ggrivers(in_basemaps) +
    geom_sf(data=ir_gauges, aes(color=bin), size=1) +
    scale_color_manual(values=cs_ir) +
    coord_sf(datum = NA, expand=F,
             xlim=c(-16000000, 18000000), ylim=c(-7000000,8700000)) +
    theme(legend.position = 'none') +
    annotation_custom(grob = ggplotGrob(histir),
                      xmin = -16000000,
                      xmax = -8000000,
                      ymin = -9000000,
                      ymax = -1000000)

  p_pr/p_ir
  return(p_pr/p_ir)
}

#------ formatscales ------------
formatscales <- function(in_df, varstoplot) {
  scales_x <- list(
    ari_ix_uav = scale_x_sqrt(expand=c(0,0)),
    cly_pc_uav = scale_x_continuous(labels=percent_format(scale=1), expand=c(0,0)),
    clz_cl_cmj = scale_x_continuous(limits=c(1,18), expand=c(0,0),
                                    breaks=seq(0,18)),
    cmi_ix_uyr = scale_x_continuous(),
    dis_m3_pyr = scale_x_sqrt(breaks=c(0, 10^2,
                                       10^(0:log10(max(in_df$dis_m3_pmn)))),
                              labels=c(0, 10^2,
                                       10^(0:log10(max(in_df$dis_m3_pmn)))),
                              expand=c(0,0)),
    dor_pc_pva = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    for_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    gla_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    kar_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    lka_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    pet_mm_uyr = scale_x_continuous(expand=c(0,0)),
    snw_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    run_mm_cyr = scale_x_continuous(expand=c(0,0)),
    swc_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    tmp_dc_uyr = scale_x_continuous(expand=c(0,0)),
    hdi_ix_cav = scale_x_continuous(expand=c(0,0)),
    hft_ix_c93 = scale_x_continuous(expand=c(0,0)),
    ORD_STRA = scale_x_continuous(expand=c(0,0)),
    gwt_m_cav = scale_x_continuous(expand=c(0,0)),
    ire_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0))
  ) %>%
    .[(names(.) %in% names(in_df)) & names(.) %in% varstoplot]
  #Only keep those variables that are actually in df and that we want to plot

  scales_y <- unlist(rep(list(scale_y_log10(expand=c(0,0))), labels = scientific_format(),
                         length(scales_x)),
                     recursive=F) %>%
    setNames(names(scales_x))
  scales_y[['dis_m3_pmn']] <- scale_y_sqrt(expand=c(0,0))

  coordcart <- lapply(varstoplot, function(var) {
    coord_cartesian(xlim=as.data.table(in_df)[, c(min(get(var), na.rm=T),
                                                  max(get(var), na.rm=T))])
  }) %>%
    setNames(names(scales_x))

  coordcart[['clz_cl_cmj']] <-  coord_cartesian(
    xlim=c(1,max(in_df$clz_cl_cmj)))
  coordcart[['kar_pc_use']] <-  coord_cartesian(
    xlim=c(0, 100))
  coordcart[['pet_mm_uyr']] <-  coord_cartesian(
    xlim=c(0, max(in_df$pet_mm_uyr)))
  coordcart[['ORD_STRA']] <-  coord_cartesian(
    xlim=c(1, 10))

  return(list(scales_x=scales_x, scales_y=scales_y, coordcart=coordcart))
}

#------ ggenvhist -------------
ggenvhist <- function(vartoplot, in_gaugedt, in_rivdt, in_predvars,
                      scalesenvhist, intermittent=TRUE) {
  print(vartoplot)
  if (intermittent) {
    vartoplot2 <- c(vartoplot, 'intermittent')
  }

  varname <- in_predvars[varcode==vartoplot, paste0(Attribute, ' ',
                                                    Keyscale,
                                                    Keystat,
                                                    ' (',unit,')')]

  penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
    geom_histogram(data=in_rivdt, bins=20, fill='lightgray') +
    geom_histogram(bins=20, fill='darkgray') +
    #aes(fill=intermittent), position = 'stack') - doesn't work with log or sqrt y scale
    # scale_fill_manual(values=c('#0F9FD6','#ff9b52'),
    #                   labels = c('Perennial', 'Intermittent'), name=NULL) +
    scalesenvhist$scales_x[[vartoplot]] +
    scalesenvhist$scales_y[[vartoplot]] +
    scalesenvhist$coordcart[[vartoplot]] +
    xlab(varname) +
    ylab('Count') +
    theme_classic() +
    theme(strip.background=element_rect(colour="white", fill='lightgray'),
          axis.title.y = element_blank(),
          axis.title = element_text(size=12))


  # if (which(vartoplot %in% varstoplot_hist)!=length(varstoplot_hist)) {
  #   penvhist <- penvhist +
  #     theme(legend.position='none')
  # }

  return(ggplotGrob(penvhist))
}

#------ layout_ggenvhist --------------------------
layout_ggenvhist <- function(in_rivernetwork, in_gaugepred, in_predvars) {
  varstoplot_hist <- c("ari_ix_uav", "cly_pc_uav", "clz_cl_cmj", "cmi_ix_uyr",
                       "dis_m3_pyr", "dor_pc_pva","for_pc_use", "gla_pc_use",
                       "kar_pc_use", "lka_pc_use", "pet_mm_uyr", "snw_pc_uyr",
                       "run_mm_cyr", "swc_pc_uyr", "tmp_dc_uyr", "hdi_ix_cav",
                       "hft_ix_c93", "UPLAND_SKM", "gwt_m_cav",  "ire_pc_use")

  scalesenvhist <- formatscales(in_df=in_rivernetwork, varstoplot=varstoplot_hist)

  # gaugeformat <- as.data.table(in_gaugepred)[, `:=`(
  #   ari_ix_cav = ari_ix_cav/100,
  #   ari_ix_uav = ari_ix_uav/100,
  #   cmi_ix_cmn = cmi_ix_cmn/100,
  #   cmi_ix_uyr = cmi_ix_uyr/100,
  #   dor_pc_pva = dor_pc_pva/100,
  #   lka_pc_cse = lka_pc_cse/10,
  #   lka_pc_use = lka_pc_use/10,
  #   tmp_dc_cmn =  tmp_dc_cmn/10,
  #   tmp_dc_cmx =  tmp_dc_cmx/10,
  #   tmp_dc_cyr =  tmp_dc_cyr/10,
  #   tmp_dc_uyr =  tmp_dc_uyr/10,
  #   gwt_m_cav = gwt_cm_cav/100
  # )]

  penvhist_grobs <- lapply(varstoplot_hist, ggenvhist,
                           in_gaugedt = gaugeformat,
                           in_rivdt = in_rivernetwork,
                           in_predvars = in_predvars,
                           scalesenvhist = scalesenvhist)
  do.call("grid.arrange", list(grobs=penvhist_grobs))
}

#------ tabulate_benchmarks ------------
tabulate_benchmarks <- function(in_bm, in_bmid) {

  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  print('Getting table content...')
  tbbm <- bm_msrtab(in_bm) %>%
    format_modelcompdat(typecomp=in_bmid)

  metrics_dat <- data.table(selection=character(),
                            type=character(),
                            learner_format=character(),
                            inner_folds=integer(),
                            inner_n_evals=integer(),
                            ntree = integer(),
                            numtrees = integer(),
                            alpha=numeric(),
                            mtry=integer(),
                            minnodesize=integer(),
                            fraction=numeric(),
                            samplefraction = numeric(),
                            `minor_weight|ratio`=numeric(),
                            minor_weight = numeric(),
                            ratio = numeric(),
                            npredictors=integer(),
                            outer_repeats=integer(),
                            outer_folds=integer()) %>%
    rbind(tbbm, use.names=TRUE, fill=TRUE)

  #Continue formatting table
  metrics_dat[is.na(ntree), ntree := numtrees]
  metrics_dat[is.na(fraction), fraction := samplefraction]
  metrics_dat[, `minor_weight|ratio` := fifelse(is.na(minor_weight),
                                                as.numeric(ratio),
                                                as.numeric(minor_weight))]


  #Create model specification tables
  setup_table <- metrics_dat[, .(selection, type, learner_format,
                                 inner_folds, inner_n_evals,
                                 alpha, mtry, minnodesize, fraction,
                                 `minor_weight|ratio`, npredictors,
                                 outer_repeats, outer_folds)] %>%
    .[!grepl('(CIF)|(MAXSTAT)', learner_format), alpha := NA] %>%
    .[grepl('MAXSTAT', learner_format), minnodesize := NA] %>%
    unique(by='learner_format')

  #Create model results table
  results_table <- metrics_dat[, .(selection, learner_format,
                                   resampling_id, outer_repeats, outer_folds,
                                   time_train=round(time_train),
                                   time_predict=round(time_predict),
                                   bacc=round(bacc, 3), threshold_class,
                                   bbrier=round(bbrier, 3), auc=round(auc, 3)
  )]

  return(list(setup=setup_table, results=results_table))
}

#------ formatmisclass_bm -------------
formatmisclass_bm <- function(in_bm, in_bmid) {
  #If path, read qs
  if (inherits(in_bm, "character")) {
    in_bm <- qs::qread(in_bm)
  }

  print('Getting table content...')
  thresh_dt <- threshold_dat(bmres =in_bm) %>%
    format_modelcompdat(in_bmid)
  return(thresh_dt)
}


#------ ggmisclass_bm -------------
ggmisclass_bm <- function(in_threshdts) {

  #
  threshplot_datall <- rbindlist(in_threshdts) %>%
    .[, list(value_mean = mean(value)),
      by=.(threshold_class, variable, learner_format, selection)]


  ggthresh_benchmark <- ggplot(threshplot_datall[variable %in% c('sens', 'spec'),],
                               aes(x=threshold_class, y=value_mean,
                                   linetype=variable,
                                   gourp=learner_format,
                                   color=learner_format)) +
    geom_line(size=1.2) +
    #scale_color_brewer(palette='Dark2') +
    scale_linetype(labels=c('Sensitivity (true positives)',
                            'Specificity (true negatives)')) +
    scale_x_continuous(expand=c(0,0), name='Threshold') +
    scale_y_continuous(expand=c(0,0), name='Value') +
    facet_wrap(~selection) +
    theme_bw() +
    theme(panel.spacing = unit(0.75, "cm"),
          legend.title = element_blank())

  return(ggthresh_benchmark)
}


#------ tabulate_globalsummary -----
#Example
# gad_tabledis <- tabulate_globalsummary(in_filestructure = filestructure,
#                                        idvars = 'gad_id_cmj',
#                                        castvar = 'dis_m3_pyr',
#                                        castvar_num = TRUE,
#                                        weightvar = 'LENGTH_KM',
#                                        valuevar = 'predbasic800cat',
#                                        valuevar_sub = 1,
#                                        binfunc = 'manual',
#                                        binarg = c(0.1, 1, 10, 100, 1000, 10000, 100000),
#                                        na.rm=T,
#                                        tidy = FALSE)

tabulate_globalsummary <- function(in_filestructure, idvars,
                                   castvar, castvar_num=TRUE,
                                   weightvar,
                                   valuevar, valuevarsub,
                                   binfunc=NULL, binarg=NULL, bintrans=NULL,
                                   na.rm=T, tidy=FALSE) {

  #Import global predictions
  rivpred <- fread(in_filestructure['out_riveratlas'])
  #Columns to import from full network
  incols <- c('HYRIV_ID', castvar, idvars, valuevar, weightvar)
  #Import global river network and join to predictions
  riveratlas <- fread_cols(file_name=in_filestructure['in_riveratlas'],
                           cols_tokeep = incols) %>%
    .[rivpred, on='HYRIV_ID']

  #If global administrative boundaries were selected for idvar, get country names
  if ('gad_id_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(file.path(in_filestructure['datdir'],
                                            'HydroATLAS',
                                            'HydroATLAS_v10_Legends.xlsx'),
                                  sheet='gad_id') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='gad_id_cmj', by.y='Country_ID') %>%
      .[, gad_id_cmj := Country_Name]
  } else if ('fmh_cl_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(file.path(in_filestructure['datdir'],
                                            'HydroATLAS',
                                            'HydroATLAS_v10_Legends.xlsx'),
                                  sheet='fmh_cl') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='fmh_cl_cmj', by.y='MHT_ID') %>%
      .[, fmh_cl_cmj := MHT_Name]
  } else if ('tbi_cl_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(file.path(in_filestructure['datdir'],
                                            'HydroATLAS',
                                            'HydroATLAS_v10_Legends.xlsx'),
                                  sheet='tbi_cl') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='tbi_cl_cmj', by.y='Biome_ID') %>%
      .[, tbi_cl_cmj := Biome_Name]
  } else if ('clz_cl_cmj' %in% incols) {
    gadnames <- readxl::read_xlsx(file.path(in_filestructure['datdir'],
                                            'HydroATLAS',
                                            'HydroATLAS_v10_Legends.xlsx'),
                                  sheet='clz_cl') %>%
      setDT

    riveratlas <- merge(riveratlas, gadnames, by.x='clz_cl_cmj', by.y='GEnZ_ID') %>%
      .[, clz_cl_cmj := GEnZ_Name]
  }

  #Bin castvar if needed
  if (!is.null(binfunc) & !is.null(binarg)) {
    riveratlas <- bin_dt(in_dt = riveratlas, binvar = castvar, valuevar = valuevar,
                         binfunc = binfunc, binarg = binarg, bintrans = bintrans,
                         na.rm = na.rm)
    castvar = 'bin_lformat'
    castvar_num = FALSE
  }


  #Compute overall number or weight for all combinations of castvar, valuevar and idvars
  #e.g. total number or river length for each flow state category for each country and river order
  statall <- riveratlas[if (na.rm) {!is.na(eval(valuevar))},
                        if (!is.null(weightvar)) sum(get(weightvar)) else .N,
                        by=c(eval(castvar), eval(valuevar), eval(idvars))] %>%
    rbind(
      .[, list('World', sum(V1)), by=c(eval(castvar), eval(valuevar))] %>%
        setnames(c('V1', 'V2'), c(idvars, 'V1'))
    )


  #Compute for each cast and id var, the percentage of the valuevar that is of the category of interest
  #(e.g. for each country and river order, percentage of river length that is intermittent)
  tidyperc <- statall[, 100*.SD[get(valuevar)==valuevarsub, sum(V1)]/sum(V1),
                      by=c(eval(castvar), eval(idvars))]

  #If cast variable is a numeric (e.g. Strahler Order), sort table in increasing order
  if (castvar_num) {
    tidyperc[, eval(castvar) := as.numeric(as.character(get(eval(castvar))))]
  }
  setorderv(tidyperc, cols=c(castvar, idvars))


  #Compute totals for each valuevar and idvar summed across castvar
  #e.g. (total length of rivers for each country, and of each flow regime)
  tidytotal <- statall[, list('Total intermittency (%)',
                              100*.SD[get(valuevar)==valuevarsub, sum(V1)]/sum(V1)),
                       by=eval(idvars)] %>%
    rbind(  statall[, list('Total stream length (10^3 km)', sum(V1)/1000),
                    by=c(eval(idvars))]) %>%
    setnames(c('V1', 'V2'), c(castvar, 'V1'))
  totalcols <- unique(tidytotal[, get(eval(castvar))]) #Name of total cols

  #Prepare casting formula for formatting
  castformula <- as.formula(paste0(paste(idvars, collapse='+'),
                                   '~',
                                   paste(castvar, collapse='+')))

  #Merge percentages and totals, then cast table if tidy != TRUE
  tidypercformat <- tidyperc[, eval(castvar) := as.factor(get(eval(castvar)))] %>%
    rbind(tidytotal)

  if (!tidy) {
    tidyperc_cast <- dcast(tidypercformat, castformula, value.var = 'V1') %>%
      .[, (totalcols) := lapply(.SD, function(x) fifelse(is.na(x), 0 , x)), #Replace NAs in total cols by 0
        .SDcols = totalcols]

    tidyperc_format <- rbind(tidyperc_cast[get(eval(idvars)) != 'World',],
                             tidyperc_cast[get(eval(idvars)) == 'World',])
  }

  return(tidyperc_format)
}
#------ compare_fr --------------------------------------
compare_fr <- function(in_filestructure, in_rivernetwork, binarg) {
  rivpred <- fread(in_filestructure['out_riveratlas']) %>%
    .[in_rivernetwork[, c('HYRIV_ID', 'HYBAS_L12', 'LENGTH_KM', 'dis_m3_pyr', 'UPLAND_SKM'),
                      with=F], on='HYRIV_ID']
  in_frdir <- file.path(in_filestructure[['resdir']],
                        'Comparison_databases/france.gdb')
  in_netpath <- file.path(in_frdir, 'network')
  in_baspath <- file.path(in_frdir, 'hydrobasins12')
  valuevarsub <- "1"

  net <- st_read(dsn = dirname(in_netpath),
                 layer = basename(in_netpath))

  bas <- st_read(dsn = dirname(in_baspath),
                 layer = basename(in_baspath)) %>%
    .[, 'HYBAS_ID', with=F]

  rivpredsub <- merge(rivpred, bas, by.x="HYBAS_L12", by.y="HYBAS_ID", all.x=F) %>%
    .[, UPLAND_SKM := round(UPLAND_SKM)]

  binlabels <- label_manualbins(binarg=binarg,
                                minval=min(net$rhtvs2_all_phi_qclass_SURF_BV))

  tidyperc_fr <- formathistab(in_dt = net,
                              castvar = "rhtvs2_all_phi_qclass_SURF_BV",
                              valuevar = "INT_RF_txt_V1",
                              valuevarsub = valuevarsub,
                              weightvar = "rhtvs2_all_phi_qclass_LONG_",
                              binfunc = 'manual',
                              binarg =  binarg,
                              binlabels = binlabels,
                              datname = 'Snelder et al. (2013) predictions') %>%
    .[, binsumlength := binsumlength/1000]

  tidyperc_riv  <- formathistab(in_dt = rivpredsub,
                                castvar = 'UPLAND_SKM',
                                valuevar = 'predbasic800cat',
                                valuevarsub = valuevarsub,
                                weightvar = 'LENGTH_KM',
                                binfunc = 'manual',
                                binarg =  binarg,
                                binlabels = binlabels,
                                datname = 'Global predictions')

  datmerge <- rbind(tidyperc_fr, tidyperc_riv) %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]

  return(
    ggcompare(datmerge, binarg)
  )
}

#------ compare_us ----------------
compare_us <- function(in_filestructure, in_rivernetwork, binarg) {
  rivpred <- fread(in_filestructure['out_riveratlas']) %>%
    .[in_rivernetwork[, c('HYRIV_ID', 'HYBAS_L12', 'LENGTH_KM', 'dis_m3_pyr',
                          'UPLAND_SKM'),
                      with=F], on='HYRIV_ID']
  in_dir <- file.path(in_filestructure[['datdir']],
                      'Comparison_databases/US')
  in_netpath_hr <- file.path(in_filestructure[['datdir']],
                             'Comparison_databases/US',
                             'NHDhr_attris.csv')
  in_netpath_mr <- file.path(in_filestructure[['datdir']],
                             'Comparison_databases/US',
                             'NHDmr_attris.csv')
  in_baspath <- file.path(in_filestructure[['resdir']],
                          'Comparison_databases/us.gdb',
                          'hydrobasins12')
  valuevarsub <- "1"

  #55800: Artificial path
  #46000: Stream/River
  #46003: Stream/River: Hydrographic Category = Intermittent
  #46006: Stream/River: Hydrographic Category = Perennial
  #46007: Stream/River: Hydrographic Category = Ephemeral

  nethr <- fread(in_netpath_hr,
               colClasses = c('character', 'character', 'numeric', 'character',
                              'integer', 'numeric', 'numeric'))
  nethr[, HUC8 := substr(as.character(ReachCode), 1, 8)]

  #Compare total

  #Check the percentage of each type of line in the NHD by HUC8 and Drainage area
  # nethr[FCode %in% c(55800, 46000, 46003, 46006, 46007),
  #     HUC8sum := .N, by=HUC8]
  # FCode_HUC8 <- nethr[FCode %in% c(55800, 46000, 46003, 46006, 46007),
  #                   as.integer(100*max(.SD[, .N]/HUC8sum)), by=.(HUC8, FCode)]
  #
  nethrbinned <- bin_dt(in_dt=nethr[TotDASqKm > 0,],
                      binvar='TotDASqKm',
                      valuevar='intermittent',
                      binfunc='manual',
                      binarg=binarg)
  # nethrbinned[FCode %in% c(55800, 46000, 46003, 46006, 46007),
  #           DAsum := .N, by=bin_lformat]
  # FCode_DA <- nethrbinned[FCode %in% c(55800, 46000, 46003, 46006, 46007),
  #                       as.integer(100*max(.SD[, .N]/DAsum)), by=.(bin_lformat, FCode)]

  #Exclude all HUC8s with more than 5% of non-identified FCode (46000)
  nethrsub <- nethrbinned[!(HUC8 %in% nethr[FCode==46000 & V1>10,unique(HUC8)]) &
                        TotDASqKm > -9999,]

  #Create intermittency variable
  nethrsub[, intermittent := fifelse(FCode %in% c(46003, 46007), '1', '0')]
  nethrsub[FCode == 55800, intermittent := -1]
  nethrsub[!(FCode %in% c(55800, 46003, 46007, 46006)), intermittent := NA]

  #Get HydroSHEDS basins that overlap with selected NHD HUC8s
  bas <- st_read(dsn = dirname(in_baspath),
                 layer = basename(in_baspath)) %>%
    .[, c('HYBAS_ID', 'HUC8'), with=F]

  #Join HydroSHEDS basins with RiverATLAS network and subselect network to match NHD selection
  rivpredsub <- merge(rivpred, bas, by.x="HYBAS_L12", by.y="HYBAS_ID", all.x=F) %>%
    .[, UPLAND_SKM := round(UPLAND_SKM)] %>%
    .[HUC8 %in% unique(nethrsub$HUC8),]

  #Get bin labels
  binlabels <- label_manualbins(binarg=binarg,
                                minval=min(nethrsub$TotDASqKm))

  #Compute bin statistics ncluding "artificial flow path"
  tidyperc_usall <- formathistab(in_dt = nethrsub[!is.na(intermittent),],
                                 castvar = "TotDASqKm",
                                 valuevar = "intermittent",
                                 valuevarsub = valuevarsub,
                                 weightvar = "LengthKM",
                                 binfunc = 'manual',
                                 binarg =  binarg,
                                 binlabels = binlabels,
                                 datname = 'U.S. National Hydrography Dataset') %>%
    .[, binsumlength := binsumlength]

  #Compute bin statistics excluding "artificial flow path"
  tidyperc_usnoartificial <- formathistab(in_dt = nethrsub[!is.na(intermittent) &
                                                           intermittent != -1,],
                                          castvar = "TotDASqKm",
                                          valuevar = "intermittent",
                                          valuevarsub = valuevarsub,
                                          weightvar = "LengthKM",
                                          binfunc = 'manual',
                                          binarg =  binarg,
                                          binlabels = binlabels,
                                          datname = 'U.S. National Hydrography Dataset') %>%
    .[, binsumlength := binsumlength]

  #Compute bin statistics for river networks
  tidyperc_riv  <- formathistab(in_dt = rivpredsub,
                                castvar = 'UPLAND_SKM',
                                valuevar = 'predbasic800cat',
                                valuevarsub = valuevarsub,
                                weightvar = 'LENGTH_KM',
                                binfunc = 'manual',
                                binarg =  binarg,
                                binlabels = binlabels,
                                datname = 'Global predictions')

  #Merge statistics
  datmerge_all <- rbind(tidyperc_usall, tidyperc_riv) %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]

  datmerge_noartificial <- rbind(tidyperc_usnoartificial, tidyperc_riv) %>%
    setorder(bin) %>%
    .[, binformat := factor(binformat, levels=unique(binformat))]

  #Create plot
  return(
    ggcompare(datmerge_all) +
      geom_bar(data=datmerge_noartificial, alpha=0.75,
               stat='identity', position='dodge')
  )
}


#------ qc_pnw ------------
qc_pnw <- function(in_filestructure, in_rivernetwork) {
  rivpred <- fread(in_filestructure['out_riveratlas']) %>%
    .[in_rivernetwork[, c('HYRIV_ID', 'HYBAS_L12', 'LENGTH_KM', 'dis_m3_pyr',
                          'UPLAND_SKM'),
                      with=F], on='HYRIV_ID']
  in_dir <- file.path(in_filestructure[['resdir']],
                      'Insitu_databases/pnw.gdb')
  in_refpts <- file.path(in_dir,
                         'StreamflowPermObs_final')
  in_fulldat <- file.path(in_dir,
                          'StreamflowPermObs_sub')

  valuevarsub <- "1"

  #Georeferenced/Snapped points to RiverATLAS network after removing duplicate observations at single sites
  refpts <- st_read(dsn = dirname(in_refpts),
                    layer = basename(in_refpts)) %>%
    setDT %>%
    setkey('OBJECTID_3')

  #All observations (excluding those that were not kept by PROSPER (aside from more recent osb, see python code for details)
  #with duplicate records for single locations (multiple observations for different dates)
  fulldat <- st_read(dsn = dirname(in_fulldat),
                     layer = basename(in_fulldat)) %>%
    setDT %>%
    setkey('OBJECTID_1')

  #Join points and compute # of observations, min and max months of observations
  refpts_full <- merge(refpts, fulldat[, c('OBJECTID', 'dupligroup'), with=F],
                       by='OBJECTID', all.y=F) %>%
    .[fulldat, on='dupligroup', allow.cartesian=T] %>%
    .[!is.na(distatlas),]

  refpts_stats <- refpts_full[, `:=`(nobs=.N,
                                     minmonth = min(Month),
                                     maxmonth = max(Month)),
                              by = dupligroup] %>%
    .[!duplicated(dupligroup),]

  #Merge points with rivernetwork by HYRIV_ID
  refpts_join <- merge(refpts_stats,
                       rivpred[, .(HYRIV_ID, HYBAS_L12, predbasic800,
                                   predbasic800cat)],
                       by='HYRIV_ID', all.y=F) %>%
    .[, refinter := fifelse(Category == 'Non-perennial', 1, 0)] %>% #Assign ephemeral and intermittent categories to 1, perennial to 0 for PNW obs
    .[, prederror := predbasic800 - refinter] #Compute prediction error

  #Scatterpoint of x = drainage area, y = predprobability - refinter
  refpts_joinmelt <- melt(
    refpts_join[, .(OBJECTID, prederror, UPLAND_SKM, dis_m3_pyr, nobs, refinter)],
    id.vars = c('OBJECTID', 'prederror', 'refinter'))


  levels(refpts_joinmelt$variable) <- c('Drainage area (km2)',
                                        'Discharge (m3)',
                                        '# of field obs.')


  colorpal <- c('#1f78b4', '#ff7f00')
  rectdf <- data.table(
    xmin=rep(0.01, 4),
    xmax=rep(Inf, 4),
    ymin=c(-1, -0.5, 0, 0.5),
    ymax=c(-0.5, 0, 0.5, 1),
    fillpal = rep(colorpal, 2)
  )

  pnw_qcplot <- ggplot(refpts_joinmelt) +
    geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax,
                               ymin=ymin, ymax=ymax, fill=fillpal),
              alpha=1/4) +
    scale_fill_manual(values=colorpal,
                      name='Predicted regime',
                      labels = c('Perennial', 'Intermittent')) +
    geom_point(aes(x=value, y=prederror, color=factor(refinter)),
               alpha=1/5) +
    geom_hline(yintercept=0, alpha=1/2) +
    new_scale_fill() +
    geom_smooth(aes(x=value, y=prederror, color=factor(refinter)),
                method='gam', formula = y ~ s(x, k=4)) +
    scale_x_sqrt() +
    geom_hline(yintercept=0) +
    scale_color_manual(values=c('#1f78b4', '#ff7f00'),
                       name='Observed regime',
                       labels = c('Perennial', 'Intermittent')) +
    coord_cartesian(expand=FALSE, clip='off') +
    labs(x='Value', y='Intermittency Prediction Residuals (IPR)') +
    theme_classic() +
    theme(legend.position = c(0.9, 0.45),
          legend.background = element_blank()) +
    facet_wrap(~variable, scales ='free_x')

  return(pnw_qcplot)
}


########## END #############
