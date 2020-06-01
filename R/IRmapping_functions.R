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
  #check <- riveratlas[cmi_ix_uyr == -9999,] #All have precipitation = 0
  in_dt2[cmi_ix_uyr == -9999, cmi_ix_uyr := 0]

  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  in_dt2[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  in_dt2[snw_pc_cmx == -9999, snw_pc_cmx:=0]

  #check <- in_dt2[is.na(sgr_dk_rav),]

  print('Number of NA values per column')
  colNAs<- in_dt2[, lapply(.SD, function(x) sum(is.na(x) | x==-9999))]
  print(colNAs)

  #Convert -9999 values to NA
  for (j in which(sapply(in_dt2,is.numeric))) { #Iterate through numeric column indices
    set(in_dt2,which(in_dt2[[j]]==-9999),j, NA)} #Set those to 0 if -9999

  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  pre_mcols <- paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)) #Monthly precipitation columns
  in_dt2[, `:=`(pre_mm_cmn = do.call(pmin, c(.SD, list(na.rm=TRUE))), #Compute minimum and maximum catchment precipitation
                pre_mm_cmx = do.call(pmax, c(.SD, list(na.rm=TRUE)))),
         .SDcols= pre_mcols] %>%
    #       min/max monthly catchment precip (while dealing with times when )
    .[, `:=`(pre_mm_cvar= fifelse(pre_mm_cmx==0, 0, pre_mm_cmn/pre_mm_cmx),
             #min/max monthly watershed discharge
             dis_mm_pvar=fifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx),
             #min monthly/average yearly watershed discharge
             dis_mm_pvaryr=fifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),
             #catchment average elv - watershec average elev
             ele_pc_rel = fifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav))]

  in_dt2[, cmi_ix_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
         .SDcols= paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly cmi
    .[, swc_pc_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
      .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] #Get minimum monthly swc


  #Scale variables based on HydroATLAS documentation
  in_dt2[, `:=`(
    ari_ix_cav = ari_ix_cav/100,
    ari_ix_uav = ari_ix_uav/100,
    cmi_ix_cmn = cmi_ix_cmn/100,
    cmi_ix_uyr = cmi_ix_uyr/100,
    dor_pc_pva = dor_pc_pva/100,
    lka_pc_cse = lka_pc_cse/10,
    lka_pc_use = lka_pc_use/10,
    tmp_dc_cmn =  tmp_dc_cmn/10,
    tmp_dc_cmx =  tmp_dc_cmx/10,
    tmp_dc_cyr =  tmp_dc_cyr/10,
    tmp_dc_uyr =  tmp_dc_uyr/10,
    gwt_m_cav = gwt_cm_cav/100
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
#'  \link{ggpd}
#'
#' @section Warning:
#' Has only been tested on \link[mlr3learners]{mlr_learners_classif.ranger}
#'
#' @section documentation to-do:
#' Can add an example down the line, add source.
#'
#' @export

extract_pd_nestedrf <- function(learner_id=1, in_rftuned, datdf, selcols, ngrid) {
  in_mod <- in_rftuned$learners[[learner_id]]
  #in_mod <- nestedresamp_ranger$learners[[1]]

  if (inherits(in_mod$learner, "GraphLearner")) {
    in_fit <- in_mod$learner$model$classif.ranger$model
  } else {
    in_fit <- in_mod$learner$model
  }

  foldperf <- extract_impperf_nestedrf(in_mod, imp=F, perf=T, pvalue=F)

  # selcols <- in_vimp_plot$data %>% #Can use that if extracting from tunredrf is expensive
  #   setorder(-imp_wmean) %>%
  #   .[colnums, variable]

  vargrid <- combn(selcols, 2, simplify=F) %>%
    do.call(rbind, .)

  ngridvec <- c(ngrid, ngrid)

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

##### -------------------- Workflow functions ---------------------------------
#------ def_filestructure -----------------
#' Parallel wrapper for comp_durfreq
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
  in_riveratlas_meta <- file.path(datdir, 'HydroATLAS', 'HydroATLAS_metadata_MLM.xlsx')
  # Output geopackage of hydrometric stations with appended predicted intermittency class
  out_gauge <- file.path(resdir, 'GRDCstations_predbasic800.gpkg')
  # River atlas attribute data1
  in_riveratlas <- file.path(resdir, 'RiverATLAS_v10tab.csv')
  # River atlas predictions table
  out_riveratlas <- file.path(resdir, 'RiverATLAS_predbasic800.csv')

  return(c(rootdir=rootdir, datdir=datdir, resdir=resdir, outgdb=outgdb,
           in_gaugep=in_gaugep, in_gaugedir=in_gaugedir, in_riveratlas_meta=in_riveratlas_meta,
           out_gauge=out_gauge, in_riveratlas=in_riveratlas,
           out_riveratlas = out_riveratlas))
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

read_gaugep <- function(in_filestructure, dist) {
  #Import gauge stations and only keep those < dist m from a HydroSHEDS reach
  return(st_read(dsn=dirname(in_filestructure['in_gaugep']),
                 layer=basename(in_filestructure['in_gaugep'])) %>%
           .[.$station_river_distance<dist,]
  )
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
    .[, c('X', 'Y') := as.data.table(sf::st_coordinates(Shape))] %>%
    comp_derivedvar #Compute derived variables, rescale some variables, remove -9999

  return(gaugestats_join)
}
#------ selectformat_predvars -----------------
selectformat_predvars <- function(in_filestructure, in_gaugestats) {
  #---- List predictor variables ----
  predcols<- c(
    'ORD_STRA',
    'dis_m3_pyr',
    'dis_m3_pmn',
    'dis_m3_pmx',
    'dis_mm_pvar',
    'dis_mm_pvaryr',
    'run_mm_cyr',
    'inu_pc_umn',
    'inu_pc_umx',
    'inu_pc_cmn',
    'lka_pc_cse',
    'lka_pc_use',
    'dor_pc_pva',
    'gwt_m_cav',
    'ele_pc_rel',
    # 'sgr_dk_rav', #Don't use stream gradient as data are missing for all of Greenland due to shitty DEM
    'clz_cl_cmj',
    'tmp_dc_cyr',
    'tmp_dc_cmn',
    'tmp_dc_cmx',
    'tmp_dc_uyr',
    'pre_mm_uyr',
    'pre_mm_cvar',
    'pre_mm_cmn',
    'pet_mm_uyr',
    'ari_ix_uav',
    'ari_ix_cav',
    'cmi_ix_uyr',
    'cmi_ix_cmn',
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
    'cly_pc_uav',
    'cly_pc_cav',
    'slt_pc_uav',
    'slt_pc_cav',
    'snd_pc_uav',
    'snd_pc_cav',
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
    'gdp_ud_cav',
    'hdi_ix_cav')


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
                                    'Elevation catchment average - watershed average'),
                          varcode=c('pre_mm_cvar',
                                    'dis_mm_pvar', 'dis_mm_pvaryr',
                                    'ele_pc_rel'))

  predcols_dt <- merge(data.table(varcode=predcols),
                       rbind(meta_format, addedvars, fill=T),
                       by='varcode', all.y=F)

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

#------ resample_learner ------------------
#Directly run resample function from mlr3 on in_task, in_learner, in_resampling
dynamic_resample <- function(task, learner, resampling, type,
                             store_models = TRUE) {
  if (is.list(learner)) {
    learner <- learner[[1]]
  }

  if ((learner$task_type == 'classif' & type=='classif') |
      (learner$task_type == 'regr' & type=='regr')) {
    resmp_rs <- resample(task, learner, resampling, store_models)
    return(resmp_rs)
  }
}



#------ combined resample results into benchmark results -------------
combine_bm <- function(in_resampleresults) {
  #When tried as_benchmark_result.ResampleResult, got "Error in setcolorder(data, slots) :
  # x has some duplicated column name(s): uhash. Please remove or rename the
  # duplicate(s) and try again.". SO use this instead
  if (length(in_resampleresults) > 1) {
    bmres_list <- lapply(in_resampleresults,
                         function(rsmpres) {
                           BenchmarkResult$new(rsmpres$data)})

    bmrbase = bmres_list[[1]]
    for (i in 2:length(bmres_list)) {
      if (in_resampleresults[[i]]$task$task_type ==
          in_resampleresults[[1]]$task$task_type) {
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

  return(bmrbase)
}

#------ benchmark_classif -----------------
benchmark_classif <- function(in_tasks,
                              insamp_nfolds, insamp_neval, insamp_nbatch,
                              outsamp_nrep, outsamp_nfolds) {

  #Create task vars
  task_classif <- in_tasks$classif

  #---------- Create learners --------------------------------------------------
  #Compute ratio of intermittent to perennial observations
  imbalance_ratio <- get_oversamp_ratio(task_classif)$ratio

  #Create basic learner
  lrn_ranger <- mlr3::lrn('classif.ranger',
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
  lrn_cforest <- mlr3::lrn('classif.cforest',
                           ntree = 500,
                           fraction = 0.632,
                           replace = FALSE,
                           alpha = 0.05,
                           mtry = round(sqrt(length(task_classif$feature_names))),
                           predict_type = "prob")

  #Create mlr3 pipe operator to oversample minority class based on major/minor ratio
  #https://mlr3gallery.mlr-org.com/mlr3-imbalanced/
  #https://mlr3pipelines.mlr-org.com/reference/mlr_pipeops_classbalancing.html
  #Sampling happens only during training phase.
  po_over <- mlr3pipelines::po("classbalancing", id = "oversample", adjust = "minor",
                               reference = "minor", shuffle = TRUE,
                               ratio = imbalance_ratio)
  #table(po_over$train(list(task_classif))$output$truth()) #Make sure that oversampling worked

  #Create mlr3 pipe operator to put a higher class weight on minority class
  po_classweights <- mlr3pipelines::po("classweights", minor_weight = imbalance_ratio)

  #Create graph learners so that oversampling happens systematically upstream of all training
  lrn_ranger_overp <- mlr3pipelines::GraphLearner$new(po_over %>>% lrn_ranger)
  lrn_cforest_overp <- mlr3pipelines::GraphLearner$new(po_over %>>% lrn_cforest)

  #Create graph learners so that class weighin happens systematically upstream of all training
  lrn_ranger_weight <- mlr3pipelines::GraphLearner$new(po_classweights %>>% lrn_ranger)
  lrn_cforest_weight <- mlr3pipelines::GraphLearner$new(po_classweights  %>>% lrn_cforest)

  #---------- Set up inner resampling ------------------------------------------
  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)

    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T),
                   lower = 5,
                   upper = floor(length(task_classif$feature_names)/2)), #Half number of features
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

  #Define performance measure
  measure_rf_class = msr("classif.bacc") #use balanced accuracy as objective function

  #Define termination rule
  evalsn = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds

  #Define hyperparameter tuner wrapper for inner sampling
  learns_classif = list(
    #Standard ranger rf without oversampling
    AutoTuner$new(learner= lrn_ranger,
                  resampling = rcv_rf,
                  measures = measure_rf_class,
                  tune_ps = regex_tuneset(lrn_ranger),
                  terminator = evalsn,
                  tuner =  tnr("random_search",
                               batch_size = insamp_nbatch)), #batch_size determines level of parallelism

    #Standard ranger rf with oversampling
    AutoTuner$new(learner= lrn_ranger_overp,
                  resampling = rcv_rf,
                  measures = measure_rf_class,
                  tune_ps = regex_tuneset(lrn_ranger_overp),
                  terminator = evalsn,
                  tuner =  tnr("random_search",
                               batch_size = insamp_nbatch)),

    #Standard ranger rf with class weights
    AutoTuner$new(learner= lrn_ranger_weight,
                  resampling = rcv_rf,
                  measures = measure_rf_class,
                  tune_ps = regex_tuneset(lrn_ranger_weight),
                  terminator = evalsn,
                  tuner =  tnr("random_search",
                               batch_size = insamp_nbatch)), #batch_size determines level of parallelism

    #CIF rf without oversampling
    lrn_cforest,

    #CIF rf with oversampling
    lrn_cforest_overp,

    #CIF rf with class
    lrn_cforest_weight
  )
  names(learns_classif) <-mlr3misc::map(learns_classif, "id")

  #---------- Set up outer resampling benchmarking -----------------------------

  #Perform outer resampling, keeping models for diagnostics later
  outer_resampling = rsmp("repeated_cv",
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)

  #Run outer resampling and benchmarking on classification learners
  nestedresamp_bmrdesign_classif <- benchmark_grid(
    tasks = task_classif,
    learners = learns_classif,
    resamplings = outer_resampling)

  nestedresamp_bmrout_classif <- benchmark(
    nestedresamp_bmrdesign_classif, store_models = TRUE)

  return(
    list(
      bm_classif = nestedresamp_bmrout_classif,
      bm_tasks = list(task_classif=task_classif),
      measure_classif = measure_rf_class
    )
  )
}
#------ benchmark_regr -----------------
benchmark_regr <- function(in_tasks,
                           insamp_nfolds, insamp_neval, insamp_nbatch,
                           outsamp_nrep, outsamp_nfolds) {

  task_regr <- in_tasks$regr
  task_regrover <- in_tasks$regover

  #---------- Create learners --------------------------------------------------
  #Create regression learner with maxstat. Represents an approximation of
  #classification learner with probabilities as the prediction type and a simili-
  #conditional-inference forest tweak to correct for the over-representation
  #of variables with many values over those with few categories
  lrn_ranger_maxstat <- mlr3::lrn('regr.ranger',
                                  num.trees=500,
                                  sample.fraction = 0.632,
                                  min.node.size = 10,
                                  replace=FALSE,
                                  splitrule = 'maxstat',
                                  importance = 'impurity_corrected',
                                  respect.unordered.factors = 'order')

  #---------- Set up inner resampling ------------------------------------------
  #Define paramet space to explore
  regex_tuneset <- function(in_lrn) {
    prmset <- names(in_lrn$param_set$tags)

    tune_rf <- ParamSet$new(list(
      ParamInt$new(grep(".*mtry", prmset, value=T),
                   lower = 5,
                   upper = floor(length(task_classif$feature_names)/2)), #Half number of features
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

  #Define performance measure
  measure_rf_reg = msr("regr.mae")

  #Define termination rule
  evalsn = term("evals", n_evals = insamp_neval) #termine tuning after 20 rounds

  #Define hyperparameter tuner wrapper for inner sampling
  learns_regr = list(
    #Regression rf with MAXSTAT
    AutoTuner$new(learner= lrn_ranger_maxstat,
                  resampling = rcv_rf,
                  measures = measure_rf_reg,
                  tune_ps = regex_tuneset(lrn_ranger_maxstat),
                  terminator = evalsn,
                  tuner =  tnr("random_search",
                               batch_size = insamp_nbatch))
  )

  #---------- Set up outer resampling benchmarking -----------------------------

  #Perform outer resampling, keeping models for diagnostics later
  outer_resampling = rsmp("repeated_cv",
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)
  #Run outer resampling and benchmarking on regression learners
  nestedresamp_bmrdesign_regr <- benchmark_grid(
    tasks = list(task_regr, task_regrover),
    learners = learns_regr,
    resamplings = outer_resampling)

  nestedresamp_bmrout_regr <- benchmark(
    nestedresamp_bmrdesign_regr, store_models = TRUE)

  return(
    list(
      bm_regr = nestedresamp_bmrout_regr,
      bm_tasks = list(task_regr=task_regr,
                      task_regrover=task_regrover),
      measure_regr = measure_rf_reg
    )
  )
}

#------ analyze_benchmark -----------------
analyze_benchmark <- function(in_bm, in_measure) {

  print(in_bm$aggregate(in_measure))
  boxcomp <- mlr3viz::autoplot(in_bm, measure = in_measure)

  print(paste('It took',
              in_bm$aggregate(msr('time_both'))$time_both,
              'seconds to train and predict with the',
              in_bm$aggregate(msr('time_both'))$learner_id,
              'model...'))

  bmdt <- as.data.table(in_bm)

  if (in_bm$task_type == 'regr') {
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

    gout <- ggmisclass(in_predictions = subpred) +
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
#------ benchmark_featsel -----------------
benchmark_featsel <- function(in_rf, in_task, in_measure,
                              pcutoff = 0.1,
                              insamp_nfolds =  NULL, insamp_nevals = NULL,
                              outsamp_nrep = NULL, outsamp_nfolds =  NULL,
                              outsamp_nfolds_sp =  NULL) {
  #pcutoff is the p_value above which features are removes

  #Apply feature/variable selection
  vimp <- weighted_vimportance_nestedrf(
    rfresamp = in_rf$resample_result(uhash=unique(as.data.table(in_rf)$uhash)),
    pvalue = TRUE) %>%
    .[,imp_wmeanper := imp_wmean/sum(imp_wmean)]

  task_featsel <- in_task$clone()$select(
    vimp[imp_pvalue <= pcutoff, as.character(varnames)])
  task_featsel$id <- paste0(in_task$id, '_featsel')

  if (is.null(outsamp_nfolds_sp)) {
    outsamp_nfolds_sp = outsamp_nfolds
  }

  #Set up outer resampling including
  outer_resampling = rsmp("repeated_cv",
                          repeats = outsamp_nrep,
                          folds = outsamp_nfolds)

  outer_resamplingsp = rsmp("repeated-spcv-coords",
                            repeats = outsamp_nrep,
                            folds = outsamp_nfolds_sp) #Create 20 folds (visualize)

  #Run outer resampling and benchmarking on classification learners
  bmrdesign_featsel <- benchmark_grid(
    tasks = list(in_task, task_featsel),
    learners = in_rf$learners$learner[[1]],
    resamplings = list(outer_resampling, outer_resamplingsp))

  bmrout_featsel <- benchmark(
    bmrdesign_featsel, store_models = TRUE)

  bm_analysis <- analyze_benchmark(in_bm = bmrout_featsel,
                                   in_measure = in_measure)

  return(list(
    bm_classif = bmrout_featsel,
    bm_tasks = list(task_classif=in_task,
                    task_classif_featsel = task_featsel),
    measure_classif = in_measure,
    bm_analysis = bm_analysis,
    vimp = vimp)
  )
}

#------ selecttrain_rf -----------------
selecttrain_rf <- function(in_rf, in_task,
                           insamp_nfolds =  NULL, insamp_nevals = NULL) {
  #Prepare autotuner for full training
  lrn_autotuner <- in_rf$clone()$learners$learner[[1]]
  subbm <- in_rf$clone()$filter(task_ids = in_task$id)

  if (!is.null(insamp_nfolds)) {
    lrn_autotuner$instance_args$resampling$param_set$values$folds <- insamp_nfolds
  }

  if (!is.null(insamp_nevals)) {
    lrn_autotuner$instance_args$terminator$param_set$values$n_evals <- insamp_nevals
  }

  #Return outer sampling object for selected model (or list of outer sampling objects)
  uhashes <- unique(as.data.table(subbm)$uhash)
  if (length(uhashes) == 1) {
    outer_resampling_output <- subbm$resample_result(uhash=uhashes)
  } else {
    outer_resampling_output <- lapply(uhashes, function(x) {
      subbm$resample_result(uhash=x)
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

#------ ggmisclass -----------------
ggmisclass <-  function(in_predictions=NULL, in_rftuned=NULL, spatial_rsp=FALSE) {
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

#------ ggpd -----------------
ggpd <- function (in_rftuned, in_predvars, colnums, ngrid, nodupli=T,
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
  varlabels <- setkey(in_predvars[, .(varcode, varname)], varcode)[
    levels(predmelt_num$variable)] %>%
    .[is.na(varname), varname := varcode]
  predmelt_num[, variable := factor(variable, labels = varlabels$varname)]

  varstoplot <- merge(in_predvars[, .(varcode, varname)],
                      data.table(varcode =
                                   c('totalYears_kept', 'yearskeptratio',
                                     'mDur', 'mFreq', 'station_river_distance',
                                     'UPLAND_SKM', 'ORD_STRA',
                                     'dis_m3_pyr', 'dor_pc_pva',
                                     'cmi_ix_uyr','ari_ix_uav')),
                      all.x =F, all.y=T)
  plotdt <- predmelt_num[variable %in% varstoplot$varcode |
                           variable %in% varstoplot$varname ,]

  uncertainty_numplot <-
    ggplot(plotdt, aes(x=value, y=preduncert, color=intermittent)) +
    geom_rect(xmin=-Inf, xmax=Inf, ymin=-0.5, ymax=0.5,
              fill='#d9d9d9', color='#d9d9d9') +
    geom_point(alpha = 1/4) +
    geom_hline(yintercept=0, alpha=1/2) +
    geom_smooth(method='gam', formula = y ~ s(x, k=3)) +
    annotate("text", x = Inf-5, y = 0.5, angle = 90,
             label = "Pred:Int, Obs:Per",
             color = '#1f78b4') +
    annotate("text", x = Inf-5, y = -0.5, angle = 90,
             label = "Pred:Per, Obs:Int",
             color = '#ff7f00') +
    scale_color_manual(values=c('#1f78b4', '#ff7f00'),
                       name='Observed regime',
                       labels = c('Perennial', 'Intermittent')) +
    #scale_x_sqrt(expand=c(0,0)) +
    coord_cartesian(clip='off') +
    facet_wrap(~variable, scales='free', labeller=label_value) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.1))


  #Plot categorical variables
  predmelt_cat <- predattri[, c('GRDC_NO', 'intermittent', 'preduncert',
                                'ENDORHEIC', 'clz_cl_cmj'), with=F] %>%
    melt(id.vars=c('GRDC_NO', 'intermittent', 'preduncert'))

  uncertainty_catplot <-
    ggplot(predmelt_cat, aes(x=as.factor(value), y=preduncert,
                             fill=intermittent, color=intermittent)) +
    geom_rect(xmin=-Inf, xmax=Inf, ymin=-0.5, ymax=0.5,
              fill='#d9d9d9', color='#d9d9d9', alpha=1/2) +
    #geom_boxplot(alpha = 0.75) +
    geom_violin(alpha=0.75, color=NA) +
    geom_hline(yintercept=0, alpha=1/2) +
    coord_cartesian(clip='off') +
    scale_fill_manual(values=c('#1f78b4', '#ff7f00'),
                      name='Observed regime',
                      labels = c('Perennial', 'Intermittent')) +
    scale_color_manual(values=c('#175885', '#9e3f00'),
                       name='Observed regime',
                       labels = c('Perennial', 'Intermittent')) +
    facet_wrap(~variable, scales='free', labeller=label_value) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.1))

  return(list(uncertainty_numplot=uncertainty_numplot,
              uncertainty_catplot=uncertainty_catplot))
}

#------ write_preds -----------------
write_preds <- function(in_filestructure, in_gaugep, in_gaugestats, in_rftuned,
                        in_predvars) {
  # ---- Output GRDC predictions as points ----
  in_gaugestats[!is.na(cly_pc_cav),  #SHould change that to not have to adjust subselection in multiple spots
                IRpredprob := in_rftuned$rf_inner$predict(in_rftuned$task)$prob[,2]]

  cols_toditch<- colnames(in_gaugep)[colnames(in_gaugep) != 'GRDC_NO']


  out_gaugep <- merge(in_gaugep,
                      in_gaugestats[, !cols_toditch, with=F], by='GRDC_NO')

  st_write(obj=out_gaugep,
           dsn=in_filestructure['out_gauge'],
           driver = 'gpkg',
           delete_dsn=T)

  # ---- Generate predictions to map on river network ----
  cols_tokeep <-  c("HYRIV_ID", in_predvars[!is.na(ID),varcode],
                    'ele_mt_cav','ele_mt_uav', 'gwt_cm_cav',
                    paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
                    paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
                    paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0)))

  riveratlas <- fread_cols(file_name=in_filestructure['in_riveratlas'],
                           cols_tokeep = cols_tokeep) %>%
    comp_derivedvar

  #Predict model (should take ~10-15 minutes for 9M rows x 40 cols — chunk it up by climate zone to avoid memory errors)
  for (clz in unique(riveratlas$clz_cl_cmj)) {
    print(clz)
    tic()
    riveratlas[!is.na(cly_pc_cav) & !is.na(cly_pc_uav) & clz_cl_cmj == clz,
               predbasic800 := in_rftuned$rf_inner$predict_newdata(as.data.frame(.SD))$prob[,2],]
    toc()
  }

  riveratlas[, predbasic800cat := ifelse(predbasic800>=0.5, 1, 0)]
  fwrite(riveratlas[, c('HYRIV_ID', 'predbasic800', 'predbasic800cat'), with=F],
         in_filestructure['out_riveratlas'])

  # --------- Return data for plotting ------------------------
  return(out_gaugep)
}

########## END #############
