% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{dynamic_resamplebm}
\alias{dynamic_resamplebm}
\title{Dynamic resampling using learners from benchmarking results}
\usage{
dynamic_resamplebm(
  in_task,
  in_bm,
  in_lrnid,
  in_resampling,
  type,
  inp_resdir = NULL,
  store_models = FALSE
)
}
\arguments{
\item{in_task}{\link[mlr3]{Task} to resample}

\item{in_bm}{\link[mlr3]{BenchmarkResult} from which to extract learner or name of
serialized file (.qs) on disk (accessed through \link[qs]{qread})}

\item{in_resampling}{\link[mlr3]{Resampling} cross-validation strategy}

\item{type}{(character) Type of learner 'classif' or 'regr'. Other values are not accepted.}

\item{inp_resdir}{(character) path to where qs file is located (excluding the name of the qs file)}

\item{store_models}{whether to keep the fitted model after the test set has been predicted. Set to TRUE if you want to further analyse the models or want to extract information like variable importance.}

\item{in__lrnid}{id of learner to extract from \link[mlr3]{BenchmarkResult} (e.g., "oversample.classif.ranger")}
}
\value{
\link[mlr3]{ResampleResult}
}
\description{
Run a train-predict-test routine using a learner from a \link[mlr3]{BenchmarkResult}:
Repeatedly apply learner on a training set of task to train a model,
then use the trained model to predict observations
of a test set. Training and test sets are defined by the resampling.
}
\details{
The dynamic aspect of this model is that it runs 'reset_tuning' on the fly
to make sure that the hyperparameter search space matches the task (e.g., if the number of candidate
predictor variables has been reduced, it adjusts mtry)
}
