% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{set_tuning}
\alias{set_tuning}
\title{Set tuning strategy}
\usage{
set_tuning(
  in_learner,
  in_measure,
  nfeatures,
  insamp_nfolds,
  insamp_neval,
  insamp_nbatch
)
}
\arguments{
\item{in_learner}{\link[mlr3]{Learner} whose hyperparameters to tune.}

\item{in_measure}{performance measure to use for hyperparameter tuning.
The hyperparameter configuration that optimizes this measure will be selected and used for final model training.}

\item{nfeatures}{total number of predictor variables in the model.}

\item{insamp_nfolds}{number of cross-validation folders used for tuning}

\item{insamp_nbatch}{number of hyperparameter configurations to evaluate at the same time. This will dictate how many
processing cores will be used in hyperparameter tuning. The greater this number, the faster the tuning, but also the
more computing intensive. This shouldn't be set higher than the number of cores on the computer used.}

\item{insamp_neva}{number of times the cross-validation must be conducted (e.g. if 2, then twice-repeated CV)}
}
\value{
a \link[mlr3tuning]{AutoTuner}.
}
\description{
Create a \link[mlr3tuning]{AutoTuner} to set the random forest hyperparameter
tuning strategy.
}
\details{
For more information on hyperparameter tuning and a table of the range of
hyperparameters and tuning strategies used, see sections IVa and IVb in the Supplementary
Information of Messager et al. 2021 at \link{https://www.nature.com/articles/s41586-021-03565-5}.
}
