% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{ggmisclass_single}
\alias{ggmisclass_single}
\title{Single ggplot of misclassification}
\usage{
ggmisclass_single(
  in_predictions = NULL,
  in_rftuned = NULL,
  spatial_rsp = FALSE
)
}
\arguments{
\item{in_predictions}{Either:
\enumerate{
\item a \link[mlr3]{PredictionClassif} or
\item a data.table of predictions for a set of CV repetitions as formatted by
\code{\link{analyze_benchmark}}.
}}

\item{in_rftuned}{Output from \link{selecttrain_rf};
list containing inner and outer resampling results + task.}

\item{spatial_rsp}{(boolean) whether to use spatial or non-spatial
cross-validation results}
}
\value{
list containing a ggplot and the threshold for which sensitivity == specificity (numeric)
}
\description{
Plot cross-validation sensitivity, specificity, and misclassification rate
based on probability threshold for classifying a reach as non-perennial.
}
\details{
in_rftuned is only needed if in_predictions is not provided.
}
