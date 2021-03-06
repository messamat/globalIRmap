% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{extract_pd_nestedrf}
\alias{extract_pd_nestedrf}
\title{Extract partial dependence (and performance) from a trained RF learner.}
\usage{
extract_pd_nestedrf(
  learner_id = 1,
  in_rftuned,
  datdf,
  selcols,
  nvariate,
  ngrid
)
}
\arguments{
\item{learner_id}{(integer) Index of the outer resampling instance to be
analyzed.}

\item{in_rftuned}{\link[mlr3]{ResampleResult} from a classification RF.}

\item{datdf}{Data from the task that was used to train RF.}

\item{selcols}{Character vector of the predictor variables to analyze.}

\item{ngrid}{(integer) Number of values of the  predictor variables over
which to compute the marginal relationship.}
}
\value{
A data.table with the following columns.
\describe{
\item{value1} - value of the first predictor variable in the pair
\item{value2} - value of the second predictor variable in the pair
\item{0} - predicted probability of 0 (e.g. probability that the river
is perennial) at value1 and value2
\item{1} - predicted probability of 1 (e.g. probability that the river
is intermittent) at value1 and value2
\item{[perf]} - performance value (same for all rows), name changes
\item{var1} - name of the first predictor variable
\item{var2} - name of the second predictor variable
}
}
\description{
Computes the marginal relationship between a subset of the predictors
(here, two variables at a time) and the model’s predictions by averaging
over the marginal distribution of the compliment of this subset of
the predictors, taking in account the interaction between the chosen
predictors.
}
\details{
Also accept learners of class \link[mlr3pipelines]{GraphLearner}. \cr
Uses \link[edarf]{partial_dependence} for computing.
}
\section{Warning}{

Has only been tested on \link[mlr3learners]{mlr_learners_classif.ranger}
}

\section{documentation to-do}{

Can add an example down the line, add source.
}

\seealso{
\code{\link{weighted_vimportance_nestedrf}},
\code{\link{ggpd_bivariate}}
}
