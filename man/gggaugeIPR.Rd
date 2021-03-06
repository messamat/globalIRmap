% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{gggaugeIPR}
\alias{gggaugeIPR}
\title{Plots of Intermittence Prediction Residuals (IPR)}
\usage{
gggaugeIPR(
  in_gpredsdt,
  in_predvars,
  spatial_rsp = FALSE,
  interthresh = 0.5,
  yearthresh
)
}
\arguments{
\item{in_gpredsdt}{model predictions and environmental characteristics of gauging stations.
Here, output from \link{bind_gaugepreds}.}

\item{in_predvars}{data.table of predictor variable codes, names and attributes.
Output from \link{selectformat_predvars}.}

\item{spatial_rsp}{(boolean) whether to use outputs from spatial (TRUE) or
non-spatial (FALSE) cross-validation.}

\item{interthresh}{(numeric) flow intermittence probability threshold above which
to classify records as non-perennial}

\item{yearthresh}{(integer) minimum year from which to analyze/consider discharge record.}
}
\value{
gtable of plots
}
\description{
Plots relating Intermittence Prediction Residuals (IPR) for gauging stations
to a variety of metadata for the gauging stations including the number of years
on record, the average flow intermittence duration, the aridity in the catchment
of the station, etc.
}
\details{
this function is used to produce Extended Data Fig. 2b-e. Not all
plots were used in the final figure though.
}
