% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{read_GRDCgauged_paths}
\alias{read_GRDCgauged_paths}
\title{Read file paths to streamflow data from GRDC gauging stations}
\usage{
read_GRDCgauged_paths(inp_GRDCgaugedir, in_gaugep)
}
\arguments{
\item{inp_GRDCgaugedir}{path to directory containing streamflow data GRDC standard files.}

\item{in_gaugep}{table containing column named \code{GRDC_NO} with the
gauge IDs that will be used to generate file path.}
}
\value{
vector of paths to GRDC-formatted streamflow time series tables, assuming
that files are called "GRDC_NO.txt", GRDC_NO being replaced with a 7-digit integer.
}
\description{
Based on selection of gauges, create a list of paths to streamflow data
associated with gauges.
}
