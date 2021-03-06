% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{readformatGRDC}
\alias{readformatGRDC}
\title{Read and pre-format GRDC data}
\usage{
readformatGRDC(path)
}
\arguments{
\item{path}{(character) path to the text file of daily discharge data in
standard GRDC format.}
}
\value{
\link[data.table]{data.table} of daily discharge data with additional columns
}
\description{
Reads text file of daily discharge data for a single GRDC station.
Creates columns for year, month, and date of last non-zero flow day +
computes yearly number of days of missing data
}
