% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{readformatGSIMmon}
\alias{readformatGSIMmon}
\title{Read and pre-format GSIM monthly data}
\source{
Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
time-series indices and homogeneity assessment. Earth System Science Data,
10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
}
\usage{
readformatGSIMmon(path)
}
\arguments{
\item{path}{(character) path to the text file of monthly discharge data in
standard GSIM format (.mon file).}
}
\value{
\link[data.table]{data.table} of monthly indices, with additional
attributes including the estimated minimum number of zero-flow days "mDur_minmo"
}
\description{
Reads text file of monthly discharge indices for a single GSIM station.
Creates columns for year, month, and season (to join with seasonal indices)
Estimate the minimum possible number of monthly zero-flow days based on indices
}
\details{
compute yearly number of days of missing data

the monthly minimum number of zero flow days n0 is estimated as follows: \cr
1. If MAX == 0, then n0 = n.available
2. else if MIN7 == 0, then n0 = 7
3. else if MIN == 0, then n0 = 1,
4. else n0 = 0\preformatted{For more information on GSIM, see [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/).
}
}
