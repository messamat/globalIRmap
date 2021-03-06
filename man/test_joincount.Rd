% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{test_joincount}
\alias{test_joincount}
\title{Join count test}
\usage{
test_joincount(in_gauges)
}
\arguments{
\item{in_gauges}{data.table of gauging stations data, including reference and
predicted flow intermittence, hydrobasin membership, and WGS84 coordinates.
Formatted internally in \link{map_basinBACC}.}
}
\value{
data.table with the p-value and standard deviate of the join count test
for reference and predicted flow intermittence at gauging stations in each
HydroBASINS level 3.
}
\description{
Test auto-correlation in observed and predicted flow intermittence class among
gauging stations using a "join count test" for each HydroBASINS level 3.
}
\details{
For each river basin that included both IRES and perennial stations and contained at least 20
gauging stations, we tested whether spatial predictions of intermittence differed further from a random
spatial distribution than the observed patterns. We did so in the following steps:
\itemize{
\item We measured the degree of clustering separately for the observed and predicted
flow intermittence class of gauging stations — by computing the join-count statistics
(Cliff & Ord, 1981) based on four nearest neighbors (see Salima & de Bellefon, 2018 for an example implementation).
\item We assessed whether the predicted spatial distribution of intermittence differed more from what
would be expected by chance (i.e., a random distribution) than the observed distribution. This
assessment was based on the standard score between the estimated join-count statistics and the joincount
statistics that would be obtained based on a random spatial distribution of flow intermittence
classes among the stations, using 1000 permutations.
}

The join-count statistics and permutations were computed with the spatial-cross validation predictions,
using the joincount.mc function from the spdep package (Bivand et al., 2009).
}
