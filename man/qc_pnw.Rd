% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRmapping_functions.R
\name{qc_pnw}
\alias{qc_pnw}
\title{Quality check model estimates in U.S. Pacific Northwest}
\usage{
qc_pnw(inp_pnwresdir, in_rivpred, predcol, interthresh = 0.5, mincutoff = 0)
}
\arguments{
\item{inp_pnwresdir}{(character) full path to directory containing formatted
data for Pacific Northwest.}

\item{in_rivpred}{data.table of model predictions for global river network. Here, output from \link{netpredformat}.}

\item{predcol}{(character) name of the column for the predicted probability
of flow intermittence generated by this study (Messager et al. 2021).}

\item{interthresh}{(numeric) flow intermittence probability threshold above which
to classify records as non-perennial.}

\item{mincutoff}{(numeric) smallest river reaches to include in comparison, in terms of long-term mean annual flow (m3/s; e.g., 0.1)}
}
\value{
list with two elements:
*"plot": scatterplots analyzing Intermittence Prediction Residuals (IPR) as a function of
modeled mean annual flow, the position of the field observation along the RiverATLAS river reach,
the number of field observations at the site, and the population density in the entire upstream area
of the corresponding reach.
*"stats": summary comparison statistics.
}
\description{
Analyze performance of model predictions of flow intermittence against PROSPER field
observations of flow state across the U.S. Pacific Northwest.
}
\details{
this function is used to format data used in producing Extended Data Figure 6.
A side effect of this function is the writing of a .gpkg of river reaches for
which there are flow observations in the U.S. Pacific Northwest with quality-checking
attributes.
}
