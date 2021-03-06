% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.gg.heatmap.list.R
\docType{methods}
\name{gg.heatmap.list}
\alias{gg.heatmap.list}
\alias{gg.heatmap.list,StefansExpressionSet-method}
\title{description of function gg.heatmap.list}
\usage{
gg.heatmap.list(dat, glist = NULL, colrs = NULL, groupCol = "GroupID",
  colCol = NULL, lowest = "gray")
}
\arguments{
\item{dat}{the StefansExpressionSet object}

\item{glist}{a list of probesets to plot (or all)}

\item{colrs}{a list of colors for the sample level boxes (or rainbow colors)}

\item{groupCol}{the column group in the samples table that contains the grouping strings}

\item{colCol}{the column group in the samples table that contains the color groups}

\item{lowest}{the color for the lowest expression value (default NULL)}
}
\description{
uses ggplot2 to plot heatmaps
}
