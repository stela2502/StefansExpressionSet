% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Grouping.addCluster.R
\docType{methods}
\name{addCluster}
\alias{addCluster}
\alias{addCluster,Clusters-method}
\title{description of function addCluster}
\usage{
addCluster(x, userGroups, col = "groupID", what = "cells")
}
\arguments{
\item{x}{the StefansExpressionSet object}

\item{userGroups}{the userGroups table}

\item{col}{the column in the userGroups table that contains the group ids default='groupID'}

\item{what}{one of 'cells' or 'genes'. Where to add the group.  default= 'cells'}
}
\description{
this function add new column to samples and annotation data NOT checking the column order
}
