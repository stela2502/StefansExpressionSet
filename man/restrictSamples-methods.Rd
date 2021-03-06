% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_restrictSamples.R
\docType{methods}
\name{restrictSamples}
\alias{restrictSamples}
\alias{restrictSamples,StefansExpressionSet-method}
\title{description of function restrictSamples}
\usage{
restrictSamples(x, column = "Analysis", value = NULL, name = "newSet",
  mode = "equals")
}
\arguments{
\item{x}{the StefansExpressionSet object}

\item{column}{which column to analyze in the samples table}

\item{value}{which value to take as 'cut off'}

\item{name}{the name of the new StefansExpressionSet}

\item{mode}{one of 'less', 'more', 'onlyless', 'equals' or 'grep'}
}
\description{
Drop the samples, that have been selected!
}
