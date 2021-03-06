% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/011_class_func.R
\docType{methods}
\name{StefansExpressionSet}
\alias{StefansExpressionSet}
\alias{StefansExpressionSet,data.frame-method}
\title{description of function StefansExpressionSet}
\usage{
StefansExpressionSet(dat, Samples, class = "StefansExpressionSet",
  name = "WorkingSet", namecol = NULL, namerow = "GeneID", outpath = "")
}
\arguments{
\item{dat}{data frame or matrix containing all expression data}

\item{Samples}{A sample description table}

\item{name}{The name of this object is going to be used in the output of all plots and tables - make it specific}

\item{namecol}{The samples table column that contains the (to be set) names of the samples in the data matrix}

\item{namerow}{This is the name of the gene level annotation column in the dat file to use as ID}

\item{outpath}{Where to store the output from the analysis}
}
\description{
this file contains all generic fnction for data export and ploting Create an StefansExpressionSet

object (S3) This object is mainly used for subsetting of the data and plotting @export
}
