% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/011_class_func.R
\docType{methods}
\name{FromCountsObj}
\alias{FromCountsObj}
\alias{FromCountsObj,FromCountsObj-method}
\title{description of function FromCountsObj}
\usage{
FromCountsObj(dat, type = "NGSexpressionSet", outpath = "",
  name = "WorkingSet")
}
\arguments{
\item{dat}{the counts object}

\item{type}{the classnam 'StefsnExpressionSet', 'NGSexpressionSet' (default) or 'SingleCellNGS'}

\item{outpath}{the outpath for this object}

\item{name}{the name of the new obejct}
}
\description{
create a new StefansExpressionSet object based on a DEseq counts object
}
