% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialProtUtils.R
\name{calcCorrelationWithScore}
\alias{calcCorrelationWithScore}
\title{calcCorrelationWithScore}
\usage{
calcCorrelationWithScore(sce, scoreName, protVal = "prot", method = "spearman")
}
\arguments{
\item{sce}{SpatialExperiment object}

\item{scoreName:}{name of score}

\item{protVal:}{name of feature, e.g. prot or substrate}

\item{method:}{spearman or pearson}
}
\value{
SingleCellExperiment objecti
}
\description{
calculates the correlation of each
element with a numeric vector of the score, such as distance or
an immune score, puts value in RowData
}
