\name{QSarray-class}
\docType{class}
\alias{QSarray-class}
\alias{dim,QSarray-method}
\alias{head,QSarray-method}
\alias{numFeatures}
\alias{numPathways}
\alias{numFeatures,QSarray-method}
\alias{numPathways,QSarray-method}
\alias{plot,QSarray-method}
\alias{print,QSarray-method}
\alias{summary,QSarray-method}

\title{Class \code{"QSarray"}}
\description{
   A list-based class which contains the results of running qusage. Generally created by \code{\link{qusage}} or \code{\link{makeComparison}}
}
\section{Objects from the Class}{
QSarray objects should not be created directly, but rather through the \code{\link{makeComparison}} function. They can also be created manually via a call to the \code{\link{newQSarray}} function, although this should be done by advanced users only.
}
\section{Components}{
  \code{QSarray} objects do not contain any slots (apart from \code{.Data}) but they should contain the following list components:
  \tabular{ll}{
    {\code{mean}} \tab { numeric vector containing mean fold changes for individual genes}\cr \cr 
    {\code{SD}} \tab {numeric vector of standard deviations for individual genes}\cr 
    {\code{dof}} \tab {numeric vector. Degrees of Freedom for each gene}\cr 
    {\code{var.method}} \tab {one of \code{("Welch's","Pooled")}, indicating the method used to calculate the variance}\cr 
    {\code{sd.alpha}} \tab {The factor each sd is multiplied by (either due to the min.variance.factor parameter in makeComparison or because of the Bayesian shrinkage of the SDs). This is used when calculating the VIF in order to correct for genes with 0 (or very small) variance.}\cr 
    {\code{labels}} \tab {The labels as input in makeComparisons, describing the group structure of the data.}\cr 
    {\code{pairVector}} \tab {A vector indicating which samples should be treated as pairs.}\cr 
    {\code{contrast}} \tab {A string describing which of the two groups in labels was compared.}\cr 
  }
  
  The following additional components are appended to the object by running \code{\link{aggregateGeneSet}} and \code{\link{calcVIF}}
  
  \tabular{ll}{ 
    {\code{pathways}} \tab {the list of genes in each gene set. Represented as a list of indices.}\cr 
    {\code{path.mean}} \tab {vector describing the mean fold change for each of the pathways provided to AggregateGeneSet}\cr 
    {\code{path.PDF}} \tab {Matrix describing the probability distributions for each of the pathways provided to AggregateGeneSet, where each column is a different gene set, and each row is a different point where this set was sampled. x-coordinates must be generated using \code{\link{getXcoords}}}\cr 
    {\code{path.size}} \tab {numeric vector containing the number of features in each pathway that mapped to the input data.}\cr 
    {\code{ranges}} \tab {the (uncorrected) range that all PDFs were calculated over. If the VIF is not used to correct the range, the x-coordinates of the PDF are the sequence of n.points from \code{path.mean-ranges} to \code{path.mean+ranges}}\cr 
    {\code{n.points}} \tab {The number of points that the PDF was calculated at. This is equal to the number of rows in \code{path.PDF}}\cr 
    {\code{vif}} \tab {the Variance Inflation Factor for each pathway, as calculated by \code{\link{calcVIF}}}
  }
}

\section{Methods}{
  \describe{
  
    \item{newQSarray}{ The constructor for the QSarray object. Should primarily be used internally by \code{\link{qusage}} or \code{\link{makeComparison}}. See \code{\link{newQSarray}} for additional details.}
  
    \item{numFeatures}{Returns the number of features (i.e. genes or probesets) in the dataset}
    \item{numPathways}{Returns the number of pathways provided to \code{\link{aggregateGeneSet}} }
    \item{dim}{dimensions of the QSarray object, as c(numFeatures, numPathways) }
    \item{print, head}{Prints a summarized version of all fields in the QSarray object.}
    \item{summary}{Prints a brief summary of the QSarray object.}
    \item{plot}{Plots the information stored in QSarray by either calling \code{\link{plotDensityCurves}} (if \code{numPathways} < 10) or \code{\link{plotCIs}} (if \code{numPathways} >= 10)}
    \item{qsTable}{Print a table with a summary of the information on the most significant gene sets in QSarray. See \code{\link{qsTable}} for more details.}
	 }
}
\author{
Christopher Bolen
}
\keyword{classes}
