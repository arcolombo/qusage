\name{fluExample}
\docType{data}
\alias{fluExample}
\alias{eset}
\alias{resp}
\title{Example gene expression set}
\description{
 This is a matrix containing microarray gene expression values taken from a publicly available dataset (GEO ID: GSE30550; Huang Y et al. PLoS Genet 2011). This dataset contains samples from 17 patients who were exposed to Influenza and had blood drawn approximately every 8 hours for a week. Patients were classified as either symptomatic (sx) or asymptomatic (asx) based on the severity of their symptoms.
 
 The portion of the dataset included here contains only two time points: a pre-exposure point (time 0) and a post-exposure time point (77h after exposure). \code{eset} is arranged in groups, with the pre-exposure samples in columns 1:17, and the post-exposure in columns 18:34. The patients are in the same order in both groups, and the response of each patient is defined by the \code{resp} variable, which is also contained in the \code{fluExample} data object.
}
\usage{
  resp
  eset
}
\format{\code{eset} is a matrix of gene expression measurements, with rows of genes and columns representing samples. \code{resp} is a character vector.}
\source{\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30550}}
\references{
 Huang Y et al. (PLoS Genet 2011)
}
\keyword{datasets}