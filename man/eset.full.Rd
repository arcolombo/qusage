\name{fluExample}
\docType{data}
\alias{fluExample}
\alias{eset.full}
\alias{flu.meta}
\title{Example gene expression set}
\description{
 This is a matrix containing microarray gene expression values taken from a publicly available dataset (GEO ID: GSE30550; Huang Y et al. PLoS Genet 2011). This dataset contains samples from 17 patients who were exposed to Influenza and had blood drawn approximately every 8 hours for a week. Patients were classified as either symptomatic (sx) or asymptomatic (asx) based on the severity of their symptoms.
 
 The portion of the dataset included here contains a total of 15 time points, including a pre-exposure point (time 0) and points approximately every 8 hours up to 108 hours. \code{eset.full} is arranged by donor, and the information for each sample is contained in a table, \code{flu.meta}. The metadata table contains 7 columns: 
 \describe{
  \item{\code{SampleID}}{The GEO sample IDs, matching the column names of eset.full}
  \item{\code{Subject}}{the subject ID}
  \item{\code{Hours}}{The hour-post-infection that the sample was collected (stored as a factor)}
  \item{\code{Hours.Numeric}}{Same as above, but stored as a numeric vector}
  \item{\code{Condition}}{The condition of the donor, either symptomatic (sx) or asymptomatic (asx)}
  \item{\code{Gender}}{The donor gender}
  \item{\code{Age}}{Age of the donor at recruitment}
 }
}
\usage{
  eset.full
  flu.meta
}
\format{\code{eset} is a matrix of gene expression measurements, with rows of genes and columns representing samples. \code{resp} is a character vector.}
\source{\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30550}}
\references{
 Huang Y et al. (PLoS Genet 2011)
}
\keyword{datasets}