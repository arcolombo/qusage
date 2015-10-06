\name{fluVaccine}
\docType{data}
\alias{fluVaccine}
\title{Gene expression sets from Flu Vaccine trials}
\description{
  \code{fluVaccine} is a list containing data from two separate Flu Vaccine treatment studies, \code{GSE59635} and \code{GSE59654}. The top-level list contains two entries, \code{esets} and \code{phenoData}.
  
  \code{fluVaccine$esets} is a list containing two expression matrices. Each matrix, which is labeled with the name of the GEO dataset it originates from, contains rows of gene expression levels and columns of samples.
  
  \code{fluVaccine$phenoData} is a list containing two data frames. The rows of each data frame correspond to the columns in the expression matrices, and each data frame contains the following columns:
  \describe{
    \item{subjectID}{The subject ID}
    \item{responder}{Factor describing the subject's response to flu vaccination, either "high" or "low".}
    \item{bloodDrawDay}{The day post-vaccination that the sample was taken.}
  }
}
\usage{
  fluVaccine
}
\format{\code{fluVaccine} is a list.}
\source{
  \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59635}
  \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59654}
}
\references{
 Thakar J et al., Aging, 2015
}
\keyword{datasets}