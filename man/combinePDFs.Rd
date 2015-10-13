\name{combinePDFs}
\alias{combinePDFs}
\title{
Combine PDFs from multiple QuSAGE comparisons
}
\description{
This function combines the results of multiple \link{qusage} runs into a single, joint PDF. The resulting PDFs will be the average of the PDFs from each individual QSarray object, weighted by the number of samples in each dataset.
}
\usage{
combinePDFs(QSarrayList, n.points=2^14)
}
\arguments{
  \item{QSarrayList}{A list of QSarray objects, each generated from the same list of geneSets}
  \item{n.points}{integer; the number of points at which to sample the convoluted t-distribution. For best results, this should be about 2-4 times higher than the n.points used in the individual QSarray objects.}
}
\details{
Like \link{aggregateGeneSet}, \code{combinePDFs} uses numerical convolution to calculate the combined PDFs for individual pathways, with each individual PDF weighted by the total number of samples used in the comparison. This method is useful for meta-analysis of multiple datasets, or for a meta comparison where the difference between two QuSAGE pdfs is of interest. 

The results of \code{combinePDFs} can be plotted (on a pathway-by-pathway basis) using the \code{\link{plotCombinedPDF}} function, or by simply calling "\code{plot}" on a QSarray object which contains the QSlist field.
}
\value{
A QSarray object containing a convolution of the original QSarrays. The new QSarray object will contain an additional field, \code{QSlist}, containing the input QSarrayList. 
}

\examples{
##create example data - a set of 500 genes normally distributed across 40 patients
  eset = matrix(rnorm(500*40),500,40, dimnames=list(1:500,1:40))
  labels = rep(c("A","B","C","D"),each=10)
   
##create a number of gene sets with varying levels of differential expression.
  geneSets = list()
  for(i in 0:10){
    genes = ((30*i)+1):(30*(i+1))
    eset[genes,labels=="B"] = eset[genes,labels=="B"] + 2 + rnorm(1)
    eset[genes,labels=="D"] = eset[genes,labels=="D"] + 1 + rnorm(1)
    
    geneSets[[paste("Set",i)]] = genes
  }
   
   
##calculate qusage results
  qsList = lapply(c("B-A","D-C"), function(comparison){
    qusage(eset,labels, comparison, geneSets)
  })
  
##combine the two QSarrays
  qsComb = combinePDFs(qsList)
  
  plot(qsComb, path.index=1)
}

