 \name{getXcoords}
 \alias{getXcoords}
 \title{Get the X coordinates for the points of the PDF}
 \description{
   Calculates the x-coordinates for the PDF of a given pathway.
 }
\usage{
getXcoords(QSarray, path.index=1, addVIF=!is.null(QSarray$vif))

}
\arguments{                                              
    \item{QSarray}{A QSarray object as output by \link{qusage} (or \link{aggregateGeneSet})}
    \item{path.index}{either an integer between 1 and numPathways(QSarray), or the name of the pathway to retrieve.}
    \item{addVIF}{a logical indicating whether to use the VIF when calculating the variance}
}
\details{
   The calculation of the x-coordinates for a PDF is not straightforward, and as such they are not included in the QSarray object initially. During the numerical convolution step, the gene set PDF is calculated at a number of points (equal to \code{QSarray$n.points}) over a range defined by:
   
  \code{c(path.mean - range, path.mean + range})
   
  However, the resulting PDF is actually the \emph{sum} of the individual gene PDFs, rather than the desired \emph{average} PDF. Therefore the range which is stored in the resulting QSarray is divided by the number of genes in the pathway, \code{QSarray$path.size}.
  
  In addition, the width of the PDF can be expanded by the Variance Inflation Factor (VIF), which is equivalent to multiplying the range of the x-coordinates by the \code{sqrt(VIF)}. If the parameter \code{addVIF=TRUE}, the VIF calculatd using the \code{calcVIF} method will be included in the calculation of the x-coordinates.
  
  In general, the x-coordinates for a pathway are calculated for each point n using the following formula:
  
  \deqn{x_n = (-1+\frac{2(n-1)}{N_{pts}-1}) \times r \times \sqrt{VIF} + \hat{\mu}_{path}}{x.n = ( seq(-1,1,length.out=N.points) * range * sqrt(VIF) ) + path.mean}
}
 
 \value{
   A numeric vector of length \code{QSarray$n.points}.
 }
 \examples{
 
  ##create example data
  eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
  labels = c(rep("A",10),rep("B",10))
   
  ##first 30 genes are differentially expressed
  eset[1:30, labels=="B"] = eset[1:30, labels=="B"] + 1
  geneSets = list(diff.set=1:30, base.set=31:60)
  
  ##Run qusage
  set.results = qusage(eset, labels, "B-A", geneSets)
  
  ##Plot the PDF (see also: plotDensityCurves() )
  x = getXcoords(set.results, 1)
  y = set.results$path.PDF[,1]
  plot(x,y, type="l")
 }
