 \name{plotDensityCurves}
 \alias{plotDensityCurves}
 \title{Plot gene set PDFs}
 \description{
   A function for plotting out the pdfs of a set of pathways.
 }
 \usage{
plotDensityCurves(QSarray,
                  path.index=1:numPathways(QSarray),
                  zeroLine=TRUE,
                  addVIF=!is.null(QSarray$vif),
                  col=NULL,
                  plot=TRUE,
                  add=FALSE,
                  xlim=NULL,ylim=NULL,
                  xlab=NULL,ylab=NULL,
                  type="l",
                  \dots)
 }
 \arguments{
    \item{QSarray}{QSarray object}
    \item{path.index}{vector describing which pathways to plot. Can either be numeric or a character vector containing the names of the pathways to plot.}
    \item{zeroLine}{a logical indicating whether to include a vertical line at 0.}
    \item{addVIF}{a logical indicating whether the VIF should be used to calculate the variance of the pathway.}
    \item{col}{the color of the curves. Can be a vector of the same length as path.index specifying the color of each individual curve.}
    \item{plot}{Logical indicating whether to create the plot. If \code{FALSE}, only the coordinates for the plot will be returned, and no new plot will be created.}
    \item{add, xlim, ylim, xlab, ylab, \dots}{parameters to be passed on to \code{plot}}

 }
 \details{
This function uses the data produced by aggregateGeneSet to plot the PDFs of the pathways in QSarray. By default, \code{plotDensityCurves} will plot a curve for each pathway in the QSarray pathway, but this behavior can be controlled by the \code{path.index} parameter. For the best plots, it is suggested that you limit the number of curves plotted to below ten. 
  
 }
 \value{
   Invisibly returns a list of the same length as path.index, where each entry is a matrix of x- and y-coordinates for that pathway.
 }
 
 \examples{
   ##create example data
  eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
  labels = c(rep("A",10),rep("B",10))
    
  ##genes 1:30 are differentially expressed
  eset[1:30, labels=="B"] = eset[1:30, labels=="B"] + 1
  geneSets = list(diff.set=1:30, base.set=31:60)
  
  ##Run qusage
  set.results = qusage(eset, labels, "B-A", geneSets)
 
  ##Plot results
  plotDensityCurves(set.results)
  
  ##plot just the first curve with a different color
  plotDensityCurves(set.results, path.index=1, col=2, lwd=2)
  
  ##plot the CDFs of the curves
  coords = plotDensityCurves(set.results, plot=FALSE)
  plot(0, type="n", xlim=c(-1,2),ylim=c(0,1),xlab="x",ylab="CDF")
  for(i in 1:length(coords)){
    points = coords[[i]]
    x = points$x
    y = cumsum(points$y)/sum(points$y)
    lines(x,y,col=i)
  }
  
 }





