 \name{plotCIs}
 \alias{plotCIs}
 \title{Plot Pathway Mean and Confidence Intervals}
 \description{
   Functions for plotting the mean and confidence intervals of a set of pathways.
 }
 \usage{
plotCIs(QSarray, 
        path.index=1:numPathways(QSarray),
        sort.by=c("mean","p","none"), 
        lowerBound=0.025,
        upperBound=1-lowerBound,
                   
        col=NULL,
        use.p.colors=TRUE,
        p.breaks=NULL,
        p.adjust.method = "fdr",
        addLegend=use.p.colors,
        lowerColorBar="none",
        lowerColorBar.cols=NULL,
                   
        addGrid=TRUE,
        x.labels=NULL, 
        cex.xaxis=1,
        shift=0.0,
                   
        add=FALSE,
        ylim=NULL, xlim=NULL,
        ylab=NULL, xlab=NULL,
        main=NULL,
        sub=NULL,
        type="p",
        \dots
        )
 }
 \arguments{
    \item{QSarray}{QSarray object}
    \item{path.index}{vector describing which pathways to plot. Can either be numeric or a character vector containing the names of the pathways to plot.}
    \item{sort.by}{One of c("mean","p","none") indicating the order that the pathways should be plotted in. If "none", the pathways will not be reordered, and the order specified in path.index will be maintained}
    \item{lowerBound, upperBound}{numeric indicating the lower and upper bounds of the confidence intervals. Default is for a 95\% confidence interval.}
    \item{col}{an optional vector indicating the color for the points. If \code{use.p.colors=FALSE} is specified, these colors will also be used for the error bars.}
    \item{use.p.colors}{logical indicating whether error bars should be colored based on the significance of the p-value.}
    \item{p.breaks}{a vector indicating where the breaks in the p-value color scheme should be. By default, breaks will be at 0.001, 0.005, 0.01, 0.05, & 0.1}
    \item{p.adjust.method}{The method to use to adjust the p-values. Must be one of the methods in \code{\link{p.adjust.methods}}.}
    \item{addLegend}{a logical specifying if a legend for the p-value color scheme be plotted}
    \item{lowerColorBar}{Options for plotting a color bar below each point. Automatically generated color bars have not yet been implemented, but custom bars can be created using the "lowerColorBar.cols" parameter.}
    \item{lowerColorBar.cols}{a vector of colors to plot as a bar below each point.}
    \item{addGrid}{Should guiding dashed lines be plotted?}
    \item{x.labels}{character vector of the same length as \code{path.index} giving the names of the pathways. By default, will use the names stored in QSarray.}
    \item{cex.xaxis}{set cex parameter seperately for x axis label}
    \item{shift}{a number between 0 and 1 decribing the amount to shift points with respects to the guiding lines and axis labels. Useful when \code{add=TRUE}}
    \item{add}{logical indicating whether a new plot should be created. If \code{FALSE}, a new plot will be generated.}
    \item{xlim, ylim, xlab, ylab, main, sub, type,\dots}{parameters to be passed on to \code{plot}}

 }
 \details{
This function uses the data produced by \code{aggregateGeneSet} to plot the means and confidence intervals of the gene sets in \code{QSarray}. By default, the gene sets will be ordered by decreasing mean, and the 95\% confidence intervals of each point (as calculated by \link{calcBayesCI}) will be added. To specify a different order, \code{sort.by} must be set to \code{"none"}, and the order specified by \code{path.index} will be used.

The points in the plot can be optionally color-coded by the significance of the (corrected) p-values. The p-values are adjusted using R's built in \code{\link{p.adjust}} method, which uses the \code{p.adjust.method} parameter to determine the algorithm being used. The colors of the points are based on the breaks specified in \code{p.breaks}. By default, more significant p-values will be plotted in bright red/green. If \code{use.p.colors} is specified and \code{addLegend=TRUE}, a legend describing the p-values will be added to the top left corner of the plot. Alternatively, if you want to specify the colors of the points individually, you can provide a vector of colors to the \code{col} parameter. 

The \code{plotCIs} function can also add a color bar along the bottom of the plot to provide additional information about the pathways. We are currently working on implementing various metrics which can be added automatically using the \code{lowerColorBar} parameter, but in the mean time, the bar can be added manually by providing a vector of colors the same length as \code{path.index} to the \code{lowerColorBar.cols} parameter.
  
 }
 
 \examples{
 
##create example data
  eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
  labels = c(rep("A",10),rep("B",10))
   
  geneSets = list()

##create a number of gene sets with varying levels of differential expression.
  for(i in 0:10){
    genes = ((30*i)+1):(30*(i+1))
    eset[genes,labels=="B"] = eset[genes,labels=="B"] + rnorm(1)
    
    geneSets[[paste("Set",i)]] = genes
  }
   
   
##calculate qusage results
  results = qusage(eset,labels, "B-A", geneSets)
   
##Plot gene set CIs
  plotCIs(results)
   
 }





