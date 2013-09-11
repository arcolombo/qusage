 \name{plotCIsGenes}
 \alias{plotCIsGenes}
 \title{Plot Gene Mean and Confidence Intervals}
 \description{
   Functions for plotting the mean and confidence intervals of the genes in a pathway.
 }
 \usage{
plotCIsGenes(QSarray, 
             path.index=1,
             gene.list=NULL,
             sort.by=NULL, 
             lowerBound=0.025,
             upperBound=1-lowerBound,
             asBand=FALSE,
                   
             col=NULL,                   
             addGrid=TRUE,
             x.labels=NULL, 
             cex.xaxis=1,
             shift=0.0,
             pathwayCI=c("band","bar","none"),
             meanCol=4,
                   
             add=FALSE,
             ylim=NULL, xlim=NULL,
             ylab=NULL, xlab=NULL,
             main=NULL,
             sub=NULL,
             \dots
             )
 }
 \arguments{
    \item{QSarray}{QSarray object}
    \item{path.index}{which pathway to plot. Can either be numeric or a character vector containing the names of the pathways to plot. Must be of length 1}
    \item{gene.list}{Character vector specifying the genes in the gene set to be plotted. If \code{sort.by='none'}, the order of these genes will be used. NAs are accepted.}
    \item{sort.by}{one of c(mean,p,none), specifying how to order the genes. If NULL and gene.list is provided, default is "none", else, default is "mean".}
    \item{lowerBound, upperBound}{numeric indicating the lower and upper bounds of the confidence intervals. Default is for a 95\% confidence interval.}
    \item{asBand}{logical indicating if CIs should be plotted as a grey band or as arrows}
    \item{col}{an optional vector indicating the color for the points.}
    \item{addGrid}{Should guiding dashed lines be plotted?}
    \item{x.labels}{character vector indicating the names of the genes to be plotted along the x-axis. By default, will use the names stored in QSarray, or \code{gene.list}, if specified.}
    \item{cex.xaxis}{set cex parameter seperately for x axis label}
    \item{shift}{a number between 0 and 1 decribing the amount to shift points with respects to the guiding lines and axis labels. Useful when add=TRUE}
    \item{pathwayCI}{A string, one of "band", "bar", or "none", determining whether to add the confidence interval for the gene set PDF to the plot. By default ("band"), a band will be plotted behind the bars for the individual genes. If "bar" is specificied, another error bar will be added before the genes' error bars. To suppress the plotting of the pathway band, specify \code{pathwayCI="none"}.}
    \item{meanCol}{color for the line indicating the mean of the pathway. Only used if \code{pathwayCI} is either \code{'band'} or \code{'bar'}}
    \item{add}{logical indicating whether a new plot should be created. If \code{FALSE}, a new plot will be generated.}
    \item{xlim, ylim, xlab, ylab, main, sub,\dots}{parameters to be passed on to \code{plot}. If \code{NULL}, defaults will be used.}

 }
 \details{
This function uses the data produced by \code{\link{makeComparison}} to plot the means and confidence intervals of the genes in an individual gene set. To only plot the means and CIs of a subset of the genes in a pathway, a list of the genes to be plotted can be specified using the \code{gene.list} parameter. By default, the genes will be ordered by decreasing mean, and the 95\% confidence intervals of each point will be added. To specify a different order, \code{sort.by} must be set to \code{"none"}, and the order specified by \code{gene.list} will be used.

The mean of the overall pathway will automatically be added as a dashed line (with color specified by \code{meanCol}), but information on the confidence interval of the aggregated pathway can optionally be plotted as well. If \code{pathwayCI} is set to either "band" or "bar", the mean and CI of the gene set will be added to the plot. Specifying "band" will add the CI as a band behind the individual points, whereas "bar" will add an additional point at the left side of the plot with the mean and CI of the pathway itself. 

 }
 
 \examples{
##create example data
   eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
   labels = c(rep("A",10),rep("B",10))
   
##first 30 genes are differentially expressed for the 2 vs. 1 comparison
   diffSet = 1:30
   eset[diffSet, labels=="B"] = eset[diffSet, labels=="B"] + 1
   
#a second gene set of non-D.E. genes
   normSet = 31:60

   geneSets = list(diffSet=diffSet, normSet=normSet)

##calculate qusage results
   results = qusage(eset,labels, "B-A", geneSets)
   
##Plot gene data from first gene set
   plotCIsGenes(results, path.index=1)
   
##Add a bar to represent the differential expression of the gene set
   plotCIsGenes(results, path.index=1, pathwayCI="bar")
 }





