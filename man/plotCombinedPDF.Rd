 \name{plotCombinedPDF}
 \alias{plotCombinedPDF}
 \title{Plot combined PDF for an individual pathway}
 \description{
   A function for plotting out the pdfs for an indivdiual pathway in the QScomb object.
 }
 \usage{
plotCombinedPDF(QScomb,
                path.index=1,
                zeroLine=TRUE,
                comb.lwd=3,path.lwd=1,
                comb.col=par("col"),
                path.col=NULL,
                legend=FALSE,legend.labs=NULL,
                add=FALSE,
                xlim=NULL,ylim=NULL,
                xlab=NULL,ylab=NULL,
                main=NULL,
                type="l",
                \dots)
 }
 \arguments{
    \item{QScomb}{\code{QScomb} object}
    \item{path.index}{index describing which pathway to plot. Can either be numeric or a string containing the names of the pathway to plot. Must be of length 1}
    \item{zeroLine}{a logical indicating whether to include a vertical line at 0.}
    \item{comb.lwd}{the lwd for the combined PDF curve.}
    \item{path.lwd}{the lwd for the individual dataset curves. Can be a vector of the same length as QScomb$QSlist specifying the lwd of each individual curve.}
    \item{comb.col}{the color of the combined PDF curve.}
    \item{path.col}{the color of the individual dataset curves. Can be a vector of the same length as QScomb$QSlist specifying the color of each individual curve.}
    \item{legend}{boolean; should a legend be added to the plot? Can also be a character vector specifying the location of the legend (i.e "topleft" vs "topright", etc. See \code{\link[graphics]{legend}} for details)}
    \item{legend.labs}{character vector; the names of each dataset in QScomb$QSlist. If not provided, the function will attempt to pull the labels from QScomb$QSlist, or will use default names if those are undefined.}
    \item{add, xlim, ylim, xlab, ylab, main, type, \dots}{parameters to be passed on to \code{plot}}

 }
 \details{
This function uses the data produced by combinePDFs to plot both the individual dataset PDFs and the combined PDF for a single pathway. By default, \code{plotCombinedPDF} will plot the PDFs for the first pathway in QScomb, but this behavior can be controlled by the \code{path.index} parameter. 
  
 }
 
 \examples{
  ##create 5 example datasets of different sizes
  esets = lapply(1:5, function(i){
    n = 20 + i*5
    eset = matrix(rnorm(500*n),500,n, dimnames=list(1:500,1:n))
    
    labels = c(rep("A",10+5*floor(i/2)),
               rep("B",10+5*ceiling(i/2))
              )
               
    ##genes 1:30 are differentially expressed
    eset[1:30, labels=="B"] = eset[1:30, labels=="B"] + rnorm(30,rnorm(1,0.5,0.5),1)
    
    return(list(eset=eset, labels=labels))
  })
    
  ##gene sets
  geneSets = list(diff.set=1:30, baseline.set=31:60)
  
  ##Run qusage on each dataset
  set.results = lapply(esets, function(dat){
    qusage(dat$eset, dat$labels, "B-A", geneSets)
  })
  
  ##run the combinePDFs function
  combined = combinePDFs(set.results)
  
  ##plot the combined PDF result for "diff.set"
  plotCombinedPDF(combined, path.index="diff.set")
  
 }





