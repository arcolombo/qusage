 \name{plotGeneSetDistributions}
 \alias{plotGeneSetDistributions}
 \title{Plot gene and gene set PDFs}
 \description{
   A function for plotting out the pdfs of all the genes in a gene set
 }
 \usage{
plotGeneSetDistributions(QSarray1, QSarray2=NULL, 
                         path.index=1, 
                         colorScheme="sdHeat",
                         alpha=1,
                         normalizePeaks=FALSE,
                         
                         addBarcode=TRUE,
                         barcode.col=NULL,
                         barcode.hei=0.2,
                         
                         groupLabel=NULL,
                         labelLoc="left",
                         xlab="Activity", 
                         ylab=NA, 
                         main=NULL,
                         lwds=c(1,3),
                         cex=1,
                         \dots)
 }
 \arguments{
    \item{QSarray1, QSarray2}{QSarray objects containing PDFs of a gene set}
    \item{path.index}{either an integer between 1 and numPathways(QSarray), or the name of the pathway to retrieve. This can be of length 1 or 2 to specify different gene sets for the top and bottom plot (see details)}
    \item{colorScheme}{This parameter specifies the color scheme to be used when plotting the individual gene PDFs. This can either be one of c("rainbow", "sdHeat") for a customized color scheme, or a vector of colors of the same length as the gene set. See the details section for more information.}
    \item{alpha}{numeric value between 0 and 1 specifying the alpha channel for the individual gene curves. Only used if \code{colorScheme} is set to one of "rainbow" or "sdHeat"}
    \item{normalizePeaks}{logical indicating whether curve heights will be normalized to the same value.}
    \item{addBarcode}{logical indicating whether a barcode-style plot should be added below the PDFs representing the means activity of each individual gene.}
    \item{barcode.col}{The color used for the bars of the barcode plot. Can be a vector of colors, or a single color which is repeated for each bar in the plot.}
    \item{barcode.hei}{a numeric value specifying the height of the barcode plot relative to the size of the PDF plot.}
    \item{groupLabel}{Vector of labels for the individual plots. If left blank, labels will be generated automatically.}
    \item{labelLoc}{vector of length 1 or 2 determining the location on the plot of where to put the label. One of "left","center", or "right"}
    \item{lwds}{a numeric vector of length 2 specifying the lwd parameters for the gene and gene set curves, respectively.}
    \item{xlab, ylab, main, cex, ...}{parameters to be passed on to \code{plot}}

 }
 \details{
  The \code{plotGeneSetDistribution} function is designed to provide a quick and intuitive look at how individual genes contribute to the overall expression of a gene set. This function plots the PDFs of each individual gene in a gene set alongside the convoluted PDF of those genes. In addition, a barcode plot representing the location of the mean fold change of each individual gene is added by default below the plot. The appearance of the curves can be controlled by the \code{colorScheme} and \code{alpha} parameters, and the barcode plot by \code{addBarcode}, \code{barcode.col}, and \code{barcode.hei}.
  
  The default \code{colorScheme}, sdHeat, will automatically color-code the gene PDFs by their standard deviations, with hotter colors being used for smaller standard deviations. This, along with \code{colorScheme="rainbow"}, are the only automatic color schemes, but \code{colorScheme} also accepts custom colors. This can be a vector of colors in any format accepted by \code{par(col)}. If the vector provided is shorter than the number of genes in the gene set, the vector will be repeated. NOTE: The order that the colors are used in is not the same as the order of genes in the original gene set. All gene sets are reordered when they are stored in the \code{QSarray$pathways} slot, and the vector provided to \code{colorScheme} will be used in this order. This also applies to any colors provided to \code{barcode.col}
  
  By default, the first pathway in the QSarray object will be plotted. If you wish to change this parameter, you can provide an alternatve pathway using the \code{path.index} parameter. This can either be an integer between 1 and \code{numPathways(QSarray1)}, or it can be a string representing the name of the pathway. 
  
  The \code{plotGeneSetDistribution} function can also be used to compare the results from two different pathways or datasets. In order to analyze two different pathways from the same QSarray object, you can provide a \code{path.index} parameter of length 2 representing the two pathways to be compared. Alternatively, a separate QSarray object can be provided as the parameter \code{QSarray2}, and the second plot will be drawn from this object. If \code{QSarray2} is provided and \code{path.index} is of length 2, the second \code{path.index} will be drawn from \code{QSarray2}.
  
  
 }
 
 \examples{
 
 ##create example data
   eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
   labels = c(rep("A1",5),rep("A2",5),rep("B1",5),rep("B2",5))
   
##first 30 genes are differentially expressed much more strongly in group "B" than in group "A"
   geneSet = 1:30
   eset[geneSet, labels=="A2"] = eset[geneSet, labels=="A2"] + 1
   eset[geneSet, labels=="B2"] = eset[geneSet, labels=="B2"] + 2
   
##calculate qusage results
   A.results = qusage(eset,labels, "A2-A1", geneSet)
   B.results = qusage(eset,labels, "B2-B1", geneSet)

##plot the gene set distribution for group A and group B side-by-side
  plotGeneSetDistributions(A.results,B.results)
  
##add labels to the right side of the plots
  plotGeneSetDistributions(A.results,B.results,groupLabel = c("A2-A1", "B2-B1"), labelLoc="right")
  
##change the colors of the curves
  plotGeneSetDistributions(A.results,B.results, colorScheme="rainbow")
 }





