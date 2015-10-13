 \name{pVal}
 \alias{pdf.pVal}
 \alias{oneWay.pVal}
 \alias{twoCurve.pVal}
 \alias{twoWay.pVal}
 \title{Calculate p-values for gene set activity}
 \description{
   Methods for calculating the significance of gene set activity, compared either to a null hypothesis (\code{pdf.pVal}), or to a separate PDF (\code{twoCurve.pVal}).
 }
\usage{
pdf.pVal(QSarray, alternative=c("two.sided","less","greater"),
         direction=FALSE, addVIF=!is.null(QSarray$vif), selfContained=TRUE)
            
twoCurve.pVal(grp1, grp2, path.index1 = 1:numPathways(grp1), path.index2 = 1:numPathways(grp2), 
              alternative=c("two.sided","less","greater"), direction=FALSE, 
              addVIF=!(is.null(grp1$vif) | is.null(grp2$vif)))
}
\arguments{                                              
    \item{QSarray, grp1, grp2}{A QSarray object as output by \link{qusage} (or \link{aggregateGeneSet})}
    \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.}
    \item{direction}{a logical indicating whether the p-values should be signed (i.e. negative fold changes return negative p-values). Ignored if \code{alternative!="two.sided"}.}
    \item{addVIF}{a logical indicating whether to use the VIF when calculating the variance}
    \item{selfContained}{a logical indicating whether the test should be self-contained or competitive. See details for more information.}
    \item{path.index1, path.index2}{numeric vectors indicating which gene sets in grp1 to compare to grp2. The length of path.index1 and path.index2 must match.}
}
\details{
   The pVal functions are designed to estimate the level of significance for the gene set activity cacluated using qusage. Because the QSarray object contains gene set information stored as a Probability Density Function (PDF), we can determine significance of an individual gene set using the \code{pdf.pVal} function by comparing the PDF to our null hypothesis (zero by default. See below).  If \code{alternative="greater"}, \code{pdf.pVal} tests whether the fold change of the gene set is greater than the null mean, and the p-value is calcuated based on the proportion of the lower tail of the PDF which is below the null hypothesis.
   
   There are two options for the null hypothesis in this method, controlled by the logical parameter "\code{selfContained}". By default, \code{pdf.pVal} performs a self-contained test, where the null hypothesis is that the mean fold change is 0. If \code{selfContained=FALSE} is specified, \code{pdf.pVal} instead performs a competitive test, where the null hypothesis is the mean fold change of all genes which are not in the pathway.
   
   An individual gene set's PDF can also be compared with a second PDF, created from either comparing a different set of samples or using a different gene set, using the \code{twoCurve.pVal} function. This function takes two QSarray objects, \code{grp1} and \code{grp2}, and by default compares the PDFs for each gene set in the two QSarray objects in order. However, this behavior can be controlled by the \code{path.index1} and \code{path.index2} parameters, which are numeric vectors specifying which gene sets should be compared. The two vectors must be the same length, and the first index in \code{path.index1} will be compared with the first index in \code{path.index2} and so on. 
}
 
 \value{
   A vector of p-values for each gene set in QSarray, or for each gene set specified with \code{path.index} when using \code{twoCurve.pVal}.
 }
 \examples{
##create example data
   eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
   labels = c(rep("A1",5),rep("A2",5),rep("B1",5),rep("B2",5))
   
   geneSets = list()

##first 30 genes are differentially expressed for the 2 vs. 1 comparison
   geneSets[["simple.diffSet"]] = 1:30
   eset[geneSets[[1]], labels=="A2"] = eset[geneSets[[1]], labels=="A2"] + 1
   eset[geneSets[[1]], labels=="B2"] = eset[geneSets[[1]], labels=="B2"] + 1
   
##second set of 30 genes different in only group B
   geneSets[["complex.diffSet"]] = 31:60
   eset[geneSets[[2]], labels=="B2"] = eset[geneSets[[2]], labels=="B2"] + 1
   
#a third gene set of non-D.E. genes
   geneSets[["normSet"]] = 61:90

##calculate qusage results
   A.results = qusage(eset,labels, "A2-A1", geneSets)
   B.results = qusage(eset,labels, "B2-B1", geneSets)

##calculate p-values for initial comparison
   pdf.pVal(A.results)
   pdf.pVal(B.results)
 
##compare the pdfs of the two groups
  twoCurve.pVal(A.results,B.results)
 }
