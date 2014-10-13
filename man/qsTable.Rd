 \name{qsTable}
 \alias{qsTable}
 \title{Summary of QSarray Results}
 \description{
   Print a table with a summary of the information on the most significant gene sets in QSarray.
 }
 \usage{
qsTable(QSarray, number=20, sort.by=c("fdr","p","logFC","none"))
 }
 \arguments{
    \item{QSarray}{A QSarray object}
    \item{number}{The number of gene sets to include in the table}
    \item{sort.by}{the metric to be used to sort the gene sets in QSarray. One of \code{c("fdr","p","logFC","none")}. You can specify just the initial letter. }
 }
 \details{
 This method will return a table with a summary of the results of qusage. 
 }
 
 \value{
   A data frame containing the following columns:
   \itemize{
     \item \code{pathway.name} - The name of the pathway
     \item \code{log.fold.change} - Average log2 fold change value of the genes in the pathway
     \item \code{p.Value} - The p-value for the gene set, as calculated using \code{pdf.pVal}
     \item \code{FDR} - The Benjamini-Hochberg False Discovery rate. Calculated using R's built-in \code{p.adjust} method.
   }
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
   
  qsTable(results)
  
##show the first 5 sets, sorted by log fold change
  qsTable(results, number=5, sort.by="l")
 }
