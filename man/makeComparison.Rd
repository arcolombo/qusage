 \name{makeComparison}
 \alias{makeComparison}
 \title{Compare Genes Between Two Groups}
 \description{
   A function to calculate comparisons between groups in a dataset.
 }
 \usage{
 makeComparison(eset, labels, contrast, pairVector=NULL, 
                var.equal = FALSE, bayesEstimation = TRUE, 
                min.variance.factor=10^-8)
 }
 \arguments{
    \item{eset}{An objet of class \link[Biobase]{ExpressionSet} containing log normalized expression data (as created by the affy and lumi packages), OR a matrix of log2(expression values), with rows of features and columns of samples}
    \item{labels}{Vector of labels representing each column of eset}
    \item{contrast}{A string describing which of the groups in \var{labels} we want to compare. This is usually of the form \sQuote{trt-ctrl}, where \sQuote{trt} and \sQuote{ctrl} are groups represented in \var{labels}.}
    \item{pairVector}{A vector of factors (usually just 1,2,3,etc.) indicating which samples are paired. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.}
    \item{var.equal}{A logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.}
    \item{bayesEstimation}{A logical variable. If true, use a bayesian framework to estimate the standard deviation (via limma's eBayes function).}
    \item{min.variance.factor}{A factor to add to the SDs to ensure that none are equal to 0. Only used if var.equal==FALSE or bayesEstimation==FALSE.}

 }
 \details{
  This function is the first step in the qusage algorithm.  It defines the t-distributions for each gene in the input dataset by calculating the fold change and standard deviation between two groups of samples.
  
  There are two primary methods to compare two groups of data, based on whether variances of the genes in the two groups should be considered equal (as specified by the the parameter var.equal).  If var.equal=F, the t-distributions are estimated using a Welch's formalism, which is implemented internally.  Else, the LIMMA package is used to calculate the t-distribution of each gene using a pooled formalism.
  
  A note on var.equal: LIMMA's linear model function can only be run when assuming equal variances.  If var.equal==TRUE, then a linear model will be created on the entire dataset at once.  One benefit of using LIMMA's pooled variance calculation is that the linear models allow for more complicated comparisons (e.g. "(A+B)-C" or similar).  This may be of interest to some users, but in order to do this, you must assume equal variances between all groups.
  
  One caveat regarding paired samples: LIMMA can not fit a linear model when the paired samples are convoluted with the groups (e.g. one set of paired (trt vs mock) samples in patients with  disease, combined with a set of paired samples from healthy controls).  If var.equal==TRUE, these groups must be run separately to correctly fit the model (e.g. run disease first, then healthy controls).
 }
 \value{
   A \link[=QSarray-class]{QSarray} object.
 } 
 \examples{
 
  ##create example data
  eset = matrix(rnorm(500*20),500,20, dimnames=list(1:500,1:20))
  labels = c(rep("A",10),rep("B",10))
   
  ##first 30 genes are differentially expressed
  eset[1:30, labels=="B"] = eset[1:30, labels=="B"] + 1
   
  ##compare the two groups
  results = makeComparison(eset, labels, "B-A")
  
## Paired Samples
  ##Group A and group B are two samples from the same set of 10 patients
  pairVector = c(1:10,1:10)
  results.paired = makeComparison(eset, labels,"B-A",pairVector=pairVector)
  
 }





