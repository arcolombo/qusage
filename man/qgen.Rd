\name{qgen}
\alias{qgen}
\title{Run qusage while incoprating generalized least squares and linear mixed model analysis at the gene level to account for repeated measures, continous covariates, and confounder adjusting.
}
\description{
A wrapper function for the three primary steps in the qusage algorithm. The first step replaces the conventional t-test framework, with additional flexibility to incorporate general linear models from the nlme package, specifically the lme and gls function calls. 
}
\usage{
qgen(eset,design,fixed,geneSets,
     contrast.factor,contrast,
     random=NULL,correlation=NULL,design.sampleid=NULL)
}
\arguments{
  \item{eset}{A matrix of log2(expression values), with rows of features and columns of samples.}
  \item{design}{A data frame consisting of sample annotation information such as the various fixed effects that wished to be included in modeling and subject identifiers for the inclusion of random effects.It is recommended that the design file have a column which consists of the column names of \var{eset}. See \var{design.sampleid}}
  \item{fixed}{A one-sided linear formula object describing the fixed-effects part of the model. The formula should always begin with a ~ operator followed by the fixed effect terms, seperated by the + operator.}
  \item{geneSets}{Either a list of pathways to be compared, or a vector of gene names representing a single gene set. See Description for more details.}
  \item{contrast.factor}{A one sided formula indicating the factor in which the user whishes to create a contrast. ~X1 indicates that a differences between the levels of X1 are of interest, where ~X1*X2 indicates that a difference between the level combinations of the the two factors are of interest. Factors included in the formula must be included in \var{fixed}.}
  \item{contrast}{A character string indicating which specific levels of \var{contrast.factor} the user wishes to compare. This is usually of the form "TrtA - TrtB".  For contrast involving level combinations of the form ~X1*X2, the levels are concatonated ie "TrtADay1 - TrtADay0".  The order of concatenation must conform to the order in \var{contrast.factor}. }
  \item{random}{An optional formula to specify random effects and is passed directly to the \link[nlme]{lme} function. For simple repeated measures study designs, the form is usually ~ 1|g where g specifies the donor or subject id's. Defaults to NULL. }
  \item{correlation}{An optional \link[nlme]{corStruct} object to specify within-group correlation structure through the residual covariance matrix and is passed directly to the \link[nlme]{lme} or \link[nlme]{gls} functions. Defaults to NULL.}
  \item{design.sampleid}{A character string indicating a variable name in \var{design} indicating the column names of \var{eset}. If set to NULL, it is assumed that the column ordering in \var{eset} corresponds to the row ordering in \var{design}. Defaults to NULL.}
  
  
}
\details{

This function runs the entire qusage method on the input data, returning a single QSarray object containing the results. Rather than conducting gene level Welch's or paired t-tests, the user can specify more general linear models inherent to \link[nlme]{lme} and \link[nlme]{gls} functions within the \link{nlme} package. This requires the specification of a design file to link the necessary covariates and repeated measures information. One consideration when running more complex linear models is that complexity combined with poorly behaved data can sometimes yield convergence issues during the optimization step.  We encourage users to run their models on individual genes through gls and lme directly to ensure their models are paramaterized correctly and as expected.  Additional filtering of probes may need to be conducted for genes that have little variation as this can cause convergence issues.

Gene sets are commonly obtained from online databases such as Broad's Molecular Signatures Database. Gene set lists can be obtained from these sites in the form of .gmt files, which can be read into R using the \code{\link{read.gmt}} function. Once the data has been read into R, the information can be passed into the qusage function as either a vector describing a single gene set, or a list of vectors representing a group of gene sets. Each pathway must be a character vector with entries matching the row names of \var{eset}. If a pathway does not contain any values matching the rownames of \var{eset}, a warning will be printed, and the function will return NAs for the values of that pathway.


}
\value{
 A \link[=QSarray-class]{QSarray} object.
}

\examples{
##Creating a design file of 20 patients (10 in conditionA/10 in conditionB 
# with 5 timepoints)
# In addition, we also create a dummy continous covariate that is partially confounded
# with condition. Condition B will have a higher mean for the covariate than condition A.

des<-data.frame(SampleID=paste("Sample",1:100,sep=""),
                Condition=rep(c("A","B"),each=50),
                Donor=rep(letters[1:20],each=5),
                Time=rep(paste("D",0:4,sep=""),20))


##Create example data - a set of 500 genes normally dstributed across 20 patients 
##with 5 timepoints
eset = matrix(rnorm(500*100),500,100, dimnames=list(1:500,1:100))
colnames(eset)= paste("Sample",1:100,sep="")

##create a number of gene sets with varying levels of differential expression between 
##conditions and between timepoints
geneSets = list()
for(i in 0:10){
  genes = ((30*i)+1):(30*(i+1))
  eset[genes,des$Condition=="B"] = eset[genes,des$Condition=="B"] + rnorm(1)
  eset[genes,des$Time=="D1"]=eset[genes,des$Time=="D1"]+rnorm(1)
  
  geneSets[[paste("Set",i)]] = genes
}

##Adding additional subject specific variability to generate repeated measures correlations
for(i in 1:500){
  eset[i,]<-eset[i,]+rep(rnorm(20,0,4),each=5)
}



##Running a linear mixed model to test for D1 vs D0 for Condition B, with a subject
## specific random effect

qs.result<-qgen(eset,des,geneSets=geneSets,
                fixed= ~Condition+Time+Time*Condition,random=~1|Donor,
                contrast.factor=~Condition*Time,contrast="BD1-BD0",
                design.sampleid="SampleID")

plot(qs.result)
}
