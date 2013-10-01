## This file contains the primary code for the qusage methodology
##
## Author: Christopher Bolen
##         Gur Yaari
## Updated: 2013-03-15
## (c) 2013 Yale University. All rights reserved.


#require(limma)

##THE BIG WRAPPER FUNCTION
## This function runs the entire qusage method on the input data, returning a single
## QSarray object containing the results of makeComparison, calcVIF, and aggregateGeneSet.
## Many of the parameters are left out of this function for simplicity, so for greater control
## each of the functions must be called separately.

qusage = function(eset,             ##a matrix of log2(expression values), with rows of features and columns of samples. OR an object of class ExpressionSet 
                  labels,           ##vector of labels representing each column of eset.
                  contrast,         ##a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'. 
                  geneSets,         ##a list of pathways to be compared. Each item in the list is a vector of names that correspond to the row names of eset.
                  pairVector=NULL,  ##A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
                  var.equal=FALSE,  ##a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.
                  filter.genes=FALSE ##a boolean indicating whether the genes in eset should be filtered to remove genes with low mean and sd.
                 ){
  cat("Calculating gene-by-gene comparisons...")
  results = makeComparison(eset, labels, contrast, pairVector=pairVector,var.equal=var.equal)
  if(filter.genes){
    results = filterGenes(results)
  }
  
  cat("Done.\nAggregating gene data for gene sets.")
  results = aggregateGeneSet(results, geneSets, silent=F)
  cat("Done.\nCalculating variance inflation factors...")
  results = calcVIF(eset, results)
  #cat("Done.\nCalculating homogeneity scores...")
  #results = calcHomogeneity(results, silent=TRUE, addVIF=FALSE)
  #cat("Done.\nCalculating correlation matrix...")  
  #results = calcPCor(eset, results)
  cat("Done.\n")
  results
}

## A wrapper function to calculate comparisons between groups in a dataset.
## 
## A note on var.equal: LIMMA's linear model function can only be run when assuming equal 
## variances. If var.equal==TRUE, then a linear model will be created on the entire dataset 
## at once. Otherwise, calcIndividualExpressions will be called for the comparison of interest,
## and means and SDs will be calculated directly.
## One benefit of using LIMMA's pooled variance calculation is that the linear models allow for more
## complicated comparisons (e.g. "(A+B)-C" or similar). This may be of interest to some users,
## but in order to do this, you must assume equal variances between all groups.
##
## One caveat regarding paired samples: LIMMA can not fit a linear model when the paired samples are
## convoluted with the groups (e.g. one set of paired (trt vs mock) samples in patients with  disease,
## combined with a set of paired samples from healthy controls). If var.equal==TRUE, these groups
## must be run separately to correctly fit the model (e.g. run disease first, then healthy controls).

makeComparison <- function(eset,       ##a matrix of log2(expression values), with rows of features and columns of samples 
                           labels,     ##vector of labels representing each column of eset.
                           contrast,   ##a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'. 
                           pairVector=NULL,  ##A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
                           var.equal = FALSE, ##a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used. 
                           bayesEstimation = TRUE, ##if true, use a bayesian framework to estimate the standard deviation (via limma's eBayes function). 
                           min.variance.factor=10^-8  ##a factor to add to the SDs to ensure that none are equal to 0. Only used if var.equal==FALSE or bayesEstimation==FALSE.
){
  
  if(is(eset, "ExpressionSet")){eset = exprs(eset)}
  ##check that input is formatted correctly
  if(length(labels)!=ncol(eset)){stop("labels length does not match columns of eset")}
  labels = as.factor(as.vector(labels))
  if(!is.character(contrast)){stop("Contrast must be a character vector of length 1.")}
  if(length(contrast)!=1){
    warning("Multiple contrasts provided. Using first contrast only.")
    contrast = contrast[1]
  }
  if(is.null(rownames(eset))){
    stop("Rownames for eset not found")
  }
  if(length(unique(rownames(eset)))!=nrow(eset) | any(rownames(eset)=="")){
    stop("The rownames of eset are invalid. Rownames must be unique and must not contain any empty values")
  }
  
  
  params = list(labels=labels, contrast = contrast)
  
  if((paired = !is.null(pairVector))){
    if(length(pairVector)!=ncol(eset)){stop("PairVector length does not match columns of eset")}
    pairVector = as.factor(as.vector(pairVector))  
    params[["pairVector"]] = pairVector
  }
  
  
  if(var.equal){
    ###################################
    ## Pooled Variance (Linear Model) method
    
    ##create design matrix
    f = "~0+labels"
    designNames = levels(labels)
    if(paired){
      f = paste(f,"+pairVector",sep="")
      designNames = c(designNames, paste("P",levels(pairVector)[-1],sep=""))
    }
    design <- model.matrix(formula(f))
    colnames(design) <- designNames
      
    ## Fit the linear model with the given deign matrix
    ## 'fit' contains info on each coefficient (i.e. column of design matrix) in the model.
    fit <- lmFit(eset, design=design)
    
    ##Contrast the coefficients against each other to get a direct comparison.
    contrast.matrix <- makeContrasts( contrasts=contrast, levels=design)
    fit2 <- contrasts.fit(fit,contrast.matrix)
    
    ##if using Bayes estimation, calculate the moderated t-statistics for each comparison
    if(bayesEstimation){
      fit2b <- eBayes(fit2)
      SD = (fit2b$coefficients/(fit2b$t))[,1]
      sd.alpha = SD/(fit2b$sigma*fit2b$stdev.unscaled)
      sd.alpha[is.infinite(sd.alpha)] = 1
      dof = fit2b$df.total
    }else{
      SD = sqrt((fit2$sigma*fit2b$stdev.unscaled)^2 + min.variance.factor)
      sd.alpha = SD/(fit2$sigma*fit2$stdev.unscaled)
      sd.alpha[is.infinite(sd.alpha)] = 1
      dof = fit2$df.residual
    }
    
    ##format 
    results = newQSarray(params, 
                      mean = fit2$coefficients[,1],
                      SD = SD,
                      sd.alpha = sd.alpha,
                      dof = dof,
                      var.method="Pooled"
                     )
  }
  if(!var.equal){
    ###################################
    ## Welch's method
    ##parse the contrast
    grps = strsplit(contrast,"-")[[1]]
    grps = sub("\\s","",grps)          ##remove whitespace
    if(length(grps)!=2){stop("Only contrasts of the form 'A-B' are allowed when var.equal is FALSE.")}
    grp.1 = labels==grps[1]  ##PostTreatment
    grp.2 = labels==grps[2]  ##Baseline
    if(sum(grp.1)==0 | sum(grp.2)==0){stop("Contrast groups do not match labels")}
    
    eset.1 = eset[,grp.1]
    eset.2 = eset[,grp.2]
    
    if(paired){
      colnames(eset.1) = pairVector[grp.1]
      colnames(eset.2) = pairVector[grp.2]
      
      if(ncol(eset.1)!=ncol(eset.2)){
        eset.1 = eset.1[,colnames(eset.1) %in% colnames(eset.2)]
        eset.2 = eset.2[,colnames(eset.2) %in% colnames(eset.1)]
      }
    }
    
    results = newQSarray(c(params, calcIndividualExpressions(eset.2,eset.1,paired=paired,min.variance.factor=min.variance.factor)))
  }
  return(results)
}

#######Calculate individual gene differential expresseion for each gene
## This function should usually be called by "makeComparison".

calcIndividualExpressions<-function(Baseline,PostTreatment,paired=FALSE,min.variance.factor=10^-6,na.rm=TRUE){
  ###Baseline is the matix of gene expressions at baseline, row names are gene names
  ###PostTreatment is the matix of gene expressions after treatment, row names are gene names
  ###paired: logical, whether the data is paired or not
  
  ##########Some error checks
  if(length(dim(Baseline))!=2 | length(dim(PostTreatment))!=2){stop("Input Matrices need to be matrices... \n")}
  if(nrow(Baseline)!=nrow(PostTreatment)){stop("Input Matrices need to have the same number of genes \n")}
  if(sum(!(rownames(Baseline)%in%rownames(PostTreatment)))){stop("Input Matrices need to have the same list of genes. Gene names should be the row names \n")}
  if(ncol(Baseline)<2 | ncol(PostTreatment)<2 ){stop("Input Matrices need to have at least two columns \n")}

  #########Reorder PostTreatment
  PostTreatment<-PostTreatment[rownames(Baseline),]
  
#  if(!abs){
    ###########Paired
    if(paired){
      if(ncol(Baseline)!=ncol(PostTreatment)){
        stop("Input Matrices need to have the same number of columns when paired flag is on \n")
      }
      if(sum(!colnames(Baseline)%in%colnames(PostTreatment))){
        stop("Input Matrices need to have the same list of samples when paired flag is on \n")
      }
      PostTreatment = PostTreatment[,colnames(Baseline)]
      ##########First calculate the differential expression for individual genes
      Sums_Base<-rowSums(Baseline, na.rm=na.rm)
      Sums_Post<-rowSums(PostTreatment, na.rm=na.rm)
      Ns<-rowSums(!is.na(Baseline-PostTreatment), na.rm=na.rm)
      if(min(Ns)!=ncol(Baseline)){warning("Some NA's in data")}
      Sigmas_Base<-rowSums((Baseline-PostTreatment-(Sums_Base-Sums_Post)/Ns)^2)/(Ns-1)
      DOF<-Ns
      if(any(DOF<3, na.rm=T)){warning("Some degrees of freedom are below minimum. They have been set to 3.")}
      DOF[DOF<3]<-3
      Mean=(Sums_Post-Sums_Base)/Ns
      SD=sqrt(Sigmas_Base/Ns)
    }
    ###########Non Paired
    if(!paired){
      ##########First calculate the differential expression for individual genes
      Sums_Base<-rowSums(Baseline, na.rm=na.rm)
      Sums_Post<-rowSums(PostTreatment, na.rm=na.rm)
      Ns_Base<-rowSums(!is.na(Baseline), na.rm=na.rm)
      Ns_Post<-rowSums(!is.na(PostTreatment), na.rm=na.rm)
      if(min(Ns_Base)!=ncol(Baseline) | min(Ns_Post)!=ncol(PostTreatment)){warning("Some NA's in data")}
      Sigmas_Base<-rowSums((Baseline-(Sums_Base)/Ns_Base)^2)/(Ns_Base-1)
      Sigmas_Post<-rowSums((PostTreatment-(Sums_Post)/Ns_Post)^2)/(Ns_Post-1)
      ROWS<-rownames(Baseline)
      DOF<-Ni(Sigmas_Post+min.variance.factor,Sigmas_Base+min.variance.factor,Ns_Post,Ns_Base)
      #calculate degrees of freedom
      if(any(DOF<3, na.rm=T)){warning("Some degrees of freedom are below minimum. They have been set to 3.")}
      DOF[DOF<3]<-3
      Mean=(Sums_Post/Ns_Post-Sums_Base/Ns_Base)
      SD=sqrt(Sigmas_Base/Ns_Base+Sigmas_Post/Ns_Post)
    }
  sd.alpha = sqrt(SD^2+min.variance.factor)/SD
  sd.alpha[is.infinite(sd.alpha)] = 1
  
  dat = newQSarray(mean=Mean,
                SD=sqrt(SD^2+min.variance.factor),
                sd.alpha = sd.alpha,
                dof=DOF,
                var.method="Welch's"
  )
  dat
}
 
 
 #####Filter out genes with SD lower than a threshold and mean lower than another threshold
 
filterGenes<-function(QSarray,  ##A QSarray object, as generated by makeComparison 
                      Min_SD=0.01, ##threshold  for SD
                      Min_Mean=0   ##threshold  for mean
                     ){
 
 if(!is.null(QSarray$pathways)){stop("too late...aggregateGeneSet already being called")}
 
 Indexes <- abs(QSarray$mean) >= Min_Mean | QSarray$SD >= Min_SD
 
 QSarray$mean <- QSarray$mean[Indexes]
 QSarray$SD <- QSarray$SD[Indexes]
 QSarray$dof <- QSarray$dof[Indexes]
 QSarray$sd.alpha <- QSarray$sd.alpha[Indexes]
 
 return(QSarray)
}
 
##Simple function to read in a .gmt file and return a list of pathways
read.gmt = function(file){
  if(!grepl("\\.gmt$",file)[1]){stop("Pathway information must be a .gmt file")}
  geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
  geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
  names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
  geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
  geneSetDB = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
  return(geneSetDB)
}
 
#######Combine individual gene differential expresseion for each pathway (Neg) ~ 1 minute

aggregateGeneSet<-function(geneResults,  ##A QSarray object, as generated by makeComparison 
                           geneSets,     ##a list of pathways to be compared, each item in the list is a vector of names that correspond to the gene names from Baseline/PostTreatment
                           n.points=2^12, ##the number of points to sample the convoluted t-distribution at.
                           silent=TRUE   ##If false, print a "." every fifth pathway, as a way to keep track of progress
                          ){
  
#   NumSDs<-c(20,20,20,10,5,2,2,2,2,1,1,rep(.5,220))*30
  NumSDs<-c(1000,abs(qt(10^-12,2:250)))
#   ,rep(abs(qt(10^-8,50)),220))  
  
  Means = geneResults$mean
  SD = geneResults$SD
  DOF=geneResults$dof
  COLS = names(Means)
  
  if(is.vector(geneSets) & !is.list(geneSets)){
    n = deparse(substitute(geneSets))
    geneSets = list(geneSets)
    names(geneSets) = n
  }
  if(is.null(names(geneSets))){names(geneSets) = 1:length(geneSets)}
  geneSets = lapply(geneSets,function(x){
    if(is.numeric(x)){
      if(any(!(x %in% 1:length(COLS)))){stop("Numeric gene set indices out of bounds")}
      return(x)
    }
    which(COLS%in%x)
  })
  
  #########First set MaxDiff to adjust to data:
  ##calculate standard deviation
  SumSigma<-sapply(names(geneSets),function(i){
      Indexes = geneSets[[i]]
      x<-sqrt(sum((SD^2*(DOF/(DOF-2)))[Indexes]))
      return(x)
  })
  
  MinDof<-sapply(names(geneSets),function(i){
      Indexes = geneSets[[i]]
      if(length(Indexes)==0){return(NA)}
      return(floor(min(DOF[Indexes])))
      })
#   MaxDiff<-pmax(NumSDs[floor(min(DOF))]*SumSigma,1)  
  MaxDiff<-pmax(NumSDs[MinDof]*SumSigma,1,na.rm=TRUE)
  PDFs = pathMeans = Sizes = NULL
  for(i in 1:length(geneSets)){
    if(!silent & i%%5==0){cat(".")}
    Indexes = geneSets[[i]]
    if(length(Indexes)!=0){
      Norm<-(2*MaxDiff[i]/{n.points-1}) #normalize 
      PDF<-sapply(Indexes,function(j){
          x = SD[j]
          MaxDiffInt<-MaxDiff[i]/x
          #y<-dt(seq(-MaxDiffInt,MaxDiffInt,length.out=n.points),DOF)  
          x1<-seq(-MaxDiffInt,MaxDiffInt,length.out=n.points)
          y<-dt(x1[1:(n.points/2)],DOF[j])
          y<-c(y,rev(y))
          y/sum(y)/Norm
      })
      Tmp<-multi_conv(PDF)
      Tmp = Tmp*(n.points-1)/MaxDiff[i]/2
    }else{
      warning(paste("Gene set: (index ",i, ") has 0 overlap with eset.",sep=""))
      Tmp = rep(NA, n.points)
    }
    PDFs = cbind(PDFs,Tmp)
    pathMeans = c(pathMeans, mean(Means[Indexes]))
    Sizes = c(Sizes, length(Indexes))
  }
  
  colnames(PDFs) = names(pathMeans) = names(Sizes) = names(geneSets)
  ##add the new data to the existing QSarray object
  geneResults$pathways = geneSets
  geneResults$path.mean = pathMeans
  geneResults$path.size = Sizes
  geneResults$ranges = MaxDiff/Sizes
  geneResults$n.points=n.points
  geneResults$path.PDF = PDFs
  return(geneResults)
}

###A function to calculate the Variance Inflation Factor (VIF) for each of the gene sets input in geneSets
## This function depends on which method you used to calculate the variance originally. 
## If you assumed pooled variance, then the variance will be calculated using LIMMA's
## "interGeneCorrelation" method (see Wu and Smyth, Nucleic Acids Res. 2012). Otherwise, the method
## will calculate the VIF separately for each group and return the average of each group's vif.

calcVIF = function(eset,         ##a matrix of log2(expression values). This must be the same dataset that was used to create geneResults
                   geneResults,  ##A QSarray object, as generated by either makeComparison or aggregateGeneSet
                   useCAMERA = geneResults$var.method=="Pooled",  ##The method used to calculate variance. See the description for more details. By default, it uses the parameter stored in the geneResults
#                    geneSets=NULL, ##a list of pathways calculate the vif for, each item in the list is a vector of names that correspond to the gene names from Baseline/PostTreatment
                   useAllData = TRUE ##Boolean parameter determining whether to use all data in eset to calculated the VIF, or whether to only use data from the groups being contrasted. Only used if useCAMERA==FALSE
                  ){
  
#   if(is.null(geneSets)){
    if(is.null(geneResults$pathways)){stop("Pathway Information not found. Please provide a list of gene sets.")}
    geneSets = geneResults$pathways
#   }else{
#     if(is.vector(geneSets) & !is.list(geneSets)){
#       n = deparse(substitute(geneSets))
#       geneSets = list(n = geneSets)
#       names(geneSets) = n
#     } 
#     geneSets = lapply(geneSets,function(x){which(names(geneResults$mean)%in%x)})
#     geneResults$pathways = geneSets
#   }
    
  if(is(eset, "ExpressionSet")){eset = exprs(eset)}
  if(class(geneResults) != "QSarray"){stop("geneResults must be a QSarray object, as created by makeComparison")}
    
  ##create design matrix
  if(useCAMERA){
    labels = geneResults$labels
    paired=F
    if(!is.null(geneResults$pairVector)){paired=T; pairVector = geneResults$pairVector}
    
    f = "~0+labels"
    designNames = levels(labels)
    if(paired){
      f = paste(f,"+pairVector",sep="")
      designNames = c(designNames, paste("P",levels(pairVector)[-1],sep=""))
    }
    design <- model.matrix(formula(f))
    colnames(design) <- designNames
  }
  
  ##run VIF calculation on each gene set
  vif = sapply(names(geneSets),function(i){
    GNames<-names(geneResults$mean)[geneSets[[i]]]
    gs.i = which(rownames(eset)%in%GNames)
#     gs.i = geneSets[[i]]
    if(length(gs.i)<2){warning("GeneSet '",i,"' contains one or zero overlapping genes. NAs produced.");return(NA)}    
    if(useCAMERA){
      return(interGeneCorrelation(eset[gs.i,],design)$vif)
    }
    else{
#       grps = sub("\\s","",strsplit(geneResults$contrast,"-")[[1]])  ##only calc vif for the groups that were compared
#       vif.grp = sapply(split(1:ncol(eset),geneResults$labels)[grps], function(j){
#         covar.mat = cov(t(eset[gs.i,j]))
#         return(sum(covar.mat)/sum(diag(covar.mat)))
#       })
#       return(mean(vif.grp,na.rm=T))
      
      ##pooled covariance matrix
      grps = split(1:ncol(eset),geneResults$labels)
      if(!useAllData){
        toInclude = sub("\\s","",strsplit(geneResults$contrast,"-")[[1]])  ##only calc vif for the groups that were compared
        grps = grps[toInclude]
      }
      covar.mat = cov(t(eset[gs.i,grps[[1]]])) * (length(grps[[1]])-1)
      if(length(grps)>1){
        for(i in 2:length(grps)){
          covar.mat = covar.mat + ( cov(t(eset[gs.i,grps[[i]]])) * (length(grps[[i]])-1) )
        } 
      }
      covar.mat = covar.mat / (ncol(eset) - length(grps))
      
      ##multiply matrix by the sd.alpha vectors
      if(!is.null(geneResults$sd.alpha)){
        a = geneResults$sd.alpha[rownames(eset)[gs.i]]
        covar.mat = t(covar.mat*a)*a
      }
      
      vif = sum(covar.mat)/sum(diag(covar.mat))
      return(vif)
    }
  })
  geneResults$vif = vif
  if(!is.null(geneResults$path.PDF)){ ##if defined, rescale the pdf with the new vif values
    geneResults$path.PDF = t(t(geneResults$path.PDF) / pdfScaleFactor(geneResults))
  }
  return(geneResults)
}


##!!----deprecated. Use pdf.pVal instead
oneWay.pVal <- function(QSarray, alternative=c("two.sided","less","greater"),
                        direction=FALSE, addVIF=!is.null(QSarray$vif), selfContained=TRUE){#, absolute=FALSE){
  warning("Hey -- I changed the name. 'pdf.pVal' is less confusing. Use that.")
  return(pdf.pVal(QSarray, alternative, direction,addVIF,selfContained))#,absolute))
}
## function for calculating a p-value for each pathway convolution as output by aggregateGeneSet.
pdf.pVal <- function(QSarray,                          ##The output of qusage (or aggregateGeneSet)
                     alternative=c("two.sided","less","greater"), ## If alternative = "two.sided", this tests whether the mean of the distribution is equal to 0. If alternative = "greater", the alternative is that the mean is greater than than 0.
                     direction=FALSE,                  ##If true (and alternative="two.sided"), p-values will be returned as eiter positive or negative if the pathway is greater or less than 0, respectively.
                     addVIF=!is.null(QSarray$vif),     ##Whether to include the VIF in the calculation of the p-values. By default, if the VIF has been calculated, it will be used
                     selfContained=TRUE                ##If false, rather than comparing to 0, it will compare the pathway mean to the mean of all genes not in the pathway.
                     #absolute=FALSE                    ##If true, computes a P value for the "absolute test".
){ 
  if(is.null(QSarray$path.PDF)){stop("convolution results not found.")}
  p = sapply(1:ncol(QSarray$path.PDF), function(i){
    if(QSarray$path.size[i]==0){return(NA)}
    
    #calculate null hypothesis
    if(!selfContained){
      path = QSarray$pathways[[i]]
      null.hyp = mean(QSarray$mean[-path])
    }else{
      null.hyp=0
#       if(!absolute){null.hyp=0}
#       else{
#         path = QSarray$pathways[[i]]
#         null.hyp=mean(getExAbs(QSarray$dof[path])*QSarray$SD[path])
#       }    
    }
    ##find location of null hypothesis and sum pdf
    x = getXcoords(QSarray,i,addVIF=addVIF)#,absolute=absolute)
    PDF_NORM<-QSarray$path.PDF[,i]/sum(QSarray$path.PDF[,i])
#     sum(QSarray$path.PDF[1:findInterval(0,x),i]) / sum(QSarray$path.PDF[,i])
    INDEX<-findInterval(null.hyp,x)
    if(INDEX==0){return(0)}; if(INDEX==length(x)){return(1)}
    sum(PDF_NORM[1:INDEX])+((null.hyp-x[INDEX])/(x[INDEX+1]-x[INDEX]))*PDF_NORM[INDEX+1]
  })
  
  dir = match.arg(alternative)
  if(dir=="two.sided"){
    p[which(p>0.5)] = -1+p[which(p>0.5)]   ### turn it into a two-tailed p
    p = ifelse(rep(direction,length(p)),p*2,abs(p*2))
    return(p)
  }
  if(dir=="less"){return(1-p)}
  return(p)
}


##!!----deprecated. Use twoCurve.pVal instead
twoWay.pVal <- function(grp1, grp2, path.index1 = 1:numPathways(grp1), path.index2 = 1:numPathways(grp2), 
                        alternative=c("two.sided","less","greater"), direction=FALSE, 
                        addVIF=!(is.null(grp1$vif) | is.null(grp2$vif))){
  warning("Hey -- I changed the name. 'twoCurve.pVal' is less confusing. Use that.")
  return(twoCurve.pVal(grp1, grp2, path.index1, path.index2, alternative, direction,addVIF))
}

## A method to compare the pathway convolutions in two QSarray objects. 
## If alternative = "two.sided", this tests whether the mean of the two distributions are equal.
## If alternative = "greater", the alternative is that grp1 is greater than group2
twoCurve.pVal<-function(grp1, grp2,
                        path.index1 = 1:numPathways(grp1),
                        path.index2 = 1:numPathways(grp2), 
                        alternative=c("two.sided","less","greater"),
                        direction=FALSE, 
                        addVIF=!(is.null(grp1$vif) | is.null(grp2$vif))
                     ){
  ##if the names of the pathways don't match, 
#   if(ncol(grp1$path.PDF)!=ncol(grp2$path.PDF) | all(colnames(grp1$path.PDF) != colnames(grp2$path.PDF))){
#     stop("Pathways in grp1 do not match pathways in grp2")
#   }
  if(length(path.index1)!=length(path.index2)){
    stop("Number of pathways in grp1 do not match number of pathways in grp2")
  }
#     if(sum(names(grp1$path.mean[path.index1])!=names(grp2$path.mean[path.index2]))){
#         warning("Some of the comparisons are made between different pathways")
#   }
  alternative=match.arg(alternative)
  x1 = sapply(path.index1,function(i){getXcoords(grp1,i,addVIF=addVIF)})
  x2 = sapply(path.index2,function(i){getXcoords(grp2,i,addVIF=addVIF)})
  Length1 = grp1$n.points
  Length2 = grp2$n.points
  Min<-apply(rbind(x1,x2),2,min)
  Max<-apply(rbind(x1,x2),2,max)
  p = sapply(1:length(path.index1), function(i){
    if(is.na(Max[i])){return(NA)}
    PDF1<-approx( x1[,i], grp1$path.PDF[,path.index1[i]],seq(Min[i],Max[i],length.out=Length1+Length2),rule=2)$y
    PDF2<-approx( x2[,i], grp2$path.PDF[,path.index2[i]],seq(Min[i],Max[i],length.out=Length1+Length2),rule=2)$y
    return(compareTwoDistsFaster(PDF1,PDF2, alternative=alternative))
  })
  ifelse(rep(direction,length(p)), p, abs(p))
}

## Calculates the x-coordinates for the PDF of a given pathway.
## path.index can either be an integer between 1 and length(path.means), or the name of the pathway.
## addVIF is a boolean determining whether information on the VIF should be used in the calculation
## of the x-coordinates.
getXcoords = function(QSarray,path.index=1, addVIF=!is.null(QSarray$vif)){ #,absolute=FALSE){
  if(length(path.index)>1){stop("path.index must be of length 1")}
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}
  sif = ifelse(addVIF,sqrt(QSarray$vif[path.index]),1)
  if(is.na(sif)){sif=1}
# if(!absolute){
   seq(-1,1,length.out=QSarray$n.points)* QSarray$ranges[path.index]* sif + QSarray$path.mean[path.index]
# }
#  else {
#  ###First calculate the new mean of the pathway based on the absolute values of the means
#  MeanAbs<-mean(abs(QSarray$mean[QSarray$pathways[[path.index]]]))
#  seq(-1,1,length.out=QSarray$n.points)* QSarray$ranges[path.index]* sif / QSarray$path.size[path.index] + MeanAbs
#  }
}


##calculate a scaling factor to multiply the PDF by (so that it actually integrates to 1)
##this is primarily used for plotting of the PDF.
pdfScaleFactor = function(QSarray, addVIF=!is.null(QSarray$vif)){
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}
  sif = sapply(QSarray$vif,function(v){ifelse(addVIF,sqrt(v),1)})
  sif[is.na(sif)] = 1
  pdfSum = colSums(QSarray$path.PDF)
  ##the scale factor is essentially the distance between points in the x coordinates times the sum of the pdf
  scaleFactor = 2 * (QSarray$ranges*sif) / (QSarray$n.points-1) * pdfSum
  scaleFactor
}


## Computes the 95% CI for a pdf
calcBayesCI <- function(QSarray,low=0.025,up=1-low,addVIF=!is.null(QSarray$vif)){
  cis = sapply(1:length(QSarray$path.mean), function(i){
    if(length(QSarray$pathways[[i]])==0){return(c(NA,NA))}
    x = getXcoords(QSarray,i,addVIF=addVIF)
    cdf = cumsum(QSarray$path.PDF[,i])
    cdf = cdf/cdf[length(cdf)]
    INDEX_LOW<-findInterval(low,cdf)
    INDEX_UP<-findInterval(up,cdf)
    return( c(  x[INDEX_LOW]+ ((low-cdf[INDEX_LOW])/(cdf[INDEX_LOW+1]-cdf[INDEX_LOW]))*(x[INDEX_LOW+1]-x[INDEX_LOW]) ,
                x[INDEX_UP] + ((up-cdf[INDEX_UP])/(cdf[INDEX_UP+1]-cdf[INDEX_UP]))*(x[INDEX_UP+1]-x[INDEX_UP])
            ) )
#     return( c(x[findInterval(low,cdf)-1] , x[findInterval(up,cdf)]) )
  })
  colnames(cis) = names(QSarray$path.mean)
  rownames(cis) = c("low","up")
  return(cis)
}

## Compute p-value of two distributions not centered at zero
compareTwoDistsFaster <-function(dens1=runif(256*8,0,1), dens2=runif(256*8,0,1),alternative="two.sided"){
  if(length(dens1)>1 & length(dens2)>1 ){
    dens1<-dens1/sum(dens1)
    dens2<-dens2/sum(dens2)
    cum2 <- cumsum(dens2)-dens2/2
    tmp<- sum(sapply(1:length(dens1),function(i)return(dens1[i]*cum2[i])))
    if(alternative=="two.sided"){
      if(tmp>0.5)tmp<-tmp-1
      return( tmp*2 )
    }
    if(alternative=="less"){return(1-tmp)}
    return(tmp)
  }
  else {
    return(NA)
  }
}

################################FFT method for convolution

multi_conv<-function(x){
  x_fft<-apply(x,2,function(x)fft(x, inverse = FALSE))
  M<-max(Mod(x_fft))
  x_fft<-x_fft/M
  Prod_fft<-apply(x_fft,1,prod)
  p1<-Re(fft(Prod_fft,inverse=TRUE))
  N<-nrow(x)
  #     if(ncol(x)%%2==0)p1<-c(p1[(N/2+1):N],p1[1:(N/2)])
  Mp1<-which.max(p1)
  Delta<-N/2-Mp1
  if(Delta>0){
    p1<-c(p1[(N-Delta):N],p1[1:(N-Delta-1)])
  }
  if(Delta<0){
    p1<-c(p1[(1-Delta):N],p1[1:(1-Delta-1)])
  }    
  p2 = abs(p1)
  return(p2/sum(p2))
}

###A helper function to calculate degrees of freedom
Ni<-function(SigmaOld,SigmaYoung,NOld,NYoung){
 # (SigmaYoung^2/NYoung+SigmaOld^2/NOld)^2/(((SigmaYoung^4/(NYoung^2*(NYoung-1))))+((SigmaOld^4/(NOld^2*(NOld-1)))))
  (SigmaYoung/NYoung+SigmaOld/NOld)^2/(((SigmaYoung^2/(NYoung^2*(NYoung-1))))+((SigmaOld^2/(NOld^2*(NOld-1)))))
}


#############################################################################################
######################### Absolute Test / Homogeneity Score functions #######################


absoluteTest = function(eset, QSarray, p.adjust.method="fdr", silent=F){
  cat("Calculating Individual Gene Pvals.")
  gene.pvals = absoluteTest.genePvals(QSarray, compareTo="z",silent=silent)
  cat("Done.\nCalculating VIF...")
  corMats = calcPCor(eset, QSarray)
  cat("Done.\n")
  
  ##aggregate p-values 
  PVALS<-(sapply(1:length(gene.pvals),function(i){
    Chi2<-(-2)*sum(log(abs(gene.pvals[[i]])))
    k = length(gene.pvals[[i]])
    corMat = corMats[[i]]
    COR = 4*k+sum(ifelse(corMat > 0,corMat*(3.25 + 0.75*corMat),corMat*(3.27 + 0.71*corMat) ))
    f = 8*(k)^2/COR
    c = COR/4/k       
    Chi2<-Chi2/c
    P<-1-pchisq(Chi2,f)
  }))
  
  return(newQSarray(QSarray, absolute.p = PVALS))
}

homogeneityScore = function(QSarray){
  cat("Calculating Individual Gene Pvals.")
  gene.pvals = absoluteTest.genePvals(QSarray, compareTo="p")
  cat("Done.\n")
  
  ##aggregate p-values 
  PVALS<-(sapply(1:length(gene.pvals),function(i){
    Chi2<-(-2)*sum(log(abs(gene.pvals[[i]])))
    k = length(gene.pvals[[i]])
    P<-1-pchisq(Chi2,2*k)
  }))
  
  homogeneity = 1/(1-log10(PVALS))
  
  return(newQSarray(QSarray, homogeneity = homogeneity))
}

###A function to calculate the P values for each gene as compared with either 
## 0, the gene set mean, or the full PDF of the gene set. 


absoluteTest.genePvals<-function(QSarray,  ##A QSarray object, as generated by aggregateGeneSet and possibly modified by calcVIF
                                 path.index=1:numPathways(QSarray), ##The pathways to calculate the pVals for.
                                 silent=TRUE,   ##If false, print a "." every fifth pathway, as a way to keep track of progress  
                                 FastApproximated=TRUE, ##fast mode which uses normal approximations for the PDF's
                                 addVIF=!is.null(QSarray$vif), ##Whether to include the VIF in the calculation of the p-values. By default, if the VIF has been calculated, it will be used
                                 NPoints=QSarray$n.points/8, ##,
                                 compareTo=c("zero","mean","pdf") ##the null hypothesis that each gene set is tested against.
                                ){
  ##check path.index
  if(is.character(path.index)){
    path.index = match(path.index, names(QSarray$pathways))
  }
  
  NumSDs<-c(1000,abs(qt(10^-8,2:250)))  
  if(is.null(QSarray$pathways)){stop("Pathway Information not found. Please run aggregateGeneSet first.")}
  geneSets = QSarray$pathways
  if(addVIF==TRUE){addVIF=!is.null(QSarray$vif)}
  Means = QSarray$mean
  SD = QSarray$SD
  DOF=QSarray$dof
  Ps = list()
  if (FastApproximated) {
    SDPath = apply(calcBayesCI(QSarray,low=0.5,up=0.8413448)[,path.index,drop=F],2,function(x)x[2]-x[1])
    if(!addVIF)SDPath = SDPath / sqrt(QSarray$vif[path.index])
  }
  
  compareTo = match.arg(compareTo)
  if(compareTo=="pdf"){
    #NPoints=QSarray$n.points/8
    for(i in path.index){
      XPath = getXcoords(QSarray,i,addVIF=addVIF)
      if(!silent & i%%5==0){cat(".")}
      Indexes = geneSets[[i]]
      if(!FastApproximated){
        XGene <- seq(-1,1,length.out=NPoints) * NumSDs[floor(min(DOF[Indexes]))] 
        Min<-min(c(XGene[1]+ Means[Indexes],XPath[1]))
        Max<-max(c(XGene[NPoints]+ Means[Indexes],XPath[QSarray$n.points]))
        X_Sample<-seq(Min,Max,length.out=NPoints)
        PDFPath<-approx( XPath, QSarray$path.PDF[,i],X_Sample,rule=2)$y
        if(length(Indexes)!=0){
          PS<-NULL
          for(j in Indexes){          
            y<-dt(XGene[1:(NPoints/2)]/ SD[j],DOF[j])
            PDFGene<-c(y,rev(y))
            PDFGene<-approx( XGene+ Means[j], PDFGene,X_Sample,rule=2)$y
            PS<-c(PS,compareTwoDistsFaster(PDFGene,PDFPath, alternative="two.sided"))                    
          }
          Ps<-c(Ps,list(PS))
        }
      }
      else{
        if(length(Indexes)!=0){
          PS<-NULL
          for(j in Indexes){  
            TMP<-pnorm( ( Means[j] - QSarray$path.mean[i]  ) / sqrt( SDPath[i]^2 + (DOF[j])/(DOF[j]-2)*SD[j]^2) )
            if(TMP>0.5) TMP <- TMP - 1
            TMP <- TMP*2
            PS<-c(PS,TMP)
          }
          Ps<-c(Ps,list(PS))
        }
      }
    }
  }
  else{
    for(i in path.index){
      if(!silent & i%%5==0){cat(".")}
      Indexes = geneSets[[i]]
      if(length(Indexes)!=0){
        PS<-NULL
        for(j in Indexes){
          SUBSTRACT=0
          if(compareTo=="mean")SUBSTRACT=QSarray$path.mean[i]
          p<-pt((Means[j]-SUBSTRACT)/SD[j],DOF[j])
          p[which(p>0.5)] = -1+p[which(p>0.5)]   ### turn it into a two-tailed p
          p = p*2
          PS<-c(PS,p)                    
        }
        Ps<-c(Ps,list(PS))
      }
    }
  }
  names(Ps) = names(geneSets)[i]
  return(Ps)
}

###A function to calculate the "homogeneity" P value FAST


absoluteTest.genePvalsFAST<-function(QSarray,  ##A QSarray object, as generated by aggregateGeneSet and possibly modified by calcVIF
                              CompareWithZero=TRUE ###Logical, if TRUE compares with mean of zero, else with mean of pathway
){
  if(is.null(QSarray$pathways)){stop("Pathway Information not found. Please provide a list of gene sets.")}
  geneSets = QSarray$pathways
  Means = QSarray$mean
  SD = QSarray$SD
  DOF=QSarray$dof
  Ps = list()
  for(i in 1:length(geneSets)){
    Indexes = geneSets[[i]]
    if(length(Indexes)!=0){
      PS<-NULL
      
      PS<-sapply(Indexes,function(j){
        SUBSTRACT=0
        if(!CompareWithZero)SUBSTRACT=QSarray$path.mean[i]
        p<-pt((Means[j]-SUBSTRACT)/SD[j],DOF[j])
        p[which(p>0.5)] = -1+p[which(p>0.5)]   ### turn it into a two-tailed p
        p = p*2
        return(p)
      })
      Ps[[i]]<-PS
    }
  }
  names(Ps) = names(geneSets)
  return(Ps)
}


###A function to calculate the Variance Inflation Factor (VIF) for each of the gene sets input in geneSets
## This function depends on which method you used to calculate the variance originally. 
## If you assumed pooled variance, then the variance will be calculated using LIMMA's
## "interGeneCorrelation" method (see Wu and Smyth, Nucleic Acids Res. 2012). Otherwise, the method
## will calculate the VIF separately for each group and return the average of each group's vif.

calcPCor = function(eset,         ##a matrix of log2(expression values). This must be the same dataset that was used to create geneResults
                    geneResults,  ##A QSarray object, as generated by either makeComparison or aggregateGeneSet
                    useCAMERA = geneResults$var.method=="Pooled",  ##The method used to calculate variance. See the description for more details. By default, it uses the parameter stored in the geneResults
#                    geneSets=NULL, ##a list of pathways calculate the vif for, each item in the list is a vector of names that correspond to the gene names from Baseline/PostTreatment
                    useAllData = TRUE ##Boolean parameter determining whether to use all data in eset to calculated the VIF, or whether to only use data from the groups being contrasted. Only used if useCAMERA==FALSE
){
  if(is.null(geneResults$pathways)){stop("Pathway Information not found. Please provide a list of gene sets.")}
  geneSets = geneResults$pathways  
  if(class(geneResults) != "QSarray"){stop("geneResults must be a QSarray object, as created by makeComparison")}
  
  ##create design matrix
  if(useCAMERA){
    labels = geneResults$labels
    paired=F
    if(!is.null(geneResults$pairVector)){paired=T; pairVector = geneResults$pairVector}
    
    f = "~0+labels"
    designNames = levels(labels)
    if(paired){
      f = paste(f,"+pairVector",sep="")
      designNames = c(designNames, paste("P",levels(pairVector)[-1],sep=""))
    }
    design <- model.matrix(formula(f))
    colnames(design) <- designNames
  }
  
  Cor = sapply(names(geneSets),function(i){
    GNames<-names(geneResults$mean)[geneSets[[i]]]
    gs.i = which(rownames(eset)%in%GNames)
    #     gs.i = geneSets[[i]]
    if(length(gs.i)<2){warning("GeneSet '",i,"' contains one or zero overlapping genes. NAs produced.");return(NA)}    
    if(useCAMERA){
      return(interGeneCorrelation(eset[gs.i,],design)$vif)
    }
    else{
      ##pooled covariance matrix
      grps = split(1:ncol(eset),geneResults$labels)
      if(!useAllData){
        toInclude = sub("\\s","",strsplit(geneResults$contrast,"-")[[1]])  ##only calc vif for the groups that were compared
        grps = grps[toInclude]
      }
      covar.mat = cov(t(eset[gs.i,grps[[1]]])) * (length(grps[[1]])-1)
      if(length(grps)>1){
        for(i in 2:length(grps)){
          covar.mat = covar.mat + ( cov(t(eset[gs.i,grps[[i]]])) * (length(grps[[i]])-1) )
        } 
      }
      covar.mat = covar.mat / (ncol(eset) - length(grps))
      
      ##multiply matrix by the sd.alpha vectors
      if(!is.null(geneResults$sd.alpha)){
        a = geneResults$sd.alpha[rownames(eset)[gs.i]]
        covar.mat = t(covar.mat*a)*a
      }
      Cor = (covar.mat)
      Cor = sapply(1:ncol(Cor),function(i)Cor[i,]/sqrt(covar.mat[i,i]))
      Cor = sapply(1:ncol(Cor),function(i)Cor[,i]/sqrt(covar.mat[i,i]))
      Cor[!is.finite(Cor)] = 0
      return(Cor)
    }
  })
  
  return(Cor)
}



approximatedNu<-c(7.2144478,4.7860183,3.4886415,2.7484892,2.2968314,2.0029888,1.8005629,1.6541477,1.5438879,1.4580784,1.4085490,1.3518729,1.3046939,1.2648251,1.2306984,1.2011619,1.1753513,1.1526056,1.1324114,1.1143633,
                  1.1024454,1.0877347,1.0743740,1.0621861,1.0510235,1.0407625,1.0312984,1.0225419,1.0144168,1.0068572,0.9999555,0.9933631,0.9871863,0.9813873,0.9759323,0.9707917,0.9659392,0.9613512,0.9570067,0.9528868,
                  0.9489987,0.9452788,0.9417375,0.9383622,0.9351415,0.9320650,0.9291232,0.9263076,0.9236102,0.9210236,0.9185485,0.9161642,0.9138723,0.9116673,0.9095447,0.9074997,0.9055282,0.9036264,0.9017905,0.9000173,
                  0.8983066,0.8966494,0.8950460,0.8934938,0.8919904,0.8905335,0.8891210,0.8877509,0.8864214,0.8851306,0.8838784,0.8826602,0.8814760,0.8803245,0.8792042,0.8781140,0.8770525,0.8760188,0.8750118,0.8740303,
                  0.8730745,0.8721414,0.8712313,0.8703431,0.8694763,0.8686299,0.8678033,0.8669957,0.8662066,0.8654354,0.8646819,0.8579745,0.8524911,0.8479247,0.8440630,0.8407547,0.8378886,0.8353818,0.8331706,0.8312056,
                  0.8294480,0.8278664,0.8264358,0.8251355,0.8239485,0.8228606,0.8218598,0.8209362,0.8200810,0.8192871,0.8185480,0.8178582,0.8172130,0.8166081,0.8160400,0.8155053,0.8150013,0.8145252,0.8140749,0.8136484,
                  0.8132437,0.8128593,0.8124937,0.8121455,0.8118135,0.8114966,0.8111938,0.8109042,0.8106269,0.8103612,0.8101063,0.8098617,0.8096266,0.8094006,0.8091831,0.8089737,0.8087719,0.8085773,0.8083896,0.8082083,
                  0.8080332,0.8078639,0.8077002,0.8075418,0.8073884,0.8072398,0.8070957,0.8069561,0.8068206,0.8066891,0.8065614,0.8064373,0.8063168,0.8061996,0.8060856,0.8059747,0.8058667,0.8057616,0.8056593,0.8055595,
                  0.8054623,0.8053675,0.8052751,0.8051849,0.8050969,0.8050110,0.8049271,0.8048452,0.8047651,0.8046869,0.8046104,0.8045356,0.8044625,0.8043909,0.8043209,0.8042524,0.8041854,0.8041197,0.8040554,0.8039924,
                  0.8039306,0.8033756,0.8029140,0.8025239,0.8021900,0.8019008,0.8016481,0.8014253,0.8012274,0.8010504,0.8008913,0.8007473,0.8006165,0.8004972,0.8003878,0.8002872,0.8001944,0.8001085,0.8000287,0.7999545,
                  0.7998852,0.7998204,0.7997597,0.7997027,0.7996490,0.7995984,0.7995507,0.7995055,0.7994627,0.7994221,0.7993835,0.7993468,0.7993119,0.7992786,0.7992468,0.7992165,0.7991874,0.7991596,0.7991330,0.7991074,
                  0.7990829,0.7990593,0.7990367,0.7990149,0.7989939,0.7989737,0.7989542,0.7989354,0.7989172,0.7988996,0.7988827,0.7988663,0.7988504,0.7988350,0.7988201,0.7988057,0.7987917,0.7987781,0.7987650,0.7987522,
                  0.7987397,0.7987277,0.7987159,0.7987045,0.7986934,0.7986826,0.7986720,0.7986618,0.7986518,0.7986421,0.7986326,0.7986233,0.7986143,0.7986055,0.7985969,0.7985885,0.7985803,0.7985722,0.7985644,0.7985568,
                  0.7985493,0.7985419,0.7985348,0.7985278,0.7985209,0.7985142,0.7985076,0.7985012,0.7984949,0.7984887,0.7984826)

getExAbs<-function(Nu){
  approx(c(seq(1,10,0.1),seq(11,100,1),seq(110,1000,10)),approximatedNu,Nu)$y
}