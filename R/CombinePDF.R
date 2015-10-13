##########################################################################################
## "CombinedPDF" is the meta-analysis of QuSAGE that combines QuSAGE results from multiple 
## data sets. In this analysis, the QuSAGE probability density functions (PDFs) computed 
## for each data set were combined into one PDF for each gene set using numerical 
## convolution.
##########################################################################################
##
## Author: Hailong Meng, Gur Yarri, Chris Bolen
## Updated: 2015-07-29
## ? 2015 Yale University. All rights reserved.
##########################################################################################

## Plot the QScomb object for a specified pathway
plotCombinedPDF = function(QScomb, ##a QScomb object
                           path.index=1,
                           zeroLine=TRUE,
                           comb.lwd=3, path.lwd=1,
                           comb.col=par("col"), path.col=NULL,
                           legend=FALSE, legend.labs=NULL,
                           add=FALSE,
                           xlim=NULL,ylim=NULL,
                           xlab=NULL,ylab=NULL,main = NULL, 
                           type="l",
                           ...){
  
  
  if(class(QScomb) != "QSarray"){
    stop("Input object must be of class 'QSarray'")
  }
  if(is.null(QScomb$QSlist)){
    stop("Input QSarray object does not appear to be the result of a call to combinePDFs")
  }
  
  ##check that path.index is of length 1 (since we can only plot one at a time)
  if(length(path.index)>1){ 
    stop("Cannot plot more than one pathway. path.index must be of length 1")
  }
  ##get pathway name to ensure consistency across individual QSarrays
  ##then convert path.index to an index
  if(!is.numeric(path.index)){
    path.name = path.index
    path.index = match(path.index, colnames(QScomb$path.PDF))
  }else{
    path.name = colnames(QScomb$path.PDF)[path.index]
  }

  ## ---------------------------------------------
  ## set defaults for undefined input variables
  
  ##calculate x and y limits
  if(is.null(xlim)){ 
    xlim = range(calcBayesCI(QScomb)[,path.index], na.rm=T)
    xlim = xlim + (xlim-mean(xlim))*2 ##double the range
    qs.xlims = lapply(QScomb$QSlist, function(QSarray){
      xlim = range(calcBayesCI(QSarray)[,path.index],na.rm=T) 
      xlim + (xlim-mean(xlim))*2
    })
    xlim = range(c(xlim, unlist(qs.xlims)))
  }
  
  comb.scaleFactor = pdfScaleFactor(QScomb)
  if(is.null(ylim)){ 
    ymax=max(QScomb$path.PDF[,path.index]/comb.scaleFactor[path.index],na.rm=T) 
    qs.ymax = lapply(QScomb$QSlist, function(QSarray){
      scaleFactor = pdfScaleFactor(QSarray)
      max(QSarray$path.PDF[,path.index]/scaleFactor[path.index],na.rm=T) 
    })
    ylim = c(0, max(c(ymax, unlist(qs.ymax))))
  }
  
  ##set default colors and line widths
  default.pal = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                  "#FFFF33","#A65628","#F781BF","#999999")
  if(is.null(comb.col)){comb.col = par("col")}
  
  if(is.null(path.col)){path.col = default.pal}
  if(length(path.col)!=length(QScomb$QSlist)){
    path.col = rep(path.col, length.out=length(QScomb$QSlist))
  }
  if(length(path.lwd)!=length(QScomb$QSlist)){
    path.lwd = rep(path.lwd, length.out=length(QScomb$QSlist))
  }

  ##default axes labels
  if(is.null(xlab)){xlab="Pathway Activity"}
  if(is.null(ylab)){ylab="Density"}
  if(is.null(main)){main = path.name} 

  #create plot region
  if(!add){
    plot(0,type='n', ylab=ylab, xlab =xlab, main=main, xlim=xlim, ylim=ylim, ...)
  }
  if(zeroLine){ abline(v=0, col = "gray", lty = "dotted",lwd=2) }
  
  ##add individual curves
  for(i in 1:length(QScomb$QSlist)){
    plotDensityCurves(QScomb$QSlist[[i]], path.index=path.name, add=T, col=path.col[i], 
                      lwd=path.lwd[i], zeroLine = F, type=type)
  }
  
  ##add main curve
  x = getXcoords(QScomb, path.index = path.index)
  y = QScomb$path.PDF[,path.index]/comb.scaleFactor[path.index]
  lines(x,y,col=comb.col, lwd=comb.lwd, type=type)

  ##add legend
  if((is.logical(legend) && legend) || is.character(legend)){
    legend_col = c(path.col, comb.col)
    if(is.null(legend.labs)){
      if(!is.null(names(QScomb$QSlist))){legend_text=c(names(QScomb$QSlist), "Meta Analysis")}
      else{legend_text = c(paste0("Dataset",1:length(QScomb$QSlist)), "Meta Analysis")}
    }else{
      legend_text = c(legend.labs, "Meta Analysis")
    }
    loc = ifelse(is.character(legend), legend, "topleft")
    legend(loc,legend=legend_text, fill=legend_col,bty="n")
  }
}

## A wrapper function to combine QuSAGE results from multiple data sets and returns a CombinedPDFResult object
combinePDFs <- function(QSarrayList,        ##list of qusage resutls from different data sets 
                        n.points = 2^14
  ){ 
  num.datasets = length(QSarrayList)
  ##check that input is formatted correctly
  # if there is no QuSAGE object in the input QuSAGE object list
  if(num.datasets==0 ){
    stop("There is no qusage result object in the input.")
  }
  # if there is only one  QuSAGE object in the input QuSAGE object list
  if(num.datasets==1 ){
    stop("There is only one qusage result object in the input.")
  }
  
  
  # if pathways of QuSAGE objects in the input QuSAGE object list do not match
  name.pathway = lapply(QSarrayList, function(x){names(x$pathways)})
  num.pathways = unique(sapply(name.pathway, length))
  if(length(num.pathways)==1){
    name.pathway = do.call(cbind,name.pathway)
  }
  if(class(name.pathway)== "list" || any(name.pathway!=name.pathway[,1])){
    stop("The pathways in input QuSAGE objects do not totally match!")
  }

  ##pull sample size from the QSarray objects 
  ## (sum in case we're combining something we already combined)
  n.samples = sapply(QSarrayList, function(q){sum(q$n.samples)})  
  if(all(n.samples==0)){ n.samples = rep(1, length(n.samples))}
  
  ##calculate all the x-coordinates
  ## extract x coords from data set j
  X.lists = lapply(1:num.pathways, function(i){
    lapply(QSarrayList, getXcoords, path.index=i)
  })
  Min = sapply(X.lists, function(x){min(unlist(x),na.rm=T)})
  Max = sapply(X.lists, function(x){max(unlist(x),na.rm=T)})
  
  ## calculate combined PDFs for each pathway
  combined.PDF = sapply(1:num.pathways, function(i){ 
    X.list = X.lists[[i]]
    X.combined=seq(Min[i],Max[i],length.out=n.points)
    
    ## transform PDF to new X coordinates
    transformed.PDF = lapply(1:num.datasets, function(j){
      if(length(QSarrayList[[j]]$pathways[[i]])==0){
        warning("Dataset (index ",j,") contains no overlapping genes with pathway (index ",i,"). Excluding dataset from combined PDF.")
        return(NULL)
      }
      approx(X.list[[j]],QSarrayList[[j]]$path.PDF[,i],X.combined,rule=2)$y
    })
    toDrop = sapply(transformed.PDF, is.null)
    
    ##check to make sure at least one dataset has a PDF
    if(all(toDrop)){
      warning("No datasets found containing genes for pathway (index ",i,"). NAs produced")
      return(rep(NA, n.points))
      
    ##else, combine the PDF
    }else{
      transformed.PDF = do.call(cbind, transformed.PDF[!toDrop])
      ## combine PDFs and use sample size as weight parameter
      PDF.Combined<- calculate_bayesGHelper( list(transformed.PDF,n.samples[!toDrop]), n.points=n.points)
      
      NORM=(Max[i]-Min[i])/(n.points-1)
      combined.PDF.tmp <-PDF.Combined/NORM
      
      ##save data
      return(combined.PDF.tmp)
    }
  })
  ##make sure PDF.Combined is an array
  if(length(dim(combined.PDF))==0) combined.PDF<-t(matrix(combined.PDF))
  
  ##set column names
  colnames(combined.PDF) = name.pathway[,1]
  
  ##calculate the PDF means
  MEAN.PDFs<- sapply(1:num.pathways,function(i){
    if(!any(is.na(combined.PDF[,i]))){
      xin = seq(Min[i],Max[i],length.out=n.points)
      return(combined.PDF[,i] %*% xin / (sum(combined.PDF[,i])))
    }else{
      return(NA)
    }
  })
  
  ##calculate ranges (by taking the smaller tail)
  RANGES<- sapply(1:num.pathways,function(i){
    if (is.na(Max[i])) { return(NA) }
    return(min( (Max[i]-MEAN.PDFs[i]), (MEAN.PDFs[i]-Min[i]) ))
  })
  
  ##resample the PDFs to the new ranges
  PDF.COMBINED.RESAMPLED<-sapply(1:num.pathways,function(i){
    if (is.na(Max[i])) { return(rep(NA,n.points*2)) }
    
    xin = seq(Min[i],Max[i],length.out=n.points)
    xout = seq(MEAN.PDFs[i]-RANGES[i],MEAN.PDFs[i]+RANGES[i],length.out=n.points)
    return(approx(xin,combined.PDF[,i],xout)$y)
  })
  colnames(PDF.COMBINED.RESAMPLED) = colnames(combined.PDF)
  
  ##build up combined PDF result object
  results = newQSarray(path.PDF = PDF.COMBINED.RESAMPLED,
                       path.mean = MEAN.PDFs,
                       ranges = RANGES,
                       n.points = n.points,
                       QSlist = QSarrayList,
                       n.samples = n.samples
                     )

  return(results)
}
 
# ## Calculate one side p value for one single PDF
# cal.PValofCombinePDFs <- function (PDF, X_combined) 
# {
# 
#     null.hyp = 0
#     x = X_combined
#     PDF_NORM <- PDF/sum(PDF)
#     INDEX <- findInterval(null.hyp, x)
#     if (INDEX == 0) {
#         return(0)
#     }
#     if (INDEX == length(x)) {
#         return(1)
#     }
#     p=sum(PDF_NORM[1:INDEX]) + ((null.hyp - x[INDEX])/(x[INDEX + 
#         1] - x[INDEX])) * PDF_NORM[INDEX + 1]
# 
#     dir = "two.sided"
#     direction = FALSE
#     if (dir == "two.sided") {
#         p[which(p > 0.5)] = -1 + p[which(p > 0.5)]
#         p = ifelse(rep(direction, length(p)), p * 2, abs(p * 2))
#         return(p)
#     } 
#     return(p)
# }

## function calculate_bayesGHelper
calculate_bayesGHelper <- function( listMatG,n.points=2^14 ){
    matG <- listMatG[[1]]
    groups <- listMatG[[2]]
    i = 1
    resConv <- matG[,i]
    denom <- groups[i]
    if(length(groups)>1){
      while( i<length(groups) ){
        i = i + 1
        resConv <- weighted_conv(resConv, matG[,i], w= ((groups[i])/denom) ,n.points=n.points)
        #cat({{groups[i]}/denom},"\n")
        denom <- denom + groups[i]
      }
    }
    return(resConv)
}

## function weighted_conv 
weighted_conv<-function(x,y,w=1,m=100,n.points=4001){
  lx<-length(x)
  ly<-length(y)
  if( lx<m | (lx*w) < m | ly<m | (ly*w)<m ){
    if(w<1){
      y1<-approx(1:ly,y,seq(1,ly,length.out=m))$y
      x1<-approx(1:lx,x,seq(1,lx,length.out=m/w))$y
      lx<-length(x1)
      ly<-length(y1)
      
    } else {
      y1<-approx(1:ly,y,seq(1,ly,length.out=m*w))$y
      x1<-approx(1:lx,x,seq(1,lx,length.out=m))$y
      lx<-length(x1)
      ly<-length(y1)
    }
    
  } else{
    x1<-x
    y1<-approx(1:ly,y,seq(1,ly,length.out=floor(lx*w)))$y
    ly<-length(y1)
  }
  tmp<-approx(x=1:(lx+ly-1),y=convolve(x1,rev(y1),type="open"),xout=seq(1,lx+ly-1,length.out=n.points))$y
  tmp[tmp<=0] = 0 
  return(tmp/sum(tmp))
}

