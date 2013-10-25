## This file contains methods for plotting the results in a QSarray object
##
## Author: Christopher Bolen
##         Gur Yaari
## Updated: 2013-03-15
## (c) 2013 Yale University. All rights reserved.


##Function for plotting out the pdfs of all the genes of a gene set
# QSarray1 -- an object of class QSarray.
# QSarray2 -- (optional) The output for a second comparison. 
#                 If provided, the function will add a second plot below the first to compare them.  
# path.index -- Which pathway in the QSarray object(s) to plot the PDF for. This can be of length 1 or 2 to specify different gene sets for the top and bottom plot. If QSarray2 is not specified and path.index is of length 2, the second plot will be drawn from QSarray1 
# normalizePeaks -- a boolean value. If true, curve heights will be normalized to the same value.
# alpha -- specifies the alpha channel for the individual gene density curves.
# addBarcode -- should the barcode plot be added below the PDFs
# barcode.col -- the color used for the bars of the barcode
# barcode.hei = height of the barcode as a fraction of the PDF plot
# groupLabel -- Vector of labels for the individual plots. If plotting two distributions, this must
#               be of length 2.
# labelLoc -- vector of length 1 or 2 determining the location on the plot of where to put the label. One of "left","center", or "right"
# xlab,ylab,main -- parameters to pass on to "plot"
# If you want to plot two curves, one below the other, geneResults2 (and optionally aggregatedData2)
# must be provided. 
plotGeneSetDistributions <- function(QSarray1, QSarray2=NULL, path.index=1, colorScheme="sdHeat",
                                     alpha=1, normalizePeaks=FALSE, addBarcode=TRUE, 
                                     barcode.col=NULL, barcode.hei = 0.2, groupLabel=NULL,labelLoc="left",
                                     xlab="Activity", ylab=NA, main=NULL,lwds=c(1,3),cex=1,...){
  ######################
  ## checking that all inputs are correct
  
  ##check path.index
  if(is.character(path.index)){
    temp = match(path.index, names(QSarray1$pathways))
    if(!is.null(QSarray2)){temp = c(temp, match(path.index, names(QSarray2$pathways)))}
    path.index = temp
  }
  if(!(length(path.index) %in% 1:2)){stop("path.index must be of either length 1 or 2")}
  
  ##get pathway genes
  if(length(path.index)==2){
    pathway1 = QSarray1$pathways[[path.index[1]]]
    pathway2 = QSarray1$pathways[[path.index[2]]]
    if(!is.null(QSarray2)){pathway2 = QSarray2$pathways[[path.index[2]]]}
  }else{
    pathway1 = pathway2 = QSarray1$pathways[[path.index]]
    if(!is.null(QSarray2)){pathway2 = QSarray2$pathways[[path.index]]}
  }
  if(length(pathway1)==0 || length(pathway2)==0){stop("Pathway in path.index has 0 genes")}
  
#   if(is.na(numPathways(QSarray1))){
#     QSarray1 = AggregateGeneSet(QSarray1,PathWayList=list(geneSet=pathway))
#   }
  if(is.null(QSarray2) & length(path.index)==2){
    twoPlots=T
    QSarray2 = QSarray1
    if(is.null(groupLabel)){groupLabel = names(QSarray1$pathways)[path.index]}
  }
  if(!is.null(QSarray2)){
    twoPlots=T
#     if(is.na(numPathways(QSarray2))){
#       QSarray2 = AggregateGeneSet(QSarray2,PathWayList=list(geneSet=pathway))
#     }
    if(is.null(groupLabel)){
      groupLabel = c(deparse(substitute(QSarray1)),deparse(substitute(QSarray2)))
    }
    if(length(groupLabel)!=2){
      if(length(groupLabel)==1 && groupLabel==""){groupLabel = c("","")}
      else{
        stop("Must provide two group labels if plotting two groups")
      }
    }  
  }else{
    twoPlots=F
    if(is.null(groupLabel)){
      groupLabel = c(deparse(substitute(QSarray1)))
    }
  }
  
  ## check color vectors
  if(length(colorScheme)>1 || length(barcode.col)>1){
    if(length(pathway1) != length(pathway2)){
      warning("Pathway 1 and Pathway 2 are different lengths. Using color vectors is not recommended.")
    }
    if(length(colorScheme)>length(pathway1)){
      stop("Length of colorScheme does not match length of pathway stored in QSarray1$pathways.")
    }
    if(length(barcode.col)>length(pathway1)){
      stop("Length of barcode.col does not match length of pathway stored in QSarray1$pathways.")
    }
  }
  
  ## set up labels
  labelLoc = match.arg(labelLoc, c("left","center","right"), several.ok=T)
  if(length(labelLoc)==1){labelLoc = rep(labelLoc,2)}
  
  #########################
  ##set up plotting area
  if(length(barcode.hei)!=1){stop("barcode.hei must be a numeric value of length 1")}
  if(twoPlots){
    if(addBarcode){
      layoutMat = matrix(c(1,4,5,7,6,2,3), ncol=1)
      hei = c(.3,1,barcode.hei,barcode.hei,1,.3,.2)
    }else{
      layoutMat = matrix(c(1,4,5,2,3), ncol=1)
      hei = c(.3,1,1,.3,.2)      
    }
  }else{
    if(addBarcode){
      layoutMat = matrix(c(1,4,5,2,3),ncol=1)
      hei = c(.3,1,barcode.hei,.2,.2)
    }else{
      layoutMat = matrix(c(1,4,2,3),ncol=1)
      hei = c(.3,1,.2,.2)
    }
  }
  layout(layoutMat,heights=hei)
  
  ##add title
  par(mar=c(0,0,0,0)); frame()
  if(!is.null(main)){
    text(0.5,0.5, main, cex=2*cex, font=2)
  }
  
  ##figure out the range for the x axis
  if(!is.null(list(...)$xlim)){
    xRange = list(...)$xlim
  }else{
    xRange = c(min(QSarray1$mean[pathway1]-3*QSarray1$SD[pathway1]),
               max(QSarray1$mean[pathway1]+3*QSarray1$SD[pathway1]))
    
    if(twoPlots){
      xRange2 = c(min(QSarray2$mean[pathway2]-3*QSarray2$SD[pathway2]),
                 max(QSarray2$mean[pathway2]+3*QSarray2$SD[pathway2]))
      xRange = range(c(xRange, xRange2))
    }
  }
  ##add the x axis
  par(mar=c(3,5,0,0))
  plot(0, xlim=xRange,axes=F, type="n",ylab="")
  X<-pretty(xRange,8)
  X<-X [X <= xRange[2] & X >= xRange[1] ]
  axis(1, X, X,cex=cex,...)
  
  ##add the x-label
  par(mar=c(0,5,0,0)); frame()
  text(0.5,0.5, xlab, cex=1.5*cex)
  
  ###############
  ## Generate the individual plots
  for(g in 1:(1+twoPlots)){
    QSarray = get(c("QSarray1","QSarray2")[g])
    pathway = get(c("pathway1","pathway2")[g])
    p.i = ifelse(length(path.index)==1,path.index, path.index[g])
    
    ##x coordinates for the main PDF
    set.x = getXcoords(QSarray, path.index=p.i)
    
    if(length(pathway)!=QSarray$path.size[p.i]){
      warning("length of 'pathway' does not match QSarray$path.size.")
    }
    
    ##plot main PDF
    if(!is.null(list(...)$ylim)){
      yRange = list(...)$ylim
    }else{
      setPeak = max(QSarray$path.PDF[,p.i]) / 
        ifelse(normalizePeaks, max(QSarray$path.PDF[,p.i])*0.25,pdfScaleFactor(QSarray)[p.i])
      
      genePeaks = dt(0,QSarray$dof[pathway])/QSarray$SD[pathway]
      if(normalizePeaks){genePeaks=1}
      gene.q = quantile(genePeaks,probs=0.9)
      
      yRange = c(0, 1.1)*max(gene.q, setPeak)
    }
    if(g==2){yRange = yRange[2:1]}
    par(mar=c(0,5,0,0))
    if(is.na(ylab)){ylab=ifelse(normalizePeaks,"","Density")}
    plot(0, xlim=xRange, ylim=yRange, type="n", axes=F,ylab=ylab,
         cex.axis=par()$cex.axis*cex,cex.lab=par()$cex.lab*cex)
    axis(2,cex=cex,...)
    
    ##add grey box for positive vs negative
    rect(0,-100,xRange[2],100, border=NA, col="#DDDDDD")
    
    ##add group label
    x = switch(labelLoc[g],
               left = xRange[1] + (xRange[2]-xRange[1])*0.05,
               center = xRange[1] + (xRange[2]-xRange[1])*0.5,
               right = xRange[2] - (xRange[2]-xRange[1])*0.05
              )
    text(x,yRange[g%%2+1] - (yRange[g%%2+1]-yRange[g])*0.1,
         groupLabel[g],pos=c(4,2)[1+(labelLoc=="right")], font=2, cex=1.5*cex)
    
    ##plot individual genes in order of their mean
    ord = order(QSarray$mean[pathway])
    
    ##Generate color scheme
    if(class(try(col2rgb(colorScheme),silent=TRUE))!="try-error"){
      ##if they provided colors
      cols = rep(colorScheme, length.out=length(pathway))
      ##if colorScheme is NA, use black
      if(length(colorScheme)==1 && is.na(colorScheme)){cols = rep("#000000", length.out=length(pathway))}
      
    }else if(colorScheme=="rainbow"){
      ##rainbow color scheme
      ##we want colors to be consistent between plots, so we define them for the pathway
      cols = rep(rainbow(10,alpha=alpha),length.out=length(pathway))
      cols = cols[order(ord)]
    } else if(colorScheme=="sdHeat"){
      ##SD Heatmap color scheme
      colFn = colorRamp(c("red","yellow","green","cyan","#00BBFF","#0088FF","#0044FF","#0000FF",
                          "#0000CC","#000099","#000077","#000055","#000033","#000022","#000011","black"))
      #peaks = log(dnorm(0,0,geneResults$data[2,pathway]))
      #peakRange = range(peaks)
      #cols = rgb(colFn((peaks-peakRange[1])/(peakRange[2]-peakRange[1]))/255)
      
      ## I decided I wanted the color range to be fixed, so I'm fixing the range of SDs from 0 to 1 
      ##  and setting anything above 1 to 1
      sds = QSarray$SD[pathway]
      sds[sds>1] = 1
      cols = rgb(colFn(sds)/255)
      cols = paste(cols, sprintf("%x",round(255*alpha)), sep="")
    }

    
    x.seq<-seq(xRange[1],xRange[2],length.out=1000)
    
    ##plot individual gene pdfs
    Max<-10
    for(i in ord){
      DT<-dt(seq(-Max,Max,length.out=10000),QSarray$dof[pathway[i]])
      mean = QSarray$mean[pathway[i]]; SD = QSarray$SD[pathway[i]]
      y.seq<-approx(x=seq(-Max,Max,length.out=10000),y=DT,xout=seq(-Max/SD,Max/SD,length.out=10000),rule=2)$y
      y.seq<-approx(x=seq(-Max,Max,length.out=10000)+mean,y=y.seq,xout=x.seq,rule=2)$y
      y.seq<-y.seq/(sum(y.seq)*(x.seq[2]-x.seq[1]))
#       divisor = ifelse(normalizePeaks, y.seq[findInterval(0,x.seq)],1)      
      divisor = ifelse(normalizePeaks, max((y.seq))/(max(yRange)/1.1/2),1)            
      y.seq<-y.seq/divisor
      lines(x.seq,y.seq,lwd=lwds[1],col=cols[i])
    }
    
    ##main pdf
    divisor = ifelse(normalizePeaks, max(QSarray$path.PDF[,p.i])*0.25,pdfScaleFactor(QSarray)[p.i])
    Y<-approx(set.x, QSarray$path.PDF[,p.i]/divisor,x.seq,rule=2)
    lines(Y,col=1, lwd=lwds[2])
    
    
    ##add barcode plot
    if(addBarcode){
      par(mar=c(0.2,5,0.2,0))
      plot(0, xlim=xRange, ylim=c(0,1), type="n", axes=F,ylab="")
      rect(0,-10,xRange[2],10, border=NA, col="#DDDDDD")
      
      ##line barcode
      if(is.null(barcode.col)){barcode.col="#222222"}
      abline(v=QSarray$mean[pathway],col=barcode.col)
      
      ##beeswarm barcode
#      require(beeswarm)
#      beeswarm(geneResults$data[1,pathway],horizontal=T,add=T,at=0.5,pch=16,cex=1)
    }
  }
  par(mar=c(5,4,4,2),mfrow=c(1,1)) ##reset the plotting window layout, so it's not quite as annoying
}


###A function for plotting the PDFs of a set of pathways. 
## This function uses the data produced by aggregateGeneSet to plot the PDFs of all the pathways in 
## QSarray. 
plotDensityCurves <- function(QSarray,       ##a QSarray object, containing the output of aggregateGeneSet.
                              path.index=1:numPathways(QSarray),    ##which pathways in QSarray to plot
                              zeroLine=TRUE, ##add a vertical line at 0
                              addVIF=!is.null(QSarray$vif),   ##use the VIF when calculating the spread of the pdf
                              col=NULL, ##colors for the curves. Defaults to R's standard colors.
                              plot=TRUE,
                              add=FALSE,  ##Boolean. If false, create a new plotting region to plot the curve to.
                              xlim=NULL,ylim=NULL, ##x- and y-limits of the plotting region.
                              xlab=NULL,ylab=NULL, ##x- and y-labels.
                              type="l",...){                ##additional parameters to pass to "plot.default"   
  
  ##check path.index
  if(is.character(path.index)){
    path.index = match(path.index, names(QSarray$pathways))
  }
  
  if(all(is.na(QSarray$path.mean[path.index]))){stop("no non-missing pathways in path.index")}
  
  scaleFactor = pdfScaleFactor(QSarray,addVIF=addVIF)
  
  if(is.null(xlim)){ xlim=range(calcBayesCI(QSarray,addVIF=addVIF)[,path.index],na.rm=T) 
                     xlim = xlim + (xlim-mean(xlim))*2
                   }
  if(is.null(ylim)){ ylim=range(t(QSarray$path.PDF[,path.index])/scaleFactor[path.index],na.rm=T) }
  
  if(is.null(xlab)){xlab="x"}
  if(is.null(ylab)){ylab="density"}
  
  if(is.null(col)){col=par("col")}

  if(!add & plot){
    plot(0, type="n", xlim=xlim, ylim=ylim,xlab=xlab,ylab=ylab, ...)
  }
  if(length(col)!=length(path.index)){
    col = rep(col, length.out=length(path.index))
  }
  if(zeroLine & plot){abline(v=0,lty=2)}
  
  retVal = list()
  
  for(i in 1:length(path.index)){
    path=path.index[i]
    x = getXcoords(QSarray,path,addVIF=addVIF)
    y = QSarray$path.PDF[,path]/scaleFactor[path]
    if(plot){
      lines(x, y, col=col[i],type=type,...)
    }
    retVal[[i]] = data.frame(x,y)
  }  
  
  ##invisibly return the x- and y- coordinates
  names(retVal) = colnames(QSarray$path.PDF)[path.index]
  invisible(retVal)
}

### Function for plotting the mean and confidence intervals of a set of pathways. 
plotCIs = function(QSarray, 
                   path.index=1:numPathways(QSarray), ##a numeric vector representing which gene sets to plot, and (if sort.by=="none") which order to plot them in. NAs will be removed if sort.by!="none", otherwise an empty column will be added.
                   sort.by=c("mean","p","none"),  ##how the pathways should be ordered.
                   lowerBound=0.025, ## lower bound for confidence interval
                   upperBound=1-lowerBound, ## upper bound for confidence interval 
                   col=NULL,     ##an optional vector indicating the color for the points. If use.p.colors=FALSE is specified, these colors will also be used for the error bars.
                   use.p.colors=TRUE, ##If true, color the points using p-values
                   p.breaks=NULL, ##a vector indicating where the breaks in the p-value color scheme should be. By default, breaks will be at 0.001, 0.005, 0.01, 0.05, & 0.1
                   p.adjust.method = "fdr", ##The method to use to adjust the p-values. Must be one of the methods in p.adjust.methods.
                   addLegend=use.p.colors, ##Should a legend for the p-value color scheme be plotted?
                   
                   lowerColorBar="none", ##Options for plotting a color bar below each point. If either "absolute" or "homogeneity" is specified, a color bar representing the values stored in the QSarray object will be added. If "none" is specified, no bar will be added. A numeric vector the same length as path.index can also be provided, and a color scheme will be created based on the values.
                   lowerColorBar.cols=NULL, ##a vector of colors to be used for the lower color bar.
                   
                   addGrid=TRUE,  ##Should guiding dashed lines be plotted?
                   x.labels=NULL, ##The labels for the individual pathways. By default, names(QSarray$path.mean)
                   cex.xaxis=1, ## set cex parameter manually for x axis label
                   shift=0.0,  ## shift of the x axis: shifts points and arrows (CI's) with respects to the guiding lines and axis labels. Useful when add=TRUE   
                   add=FALSE,    ##boolean parameter. If FALSE, a new plot is created. If TRUE, axes are not plotted
                   ylim=NULL, xlim=NULL,    ##the x- and y-limits of the plotting area
                   ylab=NULL, xlab=NULL,    ##x- and y- axis labels
                   main=NULL, ##plot title
                   sub=NULL, ##plot subtitle
                   type="p", ##plot type. one of "p"-points, "l"-line, or "b"-boths
                   ...){          ##additional parameters to pass to "plot.default"
  
  ##check path.index
  if(is.character(path.index)){
    path.index = match(path.index, names(QSarray$pathways))
  }
  
  if(all(is.na(QSarray$path.mean[path.index]))){stop("no non-missing pathways in path.index")}
  
  
  ###get values from QSarray
  means = QSarray$path.mean[path.index]
  CIs = calcBayesCI(QSarray,low=lowerBound,up=upperBound)[,path.index, drop=FALSE]
  p.vals = pdf.pVal(QSarray,direction=F)[path.index]
  p.vals = p.adjust(p.vals, p.adjust.method)
  p.vals = p.vals*sign(means)
  
  ##x labels
  if(is.null(x.labels)){
    x.labels=names(means)
  }
  
  ##get order of gene sets
  sort.by = match.arg(sort.by)
  if(sort.by=="mean"){ord = order(means, decreasing=T,na.last=NA)}
  else if(sort.by=="p"){ord = order(-log10(abs(p.vals))*sign(p.vals), decreasing=T,na.last=NA)}
  else{
    ord=1:length(means)
  }
  
  ##set plotting parameters
  originalPar = par(no.readonly=T)
  additional_args<-list(...)
  if(!"srt"%in%names(additional_args))par(srt=60)
  #if(!"pch"%in%names(additional_args))par(pch="x")
  par(...)

  ##color of the points
  if(!is.null(col)){   ##check that col is the correct length.  
    if(length(means) %% length(col) !=0 ){warning("length of col is not a multiple of input length")}
    col = rep(col, length.out=length(means))
  }else{
    col = rep(par("col"),length.out=length(means))
  }
  
  ###establish p-value color scheme
  if(is.null(p.breaks)){
    p.breaks = c(0.001,0.005,0.01,0.05,0.1)
  }else{
    p.breaks=abs(p.breaks)
    p.breaks=as.numeric(names(table(p.breaks)))
    p.breaks=p.breaks[p.breaks < 1 & p.breaks > 0]
    p.breaks = p.breaks[order(p.breaks)] ##make sure they're in ascending order
  }
  
  br.ln = length(p.breaks)
  p.breaks.twoSided<-c(-1,-rev(p.breaks),0,p.breaks,1)
  
  if(use.p.colors){
    p.colorScheme<-c(rgb(0,seq(0,1,length.out=br.ln+1),0),
                     rgb(seq(1,0,length.out=br.ln+1),0,0))
    bar.col = p.colorScheme[findInterval(p.vals,p.breaks.twoSided,rightmost.closed=T)]
    bar.col[p.vals==0] = c("#00FF00","#FF0000")[(means>0)+1][p.vals==0]
  }else{
    bar.col=col
  }
  
  if(is.null(xlab)){xlab=""}
  if(is.null(ylab)){ylab="Pathway Activity"}
  
  ##make plotting region
  if(!add){
    if(is.null(ylim)){ylim = range(CIs,na.rm=T)}
    plot(means[ord],type="n",las=1, ylim=ylim,  axes=FALSE,xlab=xlab,ylab=ylab,main=main,...)
    ##add guiding lines
    if(addGrid){abline(v=1:length(ord),col=gray(seq(0.5,1,length.out=10)),lty=2)}
    
    ##add axes
    axis(2,las=1,...)
    if(is.null(list(...)$xaxt) || list(...)$xaxt!="n"){
      Oldcex=par('cex')
      par(cex=0.5*Oldcex*cex.xaxis)
      #         axis(1, at=ord,las=2, labels=x.labels,...)
      axis(1, at=ord,las=2, labels=rep("",length(ord)),...)
      Ys<-par("usr")[3] -  par()$cxy[2]*par()$mgp[2]
      text(1:length(ord), Ys, adj = 1, labels = x.labels[ord], xpd = TRUE,...)
      par(cex=Oldcex)
    }
    abline(h=0,lty=2)
  }
  
  ##plot data
  arrows(1:length(ord)+shift, CIs[2,ord], 1:length(ord)+shift, CIs[1,ord],
         code=3,length=0.1,angle=90, col=bar.col[ord])
  points(1:length(ord)+shift,means[ord],type=type,col=col[ord],...)
  ##get x-axis labels
  if(is.null(x.labels))
    x.labels=names(QSarray$path.mean)[path.index]  
  
  
  ## add legend
  if(use.p.colors && addLegend && !add){
    ##the original color scheme 
    p.colorScheme<-c(rgb(0,seq(1,0,length.out=br.ln+1),0),
                     rgb(seq(0,1,length.out=br.ln+1),0,0))
    p.labels<-round(c(-p.breaks,1,rev(p.breaks)),3)
    n = length(p.labels)
    
    usr = par('usr')
    strsize = strwidth("W")  ##I need the height of a string oriented vertically, and it just happens that 'W' is approximately as wide as it is tall (in the default font).
    yvalues<-usr[4]-c(0,strheight("W")*1.5)
    xvalues<-usr[2]-(seq(strsize*(n+2),0,length.out=n+2)*1.2)
    
    text(mean(xvalues), yvalues[0]+strheight("X"), labels="P-values",pos=3,xpd=T)
    ##add boxes
    for(j in 1:(n+1)){
      polygon(c(xvalues[j],xvalues[j+1],xvalues[j+1],xvalues[j]),
              c(yvalues[1],yvalues[1],yvalues[2],yvalues[2]),
              col=p.colorScheme[j],border=p.colorScheme[j])
    }
    ##add labels
    for(j in 1:n){
      text(xvalues[j+1],yvalues[2]-strheight("W"),p.labels[j],adj=1,srt=90)
    }
  }
  
  ## add lower color bar
  lowBar.vals = NULL
  if(!is.character(lowerColorBar)){        ##### If plotting a random vector
    if(length(lowerColorBar)!=length(path.index)){
      warning("lowerColorBar not the same length as path.index. Color bar omitted")
    }else{
      lowBar.vals = lowerColorBar
      
      ##color scheme
      if(is.null(lowerColorBar.cols)){
        lowerColorBar.cols = rgb(seq(1,0,length.out=6),seq(1,0,length.out=6),1)
      }
      lowBar.breaks = seq(range(lowBar.vals)[1],range(lowBar.vals)[2],length.out=length(lowerColorBar.cols)+1)
    }
  }
  else if(lowerColorBar=="absolute"){     ###### If plotting absolute GST p-value
    if(is.null(QSarray$absolute.p)){
      warning("Absolute p-values not found. Color bar omitted.")
    }else{
      lowBar.vals = QSarray$absolute.p[path.index]
      lowBar.vals<-p.adjust(abs(lowBar.vals),method=p.adjust.method)
      
      ##color scheme
      lowBar.breaks = p.breaks.twoSided
      if(is.null(lowerColorBar.cols)){
        lowerColorBar.cols = c(rgb(0,seq(0,1,length.out=br.ln+1),0,0),
                               rgb(seq(1,0,length.out=br.ln+1),0,0))
      }
    }
  }
  else if(lowerColorBar=="homogeneity"){   ###### If plotting homogeneity score
    if(is.null(QSarray$homogeneity)){
      warning("Homogeneity score not found. Color bar omitted.")
    }else{
      lowBar.vals = QSarray$homogeneity[path.index]
      
      ##color scheme
      if(is.null(lowerColorBar.cols)){
        lowerColorBar.cols = rgb(seq(1,0,length.out=6),seq(1,0,length.out=6),1)
      }
      lowBar.breaks = seq(range(lowBar.vals)[1],range(lowBar.vals)[2],length.out=length(lowerColorBar.cols)+1)
    }
  }
  if(!is.null(lowBar.vals)){
    col = lowerColorBar.cols[findInterval(lowBar.vals,lowBar.breaks,rightmost.closed=T)]
    
    ###Plot box
    for(i in 1:length(ord)){
      usr = par()$usr;
      polygon(i+c(-0.5,-0.5,0.5,0.5), usr[3]+c(0.01,0.03,0.03,0.01)*(usr[4]-usr[3]),
              col=col[ord][i],border=col[ord][i])
    }
  }
  
  box(...)
  par(originalPar[c("srt","pch",names(list(...)))])
}





##Function for plotting out the pdfs of all the genes of a gene set

plotCIsGenes <- function(QSarray, ##an object of class QSarray.
                         path.index=1, ##the gene set name
                         gene.list=NULL,  ##Character vector. A predefined list of genes in the gene set to be plotted. If sort.by="none", the order of these genes will be used. NAs are accepted.
                         sort.by=NULL,  ##one of c(mean,p,none), specifying how to order the genes. If NULL and gene.list is provided, default is "none".
                         lowerBound=0.025, ## lower bound for confidence interval
                         upperBound=1-lowerBound, ## upper bound for confidence interval 
                         asBand=FALSE, ### plot CI as a grey band or as arrows
                         
#                          p.color.method="none", ##one of c("none","absolute","homogeneity")
#                          p.breaks=NULL, ##a vector indicating where the breaks in the p-value color scheme should be. By default, breaks will be at 0.001, 0.005, 0.01, 0.05, & 0.1
#                          p.adjust.method = "fdr", ##The method to use to adjust the p-values. Must be one of the methods in p.adjust.method.
#                          addLegend=p.color.method!="none", ##Should a legend for the p-value color scheme be plotted?
                         
                         col=NULL,     ##default color parameter to pass to "plot". If provided, it overrides "p.color.method"
                         addGrid=TRUE,  ##Should guiding dashed lines be plotted?
                         x.labels=NULL, ##The labels for the individual genes By default, names(QSarray$mean[gene.names]) or ORDER if not empty
                         cex.xaxis=1, ## set cex parameter manually for x axis label
                         shift=0.0,  ## shift of the x axis: shifts points and arrows (CI's) with respects to the guiding lines and axis labels. Useful when add=TRUE   
                         pathwayCI=c("band","bar","none"), ##A string, one of "band", "bar", or "none", determining whether to add the confidence interval for the gene set PDF to the plot. If "band", a band will be plotted behind the bars for the individual genes. If "bar" is specificied, another error bar will be added before the genes' error bars.
                         meanCol=4, ## color for the line that indicates the mean of the pathway
                         
                         add=FALSE,    ##boolean parameter. If FALSE, a new plot is created. If TRUE, axes are not plotted
                         ylim=NULL, xlim=NULL,    ##the x- and y-limits of the plotting area
                         ylab=NULL, xlab=NULL,    ##x- and y- axis labels
                         main=NULL, ##plot title
                         sub=NULL, ##plot subtitle
                         ...){          ##additional parameters to pass to "plot.default"
  
  ##check path.index
  if(is.character(path.index)){
    path.index = match(path.index, names(QSarray$pathways))
  }
  
  pathwayCI = match.arg(pathwayCI)
#   p.color.method=match.arg(p.color.method, c("none","absolute","homogeneity"))
  
  ###get values from QSarray
#   if(sum(!gene.names%in% names(QSarray$mean))){
#     warning("some gene names do not appear in QSarray object");
#     gene.names<-gene.names[gene.names%in% names(QSarray$mean)]
#   }
  if(sum(!is.na(path.index))==0){stop("pathway provided is not in QSarray")}
  if(sum(!is.na(path.index))>1){stop("more than one pathway provided")}
  gene.names<-names(QSarray$mean)[QSarray$pathways[[path.index]]]
  if(!is.null(gene.list)){
    gene.names = gene.list
  }
  
  if(is.na(QSarray$path.mean[path.index])){stop("pathway provided has length of 0")}
  
  
  means = QSarray$mean[gene.names]
  #ExistingGenes = 1:length(gene.names)
  #if(is.null(gene.list)){ExistingGenes=which(gene.list%in%gene.names)}
  
  ##calculate CIs for each gene
  CIs = sapply(gene.names,function(NAME){
    sd = QSarray$SD[NAME]*QSarray$sd.alpha[NAME]
    t.ci = qt(c(lowerBound,upperBound),QSarray$dof[NAME])
    return(QSarray$mean[NAME]+sd*t.ci)
  })
  #if(!length(dim(CIs)))CIs<-matrix(CIs)
  
  
  ##Get P-values for each gene
#   p.vals = sapply(gene.names,function(NAME)pt((QSarray$mean[NAME]-QSarray$path.mean[pathway.name])/(QSarray$SD[NAME]*QSarray$sd.alpha[NAME]),QSarray$dof[NAME]))
#   p.vals [ p.vals > 0.5 ] = 1-p.vals [ p.vals > 0.5 ]
#   p.vals = 2 * p.vals
#   p.vals = p.adjust(p.vals, p.adjust.method)
#   p.vals = p.vals*sign(means-QSarray$path.mean[pathway.name])
#   p.vals=NULL
#   if(p.color.method %in% c("absolute","homogeneity")){
#     p.vals = absoluteTest.genePvals(QSarray,path.index,compareTo=ifelse(p.color.method=="absolute","z","p"))[[1]]
#     p.vals = p.adjust(abs(p.vals), p.adjust.method)
#     p.vals = p.vals*sign(means-QSarray$path.mean[path.index])
#     p.vals = p.vals[gene.names]
#   }
  
  ##get order of gene sets
  if(is.null(sort.by)){sort.by = c("mean","none")[2-is.null(gene.list)]}
  sort.by = match.arg(sort.by, c("mean","p","none"))
  if(sort.by=="mean"){
    ord = order(means, decreasing=T,na.last=NA)
  }else if(sort.by=="p"){
    sd = QSarray$SD*QSarray$sd.alpha
    t.stat = QSarray$mean/(sd/sqrt(QSarray$dof+2))
    p.val = dt(t.stat,QSarray$dof)
    ord = order(-log10(p.val[gene.names])*sign(t.stat), decreasing=T, na.last=NA)
  }else{
    ord=1:length(means)
  }
  
  ##set up plotting region
  originalPar = par(no.readonly=T)
  additional_args<-list(...)
  if(!"srt"%in%names(additional_args))par(srt=60)
  if(!"pch"%in%names(additional_args))par(pch="x")
  par(...)
  
  ###establish p-value color scheme
  if(!is.null(col)){   ##if they specify "col", we won't worry about a p-value color scheme.  
    if(length(ord) %% length(col) !=0 ){warning("length of col is not a multiple of input length")}
    col = rep(col, length.out=length(ord))
  }else{
#     if(is.null(p.breaks)){
#       p.breaks = c(0.001,0.005,0.01,0.05,0.1)
#     }else{
#       p.breaks = abs(p.breaks)
#       p.breaks = as.numeric(names(table(p.breaks)))
#       p.breaks = p.breaks[p.breaks < 1 & p.breaks > 0]
#       p.breaks = p.breaks[order(p.breaks)] ##make sure they're in ascending order
#     }
#     
#     br.ln = length(p.breaks)
#     p.breaks.twoSided<-c(-1,-rev(p.breaks),0,p.breaks,1)
#     
#     if(!is.null(p.vals)){
#       p.colorScheme<-c(rgb(0,seq(0,1,length.out=br.ln+1),0),
#                        rgb(seq(1,0,length.out=br.ln+1),0,0))
#       col = p.colorScheme[findInterval(p.vals,p.breaks.twoSided,rightmost.closed=T)]
#     }else{
       col=rep(par("col"),length.out=length(means))
#     }
   }
  
  if(is.null(xlab)){xlab=""}
  if(is.null(ylab)){ylab="Gene Activity"}
  if(is.null(main)){main=names(QSarray$pathways)[path.index]}
  

  ##plot data
  if(!add){
    if(is.null(ylim)){ylim = range(CIs,na.rm=T)}
    xlim=c(1,length(ord))
    if(pathwayCI=="bar")xlim=c(0,length(ord))
    plot((1:length(ord)),means[ord],type="n",las=1, ylim=ylim,  axes=FALSE,xlab=xlab,ylab=ylab,xlim=xlim, main=main,...)
    
    ##add guiding lines
    if(addGrid){
      abline(v=1:length(ord),col=gray(seq(0.5,1,length.out=10)),lty=2)
    }
  }
  
  ## confident band for pathway PDFs
  if(pathwayCI=="band"){
    CI_Path = calcBayesCI(QSarray,low=lowerBound,up=upperBound)[,path.index]
    usr = par()$usr
    polygon(c(usr[1],usr[1],usr[2],usr[2]),c(CI_Path[1],CI_Path[2],CI_Path[2],CI_Path[1]),border=grey(0.9),col=grey(0.9))
    abline(h=QSarray$path.mean[path.index],lty=4,col=meanCol)
  }

  
  if(pathwayCI=="bar"){
    CI_Path = calcBayesCI(QSarray,low=lowerBound,up=upperBound)[,path.index]
    points(shift,QSarray$path.mean[path.index],...)
    arrows(shift, CI_Path[1], shift, CI_Path[2],
        code=3,length=0.1,angle=90, col=1,lwd=2)
    abline(h=QSarray$path.mean[path.index],lty=4,col=meanCol)
  }
  
  if(!asBand){
      points(1:length(ord)+shift, means[ord], type="p", col=col[ord],...)
      arrows(1:length(ord)+shift, CIs[2,ord], 1:length(ord)+shift, CIs[1,ord],
       code=3,length=0.1,angle=90, col=col[ord])
  }
  else{
      polygon(c(1:length(ord)+shift,length(ord):1+shift),c(CIs[2,ord],CIs[1,rev(ord)]),border=grey(0.8),col=grey(0.8))
      points(1:length(ord)+shift,means[ord],type="l",col=col[ord])
  }

  ##get x-axis labels
  if(is.null(x.labels))
    x.labels=gene.names[ord]
  ##add axes
  if(!add){
    axis(2,las=1)
    if(is.null(list(...)$xaxt) || list(...)$xaxt!="n"){
        Oldcex=par('cex')
        par(cex=0.5*Oldcex*cex.xaxis)
#         axis(1, at=ord,las=2, labels=x.labels,...)
         axis(1, at=ord,las=2, labels=rep("",length(ord)),...)
         
         Ys<-par("usr")[3] -  par()$cxy[2]*par()$mgp[2]
        text(1:length(x.labels), Ys, adj = 1,
          labels = x.labels, xpd = TRUE,...)
        if(pathwayCI=="bar"){axis(1, at=0,las=2, labels="",...)}
        if(pathwayCI=="bar"){
          text(0, Ys, adj = 1,labels = names(QSarray$pathways)[path.index], xpd = TRUE,...) 
        }
        par(cex=Oldcex)
    }
    abline(h=0,lty=1,lwd=0.5)   
  }
  ## add legend
#   if(!is.null(p.vals) && addLegend){
#     ##the original color scheme 
#     p.colorScheme<-c(rgb(0,seq(1,0,length.out=br.ln+1),0),
#                      rgb(seq(0,1,length.out=br.ln+1),0,0))
#     p.labels<-round(c(-p.breaks,1,rev(p.breaks)),3)
#     n = length(p.labels)
#     
#     usr = par('usr')
#     strsize = strwidth("W")  ##I need the height of a string oriented vertically, and it just happens that 'W' is approximately as wide as it is tall (in the default font).
#     yvalues<-usr[4]-c(0,strheight("W")*1.5)
#     xvalues<-usr[2]-(seq(strsize*(n+2),0,length.out=n+2)*1.2)
#     
#     text(mean(xvalues), yvalues[0]+strheight("X"), labels="P-values",pos=3,xpd=T)
#     ##add boxes
#     for(j in 1:(n+1)){
#       polygon(c(xvalues[j],xvalues[j+1],xvalues[j+1],xvalues[j]),
#               c(yvalues[1],yvalues[1],yvalues[2],yvalues[2]),
#               col=p.colorScheme[j],border=p.colorScheme[j])
#     }
#     ##add labels
#     for(j in 1:n){
#       text(xvalues[j+1],yvalues[2]-strheight("W"),p.labels[j],adj=1,srt=90)
#     }
#   }
   ###Plot box
  box(...)
  par(originalPar[c("srt","pch",names(list(...)))])
}







