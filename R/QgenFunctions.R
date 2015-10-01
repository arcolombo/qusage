
#Simple function to check if a design matrix from a linear model is well conditioned.
#This is needed for the adjustment to the residuals if a random effect is present with few number of replicates
#Sometimes, depending on the model, the random effect can't be treated as a fixed effect as the design matrix of
#fixed effects becomes singular.  For longitudinal studies this is the case when subjects are nested within a single
#level of another fixed effect factor

nonsingular.check <- function(m) class(try(solve(m),silent=T))=="matrix"




#The qgen function.



qgen<-function(eset,design,fixed,geneSets,contrast.factor,contrast,
               random=NULL,correlation=NULL,design.sampleid=NULL){
  
  
  fixed=paste(as.character(fixed),sep="",collapse="")
  
  #If sample id variable is given we will do additional checks to make sure the expression
  #and design file match up. Otherwise it is assumed that the order of the samples in the 
  #design file is the same as the columns of eset.
  if(!is.null(design.sampleid)){
    #Removing samples in eset that do not correspond with the design and vice versa
    if(all(design[,design.sampleid] %in% colnames(eset))){
      cat("All samples present in design")
    } else{
      cat("Descrepancy of samples between eset and design. Samples not included in both are removed.")
    }
    if(all(colnames(eset) %in% design[,design.sampleid])){
      cat("\nSamples in Design and Expression Match")
    }else{
      cat("\nDescrepancy of samples between eset and design.  Samples not included in both are removed.")
    }
    #ifelse(all(design[,design.sampleid] %in% colnames(eset)),cat("All samples present in design"),cat("Descrepancy of samples between eset and design. Samples not included in both are removed."))
    #ifelse(all(colnames(eset) %in% design[,design.sampleid]),cat("Samples in Design and Expression Match"),cat("Descrepancy of samples between eset and design.  Samples not included in both are removed."))
    
    #Keep samples that are included in both eset and design
    index<-which(colnames(eset) %in% design[,design.sampleid])
    eset<-eset[,index]
    index2<-which(design[,design.sampleid] %in% colnames(eset))
    design<-design[index2,]
    
    #Order Design and expression so they are compatible for model merging
    eset<-eset[,order(colnames(eset))]
    design<-design[order(design[,design.sampleid]),]
    
  }
  
  #If samples have no sample id it is assumed they are in correct order.  We simply check 
  #if the row and column numbers match, as that is all we can do to find if there is an issue.
  
  if(is.na(design.sampleid)){
    warning("No sample identifier given. Rows in design are assumed to correspond to the columns of eset.")
    if(dim(eset)[2]!= dim(design)[1]){
      stop("Number of rows in design do not match number of samples in eset.") 
    }  
  }
  
  
  #Removing rows in eset that are not in the genesets to save computational time
  geneindex<-which(rownames(eset) %in% unlist(geneSets))
  eset<-eset[geneindex,]
  
  
  
  #Converting contrast information into actual contrast vector and checking factor name compatibility  
  contrast.factor.names<-gsub(" ","",unlist(strsplit(as.character(contrast.factor),split="*",fixed=T))[-1])
  
  contrast.factor.2<-vector("list", length(contrast.factor.names)) 
  for(i in 1:length(contrast.factor.names)){
    contrast.factor.2[[i]]<-levels(design[,contrast.factor.names[i]])
  }
  new.factor.levels<-do.call(paste,c(do.call(expand.grid,contrast.factor.2),sep=""))
  if(!all(unlist(lapply(design[,contrast.factor.names],is.factor)))){
    stop("Variable included in contrast.factor is not of type 'factor'. ")
  }
  if(!all((new.factor.levels)==make.names(new.factor.levels,unique=TRUE))){
    stop("Factor level combinations specified in contrast.factor do not have valid R names.")
  } 
  contrast<-list(comparison=as.vector(do.call(makeContrasts,list(contrasts=contrast,levels=new.factor.levels))))
  
  
  
  #Determining which modeling function to use GLS or LME and what type of residuals should be calculated for qusage.
  ModelInd<-0
  ResidualMethod="Conditional"
  
  if(!is.null(random)){
    ModelInd<-1
    #checking if the random effect is "more fancy" than a random effects regarding a single repeated measure
    #aka no nested random effects. Also determining if the the number of replicates is low to take the addition residual
    #approach form the Q-gen paper
    
    val<-length(grep("/",as.character(random)[2]))
    if(val>0){ResidualMethod="Conditional"}
    if(val==0){
      rep.var<-strsplit(as.character(random)[2],"| ",fixed=T)[[1]][2]
      if(median(table(design[,rep.var]))>4){ResidualMethod="Conditional"}
      if(median(table(design[,rep.var]))<5){
        mm<-model.matrix(formula(paste(fixed,rep.var,sep="+")),data=design)
        ResidualMethod=ifelse(nonsingular.check(t(mm)%*%mm),"Adjusted","Conditional")
        if(ResidualMethod=="Conditional"){
          warning("Number of random effect replicates are low.  Adjusted Vif estimation can not be conducted due to nested factors. The conditional residuals will be used and could potentially be overly optimistic.")
        }
      }
    }
    
    
    
  }
  
  #Running row by row modeling using GLS or LME depending if random effects are in the model or not
  
  #Initializing the main components needed for qusage analyis, logFc, Standard error, and d.o.f. 
  #as well as the residual matrix
  #lmm.result<-vector("list",dim(eset)[1])
  
  mean<-vector("numeric",length=dim(eset)[1])
  SD<-vector("numeric",length=dim(eset)[1])
  dof<-vector("numeric",length=dim(eset)[1])
  names(mean)<-rownames(eset)
  names(SD)<-rownames(eset)
  names(dof)<-rownames(eset)
  converge.ind<-c()
  
  residual.matrix<-matrix(, nrow = nrow(eset), ncol = ncol(eset))
  if(ResidualMethod=="Adjusted"){
    cat("Number of random effect replicates are low.  The adjusted residuals will be used for VIF estimation")
  }
  
  cat("\nRunning row by row modeling.")
  cat("\nPercent Complete:")
  percent=0.1
  #GLS
  if(ModelInd==0){
    for (i in 1:nrow(eset)){
      if((i/nrow(eset))>percent){
        cat(paste(percent*100,"%...",sep="")) 
        percent<-percent+.1
      }
      
      design$y<-eset[i,]                        
      modresult<-try(do.call(gls,
                         list(model=formula(paste("y",fixed,sep="")),data=design,
                              correlation=correlation,control=glsControl(opt='optim'))
                         ),
                     silent=T
                     )
      if(class(modresult)=="try-error"){
        if(grepl("converg", modresult)){
          converge.ind<-c(converge.ind,i)
        }else{
          stop(modresult)
        }
      }else{
        contrast.result<-summary(contrast(lsmeans(modresult,contrast.factor),contrast,by=NULL))
        mean[i]<-contrast.result$estimate
        SD[i]<-contrast.result$SE
        dof[i]<-contrast.result$df
        residual.matrix[i,]<-residuals(modresult)
      }
    }
  }
  
  
  #LME
  if(ModelInd==1){
    for (i in 1:nrow(eset)){
      if((i/nrow(eset))>percent){
        cat(paste(percent*100,"%...",sep="")) 
        percent<-percent+.1
      }
      design$y<-eset[i,]
      modresult<-try(lme(formula(paste("y",fixed,sep="")),data=design,random=random,
                         correlation=correlation,control=lmeControl(opt='optim')),
                     silent = T)
      if(class(modresult)=="try-error"){
        if(grepl("converg", modresult)){
          converge.ind<-c(converge.ind,i)
        }else{
          stop(modresult)
        }
      }else{
      
        contrast.result<-summary(contrast(lsmeans(modresult,contrast.factor),contrast,by=NULL))
        
        mean[i]<-contrast.result$estimate
        SD[i]<-contrast.result$SE
        dof[i]<-contrast.result$df
        
        if(ResidualMethod=="Conditional"){residual.matrix[i,]<-residuals(modresult)}
        if(ResidualMethod=="Adjusted"){
          residual.matrix[i,]<-residuals(lm(formula(paste("y",fixed,"+",rep.var,sep="")),data=design))
        }
      }
    }
  }
  
  
  
  rownames(residual.matrix)<-rownames(eset)
  colnames(residual.matrix)<-rep("Resid",dim(eset)[2])
  
  if(length(converge.ind)>0){
    residual.matrix<-residual.matrix[-converge.ind,]
    mean<-mean[-converge.ind]
    SD=SD[-converge.ind]
    dof=dof[-converge.ind]
    warning(paste(length(converge.ind)," gene(s) for the linear model did not converge and were removed. The first omission is found at  ",rownames(eset)[converge.ind[1]],sep=""))
  }
  
  results<-newQSarray(list(mean=mean,SD=SD,dof=dof,labels=colnames(residual.matrix)))
  #final.result<-list(residual.matrix=residual.matrix,QSobject=QSobject,res.method=ResidualMethod)
  
  overlap = sapply(geneSets, function(x){ sum(x %in%
                                                rownames(eset)) })
  geneSets=geneSets[overlap>0]
  
  
  
  cat("\nAggregating gene data for gene sets.")
  results<-aggregateGeneSet(results, geneSets,n.points=2^14)
  cat("Done. \nCalculating VIF's on residual matrix.")
  results<-calcVIF(residual.matrix,results,useCAMERA=FALSE)
  cat("\nQ-Gen analysis complete.")
  results
}