n_reps  = 1500

#some guys list combine thing from
# http://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


my.unlist<-function(data){
  return(list(data[,1],data[,2]))
}

plot.with.CI <- function(inVec,plotDataUpList,plotDataDownList){
  
  inVec<-as.data.frame(t(as.numeric(inVec)))
  names(inVec)<-header
  

  tauUpLine<-sort(c(0,seq(0,max(tauUp1),1/20),downStartTimes)) #19 Oct 2015 - added extra '0' to this list. Problem is that calcFunctionXXPop sets first two timepoints to have B+K+ <-0 because it assumes it's working with the data which has two time points at t = 0, but fucks up if you don't (only plots). This extra 0 means this truncation is correct.
  
  #recalcualte position in upTau of downstart times (use this to cheat a bit when calculating initial conditions - not so important here when only called once but meh!)
  downStartUpPosLine<-list()
  for(i in 1:nDown){
    tau<-eval(join(c(tauL,i)))
    downStartUpPosLine[[i]]<-match(min(tau),tauUpLine)
  }
  downStartUpPosLine<-unlist(downStartUpPosLine)
  
  tauDownLineList<-list()
  for(i in 1:nDown){
    tau<-tauDownList[[i]]
    tauDownLineList[[i]]<-seq(min(tau),max(tau),(max(tau)-min(tau))/20)
  }

  ############################################################################################################
  if(model=="kinHet"||model=='kinHetExtended'||model=='kinHetExtended2'){
    
    pars<-calcParsKinHet(inVec)
    parsA<-pars[[1]]
    parsB<-pars[[2]]
    
    best.fit.vals<-calcFunctionDoublePop(parsA,parsB,tauUpLine,tauDownLineList,downStartUpPosLine)
    
    boot.folder<-paste(outputFolder,"/",cell,"/bootstrap/B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"/",sep="")
  
    boot.frac.list<-rep( list(list()), 2*(nUp+nDown ) ) 
 
    boot.vals<-foreach(j=1:n_reps,.combine='comb', .multicombine=TRUE,.init=boot.frac.list) %do% {   .export=c('bBoxMax','kBoxMax','sigV','model','source.switch','calcParsKinHet','header','nuV')
      library(xlsx) 
      boot.file<-paste(boot.folder,"B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"Replicate",toString(j),"bestfit.xlsx",sep="")
      vals<-read.xlsx(boot.file,sheetIndex=1)
      import.header<-unlist(lapply(vals[[1]],as.character)) #only needs to be defined once but meh, if statement just as costly
      vals<-as.data.frame(rbind(vals[[2]]))
      names(vals)<-import.header
      all.pars<-calcParsKinHet(vals[header])
      boot.parsA<-all.pars[[1]]
      boot.parsB<-all.pars[[2]]
      
      boot.function.vals<-calcFunctionDoublePop(boot.parsA,boot.parsB,tauUpLine,tauDownLineList,downStartUpPosLine)
      boot.fracs<-lapply(boot.function.vals,computeFinalFracs)
      for(kk in 1:nUp){
        boot.fracs[[kk]][1:2,]<-0 
      }    
      rm(kk)

      output<-foreach(i2 = 1:(nUp+nDown )) %do% { #only four slopes but each has 2 curves for the two brduki67 pops
        my.unlist(boot.fracs[[i2]])
      } 
      
      output<-unlist(output,recursive=F)
      return(output)#returns list of lists (1 for each slope+Ki67Brdu+- pop, so 8 total)
    }     

    all.curve.CI<-foreach(i3 = 1:(2*(nUp+nDown ))) %do% {
       temp.frame<-do.call(cbind,boot.vals[[i3]])
       temp.CI<-t(apply(temp.frame, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE))
       colnames(temp.CI)<-c('lower','upper')
       temp.CI
       
    }

  }  
  
  ############################################################################################################
  if(model=="tempHet"){
    
    pars<-calcParsTempHet(inVec)
    
    best.fit.vals<-calcFunctionSinglePop(pars,tauUpLine,tauDownLineList,downStartUpPosLine)
    
    boot.folder<-paste(outputFolder,"/",cell,"/bootstrap/B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"/",sep="")
    
    boot.frac.list<-rep( list(list()), 2*(nUp+nDown ) ) 
    
    boot.vals<-foreach(j=1:n_reps,.combine='comb', .multicombine=TRUE,.init=boot.frac.list, .export=c("bBoxMax","kBoxMax","sigV","model","source.switch","phiGlobal","calcParsTempHet","header","nuV")) %do% {  
      library(xlsx) 
      library(nleqslv)
      boot.file<-paste(boot.folder,"B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"Replicate",toString(j),"bestfit.xlsx",sep="")
      vals<-read.xlsx(boot.file,sheetIndex=1)
      import.header<-unlist(lapply(vals[[1]],as.character)) #only needs to be defined once but meh, if statement just as costly
      vals<-as.data.frame(rbind(vals[[2]]))
      names(vals)<-import.header
      all.pars<-calcParsTempHet(vals[header])
      
      boot.function.vals<-calcFunctionSinglePop(all.pars,tauUpLine,tauDownLineList,downStartUpPosLine)
      boot.fracs<-lapply(boot.function.vals,computeFinalFracs)
      for(kk in 1:nUp){
        boot.fracs[[kk]][1:2,]<-0 
      }    
      rm(kk)
      output<-foreach(i2 = 1:(nUp+nDown )) %do% { #only four slopes but each has 2 curves for the two brduki67 pops
      my.unlist(boot.fracs[[i2]])
      } 
      output<-unlist(output,recursive=F)
      return(output)#returns list of lists (1 for each slope+Ki67Brdu+- pop, so 8 total)
    }     
    
    all.curve.CI<-foreach(i3 = 1:(2*(nUp+nDown ))) %do% {
      temp.frame<-do.call(cbind,boot.vals[[i3]])
      temp.CI<-t(apply(temp.frame, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE))
      colnames(temp.CI)<-c('lower','upper')
      temp.CI
    }
    
  }

  ############################################################################################################

  fracs<-lapply(best.fit.vals,computeFinalFracs)
  for(kk in 1:nUp){
    fracs[[kk]][1:2,]<-0 
  }    
  rm(kk)  
  
  finalFracListKm <- list()
  finalDataListKm <- list()
  finalFracListKp <- list()
  finalDataListKp <- list()
  list() 
  upFracs<-fracs[seq(1:nUp)] 
  
  for(i in 1:nUp){
    finalFracListKm[[i]]<-data.frame(t=tauUpLine,kMinusFrac=upFracs[[i]][,1],lower=all.curve.CI[[2*i-1]][,1],upper=all.curve.CI[[2*i-1]][,2]) #super inefficient to keep calling computeFinalFracs, but okay for now as a hack
    finalDataListKm[[i]]<-data.frame(t=tauUp1,kMinusFrac=plotDataUpList[[i]][,1])
    
    finalFracListKp[[i]]<-data.frame(t=tauUpLine,kPosFrac=upFracs[[i]][,2],lower=all.curve.CI[[2*i]][,1],upper=all.curve.CI[[2*i]][,2])
    finalDataListKp[[i]]<-data.frame(t=tauUp1,kPosFrac=plotDataUpList[[i]][,2])
  }
  
  
  downFracs<-fracs[-seq(1:nUp)] 

  for(i in 1:nDown){

    finalFracListKm[[i+nUp]]<-data.frame(t=tauDownLineList[[i]],kMinusFrac=downFracs[[i]][,1],lower=all.curve.CI[[2*nUp+2*i-1]][,1],upper=all.curve.CI[[2*nUp+2*i-1]][,2])#super inefficient to keep calling computeFinalFracs, but okay for now as a hack
    finalDataListKm[[i+nUp]]<-data.frame(t=tauDownList[[i]],kMinusFrac=plotDataDownList[[i]][,1])

    finalFracListKp[[i+nUp]]<-data.frame(t=tauDownLineList[[i]],kPosFrac=downFracs[[i]][,2],lower=all.curve.CI[[2*nUp+2*i]][,1],upper=all.curve.CI[[2*nUp+2*i]][,2])
    finalDataListKp[[i+nUp]]<-data.frame(t=tauDownList[[i]],kPosFrac=plotDataDownList[[i]][,2])
  }
 


  for(i in 1:length(finalFracListKm)){
    colnames( finalFracListKm[[i]])<-c("t","best","lower","upper")
    colnames( finalDataListKm[[i]])<-c("t","best")
  }

  for(i in 1:length(finalFracListKp)){
    colnames( finalFracListKp[[i]])<-c("t","best","lower","upper")
    colnames( finalDataListKp[[i]])<-c("t","best")
  }

  names(finalFracListKm)<-labelsKm
  names(finalDataListKm)<-labelsKm
  
  names(finalFracListKp)<-labelsKp
  names(finalDataListKp)<-labelsKp

  ############################################################################################################
  cols <- c("blue", "red","green")
  global.theme<-theme_bw()+theme(legend.position="none")#theme(axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=12)) +theme(plot.title = element_text(family = "Helvetica", color="#666666", face="bold", size=15, hjust=0)) 

  combinedLineKm<-melt(finalFracListKm,id.vars=c("t","best","lower","upper"))
  combinedDataKm<-melt(finalDataListKm,id.vars=c("t","best"))
  
  combinedLineKp<-melt(finalFracListKp,id.vars=c("t","best","lower","upper"))
  combinedDataKp<-melt(finalDataListKp,id.vars=c("t","best"))
  
  combinedLineKm[,1]<-combinedLineKm[,1]*7
  combinedDataKm[,1]<-combinedDataKm[,1]*7
  
  combinedLineKp[,1]<-combinedLineKp[,1]*7
  combinedDataKp[,1]<-combinedDataKp[,1]*7
 
 
  combinedPlotKm<-ggplot(combinedLineKm,aes(t,best,colour=L1))+geom_line()+geom_ribbon(aes(ymin=lower,ymax=upper,colour=L1),alpha=0.2)+ylab("Brdu+ frac in Ki67-")+ xlab("t (days)")+ ylim(0,1)+geom_point(data=combinedDataKm)+global.theme
 
  combinedPlotKp<-ggplot(combinedLineKp,aes(t,best,colour=L1))+geom_line()+geom_ribbon(aes(ymin=lower,ymax=upper,colour=L1),alpha=0.2)+ylab("Brdu+ frac in Ki67+")+ xlab("t (days)")+ ylim(0,1)+geom_point(data=combinedDataKp)+global.theme+ggtitle(paste(modelLabel," B",toString(bBoxMax),"K",toString(kBoxMax),sep=""))
 
 
  write.xlsx(combinedLineKm,paste(fileE,"Ki67MinusBoot.xlsx",sep=""))
  write.xlsx(combinedLineKp,paste(fileE,"Ki67PlusBoot.xlsx",sep=""))
  ###################################################
  #make image
  
  #export image

  pdf(file=paste(fileE,"BootPlots.pdf",sep=""))
  grid.arrange(combinedPlotKp,combinedPlotKm, ncol = 1)  
  dev.off()


  return()
}
