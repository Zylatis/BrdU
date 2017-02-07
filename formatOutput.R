######################################################################################
# Set of functions to format the output of the fits to be exported in xlsx format
# and used by the analysis files
######################################################################################



#need to add back in parameters that are set to constants, i.e. ones that GenSA doesn't see and thus don't emerge in doFit output
######################################################################################
if(eval(sigV) != -2&& !("sig" %in% header)){ #don't want to keep adding if re-running small segments in de-bug mode!
  doFit<-c(doFit,sig=eval(sigV))
  header<-c(header,"sig")
}

if((source.switch=="delayStep"||source.switch=="expSwitch")&& nuVstat =="fixed"){ 
  doFit<-c(doFit,nu=eval(nuV))
  header<-c(header,"nu")
}

######################################################################################

kinHetFormat<-function(uD,dD,doPlot){  
  
  ######################################################################################
  #back-calculate alpha and delta and format results for output
  ######################################################################################
  doFitList<-as.list(doFit)
  
  #prbably should output q here but anyway...done in analysis
  #qVal<-with(doFitList,calc.q(phiGlobal,phiA,phiB))

  allPars<-calcParsKinHet(doFit)
  
  setA<-allPars[[1]][c("a","del")]
  setB<-allPars[[2]][c("a","del")]

  output<-c(computeSSR(doFit,tauUp1,tauDownList,uD,dD),doFit,phiGlobal,setA,setB)
  names(output)<-c("lnL",header,"phiGlobal","alphaA","delA","alphaB","delB")

  if(sigV==0){
    output<-c(output,0)
    names(output)<-c("lnL",header,"phiGlobal","alphaA","delA","alphaB","delB","zeta")
  }
  
  # if we have not fitted phiB (a_B=0) then it will not appear in doFit so we have to restore in manually.
  if(zero.alphaB==1){
    output<-c(computeSSR(doFit,tauUp1,tauDownList,uD,dD),doFit,allPars[[2]]["phi"],phiGlobal,setA,setB)
    names(output)<-c("lnL",header,"phiB","phiGlobal","alphaA","delA","alphaB","delB")
  }
  
  write.xlsx(cbind(output),paste(fileE,"bestfit.xlsx",sep=""))
  
  write.xlsx(toggles,paste(fileE,"SAtoggles.xlsx",sep=""))
  #write.table(startx,file=paste(fileE,"startx.txt",sep=""))
  if(doPlot==1){
    makePlot(doFit,uD,dD) #assuming only one upslope (can always add new data to same curve - only parameterized by 'time on' time)
  }

}
kinHetExtendedFormat<-kinHetFormat
kinHetExtended2Format<-kinHetFormat
######################################################################################

tempHetFormat<-function(uD,dD,doPlot){
  
  ######################################################################################
  #back-calculate alpha and delta and format results for output
  ######################################################################################
  doFitList<-as.list(doFit)
  
  allPars<-calcParsTempHet(doFit)
  knockedPars<-allPars[c("a","del")]

  output<-c(computeSSR(doFit,tauUp1,tauDownList,uD,dD),doFit,phiGlobal,knockedPars)
  names(output)<-c("lnL",header,"phiGlobal","alpha","del")
  write.xlsx(cbind(output),paste(fileE,"bestfit.xlsx",sep=""))
  
  write.xlsx(toggles,paste(fileE,"SAtoggles.xlsx",sep=""))
  #write.table(startx,file=paste(fileE,"startx.txt",sep=""))
  if(doPlot==1){
    makePlot(doFit,uD,dD) #assuming only one upslope (can always add new data to same curve - only parameterized by 'time on' time)
  }
}
