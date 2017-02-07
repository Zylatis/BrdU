
require(stats)
library(nleqslv)
require(xlsx)
require(gdata)
library(deSolve)
library(fBasics)
library(GenSA)
library(ggplot2)
library(reshape)
library(gridExtra)
pointSize<-0.025 #pointsize for kpos and kminus frac dots
#q<-expression((phiGlobal - phiB)/(phiA - phiB))
startx <- c(1,1.1) #starting points for condition.solve() fn which uses the nleqslv package to solve for 'a' and 'delta' in a given model.

############################################################################################################

#function used to transform data when calcualting residuals. This is the same for both datasets as they're both bounded between 0 and 1
Tform <- function(x) {
  asin(sqrt(x))
}

calc.q<-function(phiG,phiA,phiB){
  (phiG-phiB)/(phiA-phiB)
}

sigmaB.sub<-function(zeta,sigmaG,phi,phiA,phiGlobal){
  (1-zeta)*sigmaG*(phiA-phi)/(phiA-phiGlobal)
}
############################################################################################################

# function used when summing over various K+ or B+ subcompartments to make final two functions
# faster way to do this with inbuilt functions? vectorized etc... 
addColumns<-function(vec){
  result<-rep(0,nrow(vec))
  for(i in 1:ncol(vec)){ 
    result<-vec[,i] + result
  }
  return(result)
}

############################################################################################################
# Mapping from Mathematica's output of certain functions
Sqrt<-function(x){ return(sqrt(x))}
HeavisideTheta<-function(x){ Heaviside(x) }

############################################################################################################
# function used to solve for a and delta
# the vector 'conditions' contains the LHS of the two equations which are set to = 0 and solved for 'a' and 'delta' of that [sub]population.
# solving these equations corresponds to imposing two of the 'constant box size' conditions (the rest are imposted by the definition of the initial conditions in 'IC')
condition.solve<-function(x,inputPars){ # NOTE: 26th Nov did replace-all to rename this conditionSolve->condition.solve so if other files throw errors this is the first port of call! Checked most seems okay...
  a<-x[1]
  del<-x[2]
  output<-numeric(length(x))
  output[1]<-with(as.list(inputPars),eval(conditions[[1]]))
  output[2]<-with(as.list(inputPars),eval(conditions[[2]]))
  output
}


############################################################################################################

#calculate residuals based on solved values of the BrdU+ pos frac in K+/- populations (pooled over whatever subpops exist)
calcResiduals<-function(fracs,upData,downData){

  outputList<-vector("list", 2*(nUp+nDown))
  residualSetUp <- vector("list", nUp)
  residCounter<-1
  for(i in 1:nUp){
    residualSetUp[[i]] <- Tform(upData[[i]])-Tform(fracs[[i]]) #residual = data-function!! important when seeing what to add back for bootstrap
     outputList[[residCounter]]<- residualSetUp[[i]][,1]
     outputList[[residCounter+1]]<- residualSetUp[[i]][,2]
    residCounter<-residCounter+2
  }
  
  residualSetDown <- vector("list", nDown)
  for(i in 1:nDown){
    residualSetDown[[i]]<-Tform(downData[[i]])-Tform(fracs[[i+nUp]]) #offset as 'fracs' includes upslope data too!
    outputList[[residCounter]]<- residualSetDown[[i]][,1]
    outputList[[residCounter+1]]<- residualSetDown[[i]][,2]
    residCounter<-residCounter+2
  }
  
  return(outputList)
}

############################################################################################################
# SORT OF HACK TO SET ALPHA_B = 0 IN KINHET CASE BASED ON TOGGLE IN MAIN.R
# toggle is zero.alphaB
############################################################################################################

  calcParsKinHet<-function(inVec){
    inVec<-as.data.frame(t(as.numeric(inVec)))
    names(inVec)<-header
    
    inVecList<-as.list(inVec)
    qVal<-with(inVecList,calc.q(phiGlobal,phiA,phiB))
    
    #parameters used by both sub-populations
    commonPars<-cbind(
      b=with(inVecList,b),
      eps=with(inVecList,eps),
      nu=with(inVecList,eval(nuV))
    )
    #sig_G*N_0*\zeta = sig_A*N_A => sig_A = sig_G*zeta/N_A
    # Here we deal with the scale factor N0 by setting the total population size = 1 then population A is q*N0=q, and B is (1-q). Thus when we combine the 'counts' from populations A and B
    # they are on the same 'scale' and can indeed be added to get a total size.
    parsA<-as.data.frame(cbind(commonPars,phi=with(inVecList,phiA),sig=with(inVecList,zeta*eval(sigV)/qVal), N0=qVal))
    parsB<-as.data.frame(cbind(commonPars,phi=with(inVecList,phiB),sig=with(inVecList,(1-zeta)*eval(sigV)/(1-qVal)), N0=1-qVal))
    
    knockedParsA<-(nleqslv(startx,condition.solve,jac=NULL,parsA,control=list(ftol=solTolerance,allowSingular=TRUE,maxit=120)))[[1]]  #allow singular jacobian try 12/9/15
    knockedParsB<-(nleqslv(startx,condition.solve,jac=NULL,parsB,control=list(ftol=solTolerance,allowSingular=TRUE,maxit=120)))[[1]] 
    
    #2Dec removed try catch on knockedPars 
    alphaNdeltaA <- cbind(a = knockedParsA[[1]], del = knockedParsA[[2]])
    alphaNdeltaB <- cbind(a = knockedParsB[[1]], del = knockedParsB[[2]])
    
    parsA<-cbind(parsA,alphaNdeltaA)
    parsB<-cbind(parsB,alphaNdeltaB)
    return(list(parsA,parsB))
  }
 
############################################################################################################

calcParsTempHet<-function(inVec){
  inVec<-as.data.frame(t(as.numeric(inVec)))
  names(inVec)<-header
  
  inVecList<-as.list(inVec)

  #sig_G*N_0*\zeta = sig_A*N_A => sig_A = sig_G*zeta/N_A
  pars<-data.frame(b=with(inVecList,b),eps=with(inVecList,eps),nu=with(inVecList,eval(nuV)),phi=phiGlobal,N0=1,sig=sigV,kappa=with(inVecList,kappa))
  
  knockedPars<-(nleqslv(startx,condition.solve,jac=NULL,pars,control=list(ftol=solTolerance,allowSingular=TRUE,maxit=120)))[[1]]
  
  alphaNdelta <- cbind(a = knockedPars[[1]], del = knockedPars[[2]])
  
  pars<-cbind(pars,alphaNdelta)

  return(pars)
}
############################################################################################################

computeSSRkinHet <- function(inVec,tUp,tDown,upDlist,downDlist){
  pars<-calcParsKinHet(inVec)
  parsA<-pars[[1]]
  parsB<-pars[[2]]
  #crude fix for avoiding b<bmin vals, seems strongly correlated with negative values in counts (source too large, box size constancy etc) but no counter examples seen yet
  #now importing imaginary and real parts of knocked out pars separately because R is a jerkface.
  if(parsA[["a"]]<0||parsA[["del"]]<0||parsB[["a"]]<0||parsA[["del"]]<0||parsA["a"]>10||parsA["del"]>10||parsB["a"]>10||parsB["del"]>10||parsB["sig"]<0||parsB["phi"]>phiGlobal){
    SSR<-10^10
  } else {    
    functionVals<-calcFunctionDoublePop(parsA,parsB,tUp,tDown,downStartUpPos)

    fracs<-lapply(functionVals,computeFinalFracs)

    fracsSSR<-fracs
    fracs[[1]][1:2,]<-0                         #force t = 0 upslope fracs to be zero, avoids annoying small negative numbers that can only be removed with costly if statements
    
    residualList<-calcResiduals(fracs,upDlist,downDlist) # feed in actual experimental data, this changes when we do bootstraps!

    #nrows now defined in getData.R
    SSR<-0
    check<-NULL
  #  browser()
    for(i in ssrStart:length(residualList)){ #only include SSR contributions from certain sets
      SSR<- SSR+ nrows[i]*(log(sum(residualList[[i]]^2)))
      check[i]<-nrows[i]*(log(sum(residualList[[i]]^2)))
    }
  }
  if(is.nan(SSR)){
    SSR = 10^10
  }
#browser()
  return(SSR)
}


############################################################################################################

computeSSRtempHet <- function(inVec,tUp,tDown,upDlist,downDlist){
  pars<-calcParsTempHet(inVec)
  #crude fix for avoiding b<bmin vals, seems strongly correlated with negative values in counts (source too large, box size constancy etc) but no counter examples seen yet
  #now importing imaginary and real parts of knocked out pars separately because R is a jerkface.
  if(pars[["a"]]<0||pars[["del"]]<0||pars["a"]>10||pars["del"]>10){
    SSR<-10^10
  } else {   

    functionVals<-calcFunctionSinglePop(pars,tUp,tDown,downStartUpPos)

    fracs<-lapply(functionVals,computeFinalFracs)
    fracsSSR<-fracs
    fracs[[1]][1:2,]<-0                         #force t = 0 upslope fracs to be zero, avoids annoying small negative numbers that can only be removed with costly if statements
    
    residualList<-calcResiduals(fracs,upDlist,downDlist) # feed in actual experimental data, this changes when we do bootstraps!
    
    #nrows now defined in getData.R
    SSR<-0
    check<-c()
    for(i in ssrStart:length(residualList)){ #only include SSR contributions from certain sets
      SSR<- SSR + nrows[i]*(log(sum(residualList[[i]]^2)))
    }
    if(is.nan(SSR)){
      SSR = 10^10
    }
  }
  
  return(SSR)
}

############################################################################################################


makePlot <- function(inVec,plotDataUpList,plotDataDownList){
  inVec<-as.data.frame(t(as.numeric(inVec)))
  names(inVec)<-header
  
  inVecList<-as.list(inVec)
  
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
  if(model=="kinHet"||model=="kinHetExtended"||model=="kinHetExtended2"){
        
    pars<-calcParsKinHet(inVec)
    parsA<-pars[[1]]
    parsB<-pars[[2]]
    
    functionVals<-calcFunctionDoublePop(parsA,parsB,tauUpLine,tauDownLineList,downStartUpPosLine)
    
  }
  
  ############################################################################################################
  if(model=="tempHet"){
                      
                      pars<-calcParsTempHet(inVec)
                      
                      functionVals<-calcFunctionSinglePop(pars,tauUpLine,tauDownLineList,downStartUpPosLine)
                      
  }
  ############################################################################################################
  
  fracs<-lapply(functionVals,computeFinalFracs)

  finalFracListKm <- list()
  finalDataListKm <- list()
  finalFracListKp <- list()
  finalDataListKp <- list()
  
  upFracs<-fracs[seq(1:nUp)] 
  for(i in 1:nUp){
    finalFracListKm[[i]]<-data.frame(t=tauUpLine,kMinusFrac=upFracs[[i]][,1]) #super inefficient to keep calling computeFinalFracs, but okay for now as a hack
    finalDataListKm[[i]]<-data.frame(t=tauUp1,kMinusFrac=plotDataUpList[[i]][,1])
    
    finalFracListKp[[i]]<-data.frame(t=tauUpLine,kPosFrac=upFracs[[i]][,2])
    finalDataListKp[[i]]<-data.frame(t=tauUp1,kPosFrac=plotDataUpList[[i]][,2])
  }

  downFracs<-fracs[-seq(1:nUp)] 
  for(i in 1:nDown){
    finalFracListKm[[i+nUp]]<-data.frame(t=tauDownLineList[[i]],kMinusFrac=downFracs[[i]][,1]) #super inefficient to keep calling computeFinalFracs, but okay for now as a hack
    finalDataListKm[[i+nUp]]<-data.frame(t=tauDownList[[i]],kMinusFrac=plotDataDownList[[i]][,1])
    
    finalFracListKp[[i+nUp]]<-data.frame(t=tauDownLineList[[i]],kPosFrac=downFracs[[i]][,2])
    finalDataListKp[[i+nUp]]<-data.frame(t=tauDownList[[i]],kPosFrac=plotDataDownList[[i]][,2])
  }

  for(i in 1:length(finalFracListKm)){
    colnames( finalFracListKm[[i]])<-c("t",labelsKm[i])
    colnames( finalDataListKm[[i]])<-c("t",labelsKm[i])
  }

  for(i in 1:length(finalFracListKp)){
    colnames( finalFracListKp[[i]])<-c("t",labelsKp[i])
    colnames( finalDataListKp[[i]])<-c("t",labelsKp[i])
  }

  names(finalFracListKm)<-labelsKm
  names(finalDataListKm)<-labelsKm
  
  names(finalFracListKp)<-labelsKp
  names(finalDataListKp)<-labelsKp

  ############################################################################################################
  cols <- c("blue", "red","green")
  global.theme<-theme(axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=12)) +theme(plot.title = element_text(family = "Helvetica", color="#666666", face="bold", size=15, hjust=0)) +theme(legend.position="none")
  
  
  combinedLineKm<-melt(finalFracListKm,id.vars="t")
  combinedDataKm<-melt(finalDataListKm,id.vars="t")

  combinedLineKp<-melt(finalFracListKp,id.vars="t")
  combinedDataKp<-melt(finalDataListKp,id.vars="t")
  
  combinedLineKm[,1]<-combinedLineKm[,1]*7
  combinedDataKm[,1]<-combinedDataKm[,1]*7

  combinedLineKp[,1]<-combinedLineKp[,1]*7
  combinedDataKp[,1]<-combinedDataKp[,1]*7

  combinedPlotKm<-ggplot(combinedLineKm,aes(t,value,colour=L1))+geom_line()+ylab("Brdu+ frac in Ki67-")+ xlab("t (days)")+ ylim(0,1)+geom_point(data=combinedDataKm)+global.theme
  combinedPlotKp<-ggplot(combinedLineKp,aes(t,value,colour=L1))+geom_line()+ylab("Brdu+ frac in Ki67+")+ xlab("t (days)")+ ylim(0,1)+geom_point(data=combinedDataKp)+global.theme+ggtitle(paste(modelLabel," B",toString(bBoxMax),"K",toString(kBoxMax),sep=""))

  #look at images on TRANSFORMED scaled using Tform(x)

  combinedLineKm.transformed<-combinedLineKm
  combinedLineKp.transformed<-combinedLineKp
  
  combinedDataKm.transformed<-combinedDataKm
  combinedDataKp.transformed<-combinedDataKp
  
  combinedLineKm.transformed[,3]<-Tform(combinedLineKm.transformed[,3])
  combinedLineKp.transformed[,3]<-Tform(combinedLineKp.transformed[,3])
  
  combinedDataKm.transformed[,3]<-Tform(combinedDataKm.transformed[,3])
  combinedDataKp.transformed[,3]<-Tform(combinedDataKp.transformed[,3])
    
  combinedPlotKm.transformed<-ggplot(combinedLineKm.transformed,aes(t,value,colour=L1))+geom_line()+ylab("Brdu+ frac in Ki67-")+ xlab("t (days)")+ ylim(0,1)+geom_point(data=combinedDataKm.transformed)+global.theme
  combinedPlotKp.transformed<-ggplot(combinedLineKp.transformed,aes(t,value,colour=L1))+geom_line()+ylab("Brdu+ frac in Ki67+")+ xlab("t (days)")+ ylim(0,1.5)+geom_point(data=combinedDataKp.transformed)+global.theme+ggtitle(paste(modelLabel," B",toString(bBoxMax),"K",toString(kBoxMax),sep=""))
 
  
###################################################
  #make image
  
  #finalPlot<-grid.arrange(combinedPlotKm, combinedPlotKp, ncol = 1, top=paste("B",toString(bBoxMax),"K",toString(kBoxMax)," ",model,sep=""))
  #export image
  pdf(file=paste(fileE,"Plots.pdf",sep=""), useDingbats=FALSE)
  grid.arrange(combinedPlotKp,combinedPlotKm, ncol = 1) 
  dev.off()

  pdf(file=paste(fileE,"TransformedPlots.pdf",sep=""), useDingbats=FALSE)
  grid.arrange(combinedPlotKp.transformed,combinedPlotKm.transformed, ncol = 1) 
  dev.off()
  return()
}

############################################################################################################

findFit<-function(tUp,tDown,upDlist,downDlist){  
#  traceFile<-paste(fileE,"TraceMatrix.txt",sep="")
#  trace.fn=NULL #traceFile - changed on 20/01/16 so that bootstraps stop saving trace files which are ~4mb each. Ideally we want to have this only off for bootstrap and on for best fit but same function call for all...

  out <- GenSA(par=start,lower = lower, upper = upper, fn = computeSSR, control=list(verbose=TRUE,temperature=temp,maxit=maxit),tUp,tDown,upDlist,downDlist)
  fitPars<-out["par"][[1]]
  return(fitPars)
}


#################################################################################################################################

#handy debug stuff

#################################################################################################################################
testConstancy<-function(inVec,tUp,tDown){
  pars<-calcParsKinHet(inVec)
  parsA<-pars[[1]]
  parsB<-pars[[2]]

  allPops<-calcFunctionDoublePop(parsA,parsB,tUp,tDown,downStartUpPos)
  
  test<-toString(functions)
  allFns<-unlist(strsplit(test,split=", "))
  KpList<-c(grep("BmKp+",allFns, perl=TRUE,value=TRUE),grep("BpKp+",allFns, perl=TRUE,value=TRUE))  
  KmList<-c(grep("BmKm+",allFns, perl=TRUE,value=TRUE),grep("BpKm+",allFns, perl=TRUE,value=TRUE))  

  kNeg<-c(paste("BmKm",sep=""), foreach(bp = 1:bBoxMax,.combine='c') %do% {
    paste("BpKm",toString(bp),sep="")
  }
  )
  print(kNeg)
  outkNeg<-foreach(i = 1:length(allPops)) %do% {
    val<-rowSums(allPops[[i]][kNeg])
    val<-as.data.frame(val)
    names(val)<-paste("Summed(Brdu) Km",sep="")
    val
    }

  #check single K+ stage summed over Brdu to ensure constant over time
  all.compartments<-foreach(kp = 1:kBoxMax,.combine='rbind') %do% {
    allKp<-c(paste("BmKp",toString(kp),sep=""), foreach(bp = 1:bBoxMax,.combine='c') %do% {
      paste("BpKp",toString(bp),toString(kp),sep="")
    }
    )
    allKp
  }

  labels<-foreach(kp = 1:kBoxMax,.combine='rbind') %do% {
   paste("Summed(Brdu) Kp",toString(kp),sep="")
  }
  
  out<-foreach(i = 1:length(allPops)) %do% {
    const <- foreach(comp = 1:nrow(all.compartments),.combine='cbind') %do% {
      rowSums(allPops[[i]][all.compartments[comp,]])
    }
    const<-as.data.frame(const)
    names(const)<-labels[,1]
    const
  }
  #names(out)<-all.compartments
final<-  foreach(i=1:length(out)) %do% {
    cbind(outkNeg[[i]],out[[i]])
  }
  #print(out[[1]])
  #print(outkNeg[[1]])
  print(final)

  
}
#################################################################################################################################
test.nleqslv<-function(inVec){
  inVec<-as.data.frame(t(as.numeric(inVec)))
  names(inVec)<-header
  
  inVecList<-as.list(inVec)
  qVal<-with(inVecList,calc.q(phiGlobal,phiA,phiB))
  #sig_G*N_0*\zeta = sig_A*N_A => sig_A = sig_G*zeta/N_A
  commonPars<-cbind(b=with(inVecList,b),eps=with(inVecList,eps),nu=with(inVecList,eval(nuV)))
  
  parsA<-as.data.frame(cbind(commonPars,phi=with(inVecList,phiA),sig=with(inVecList,zeta*eval(sigV)/qVal), N0=qVal))
  parsB<-as.data.frame(cbind(commonPars,phi=with(inVecList,phiB),sig=with(inVecList,(1-zeta)*eval(sigV)/(1-qVal)), N0=1-qVal))
  tryCatch({
    knockedParsA<-(nleqslv(startx,condition.solve,jac=NULL,parsA,control=list(ftol=solTolerance,allowSingular=TRUE,maxit=120)))  #allow singular jacobian try 12/9/15
    knockedParsB<-(nleqslv(startx,condition.solve,jac=NULL,parsB,control=list(ftol=solTolerance,allowSingular=TRUE,maxit=120))) 
  }, 
  error = function(c) { browser()
  }
  )
  return(list(knockedParsA,knockedParsB))
}

# doFit<-c(2.0444666903,0.7656403472,0.6863492737,0.0901125212,0.2258493704,0.2430867427)
# testConstancy(doFit,tauUp1,tauDownList)

#################################################################################################################################

calcParsKinHetExtended<-calcParsKinHet
computeSSRkinHetExtended<-computeSSRkinHet

calcParsKinHetExtended2<-calcParsKinHet
computeSSRkinHetExtended2<-computeSSRkinHet

#################################################################################################################################
