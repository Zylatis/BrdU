#############################################################################################
# Set of functions used in the analysis of data from fits
# It is here that we define parameters that are not directly extracted from the fits, i.e. 
# population weighted parameters
############################################################################################# 

weightedSum<-function(qV,valA,valB){
  qV*valA+(1-qV)*valB
}


calc.q<-function(phiG,phiA,phiB){
  (phiG-phiB)/(phiA-phiB)
}

paste.vectors<-function(name.vector){
  c(paste(name.vector,".lower",sep=""),paste(name.vector,".upper",sep=""))
}

global.taus<-function(qVal,aVal,bVal){
  qVal/aVal+(1-qVal)/bVal
}
############################################################################################# 
global.pars<-c("alphaG","deltaG","ki67.lifetime","global.interdivisiontime","global.lifetime","source.frac",'nu','invNu')

if(model=="kinHet"||model=='kinHetExtended'||model=="kinHetExtended2"){
  other.pars<-c("q","interdivision.time.A","interdivision.time.B","lifetime.A","lifetime.B")
  implicit.plot.pars<-c("zeta",'eps', 'phiA','phiB') #these are parameters that are direct outputs of the fitting procedure and that we want to plot (see combineBootData.R for where this is added in)
}

if(model=="tempHet"){
  other.pars<-c("lossrate.ratio","invDelta","invKappa")
  implicit.plot.pars<- c('eps')
}


get.parameters<-function(header,output,model,sigVal){
  names(output)<-header
  phiGlobalpos<-match("phiGlobal",header)
  phiGlobal<-output[,phiGlobalpos]
  ki67.lifetime<-1/output["b"][,1]
  header<-c(header,global.pars,other.pars)
  eps<- output['eps'][,1]
  ############################################################################################ # kinHet/tempHet additional calcs if statement
  if(model=="kinHet"||model=='kinHetExtended'||model=='kinHetExtended2'){


    
    phiApos<-match("phiA",header)
    phiBpos<-match("phiB",header)
    alphaApos<-match("alphaA",header)
    alphaBpos<-match("alphaB",header)
    deltaApos<-match("delA",header)
    deltaBpos<-match("delB",header)
    epsPos<-match("eps",header)
    #Calculate and add in q
    qList<-calc.q(phiGlobal,output[,phiApos],output[,phiBpos])
    
    alphaA<-output["alphaA"][,1]
    alphaB<-output["alphaB"][,1]

    phiA<-output['phiA'][,1]
    phiB<-output['phiB'][,1]
    

    if(source.switch == 'immediate'){
      nuVal<-0
      invNu<-0
    } else {
      nuVal<-output['nu'][,1]
      invNu <- 1/nuVal
    }
    alphaG<-weightedSum(qList,alphaA,alphaB)
    source.contrib<-sigVal/(sigVal+alphaG)
    global.interdivisiontime<-global.taus(qList,alphaA,alphaB)
    
    if(model=='kinHet'){
      deltafrac <- 1
    }
    if(model=='kinHetExtended'){
      deltafrac <- 10
    }
    if(model=='kinHetExtended2'){
      deltafrac <- 0.1
    }
    
    deltaAki67Neg<-output["delA"][,1]
    deltaAki67Pos<-deltafrac*deltaAki67Neg
    lifetimeA <- global.taus(phiA, deltaAki67Pos, deltaAki67Neg)
    
    deltaBki67Neg<-output["delB"][,1]
    deltaBki67Pos<-deltafrac*deltaBki67Neg
    lifetimeB <- global.taus(phiB, deltaBki67Pos, deltaBki67Neg)
    
    global.lifetime<-weightedSum(qList,lifetimeA,lifetimeB)
    
    
    deltaA<- weightedSum(phiA, deltaAki67Pos, deltaAki67Neg)
    deltaB<- weightedSum(phiB, deltaBki67Pos, deltaBki67Neg)
    deltaG<-weightedSum(qList,deltaA,deltaB)

    #transformation from days to weeks not done here, only in Analysis.R; thus everything up to Analysis.R assumes weeks, hence the 1/alpha not 7/alpha
    output<-cbind(output,alphaG,deltaG,ki67.lifetime,global.interdivisiontime,global.lifetime,source.contrib,nuVal,invNu,qList,1/alphaA,1/alphaB,lifetimeA,lifetimeB)
    names(output)<-header
  }
  
  

  if(model=="tempHet"){ 

    kappaPos<-match("kappa",header)
    deltaPos<-match("del",header)
    
    alphaPos<-match("alpha",header)        #re-name alpha->alphaG so in line with kinHet and we can compute sigG/alphaG in same code
    alphaG<-output[,alphaPos]

    if(source.switch == 'immediate'){
      nuVal<-0
      invNu<-0
    } else {
      nuVal<-output['nu'][,1]
      invNu <- 1/nuVal
    }
    deltaG<-phiGlobal*output[,kappaPos]+(1-phiGlobal)*output[,deltaPos] # delta_G = phi*kappa+(1-phi)*delta weighted sum of ki67+/- death rates
    invKappa = 1/output[,kappaPos]
    invDelta = 1/output[,deltaPos]
    source.contrib<-sigVal/(sigVal+alphaG)
    
    global.interdivisiontime<-1/output[,alphaPos]
    global.lifetime<-global.taus(phiGlobal,output[,kappaPos],output[,deltaPos]) #weighted average of lifetimes based on ki67+% 
    rel.loss.rate<-output[,deltaPos]/output[,kappaPos]
    
    output<-cbind(output,alphaG,deltaG,ki67.lifetime,global.interdivisiontime,global.lifetime,source.contrib,nuVal,invNu,rel.loss.rate,invDelta,invKappa)
    names(output)<-header
  }
  ############################################################################################ # kinHet/tempHet additional calcs if statement


  return(list(output,header))
}
