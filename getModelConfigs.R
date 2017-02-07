################################################################
# This file contains the search range settings for the best 
# fit routine.

# Each cell and heterogeneity type has its own min and max
# box ranges and input source size settings

# Source sizes are in units of % of memory pool PER WEEK
# Currently allowed models: kinHet, kinHetExtended, kinHetExtended2
# tempHet

# Refer to MakeRFunctions.nb for exact definitions.
################################################################

################################################################
# USER MODIFICATION AREA BEGINS HERE
################################################################
if(cell=="4Tem"&&(model=="kinHet")){
  kMin<-14#7
  bMin<-1
  
  kB<-16
  bB<-2#3
	
  sourceList<-c(0,5,7,10)/100
}

################################################################
if(cell=="4Tem"&&(model=="kinHetExtended"||model=="kinHetExtended2")){
  kMin<-13
  bMin<-2
  
  kB<-17
  bB<-2
  
  sourceList<-c(0,5,7,10,15,20,25,30)/100
}


################################################################

if(cell=="4Tcm"&&(model=="kinHet")){
   kMin<-13
   bMin<-1
   
   kB<-17
   bB<-2#3
   
   sourceList<-c(0,5,10.6,10,15,20,25,30)/100
  
}


if(cell=="4Tcm"&&(model=="kinHetExtended"||model=="kinHetExtended2")){
  kMin<-13
  bMin<-2
  
  kB<-17
  bB<-2
  
  sourceList<-c(0,5,10.6,10,15,20,25,30)/100
  
}

################################################################

if(cell=="4Tem"&&model=="tempHet"){
  kMin<-1
  bMin<-1	
  
  kB<-5
  bB<-2
  sourceList<-c(0,5,7,10,15,20,25,30)/100
}

################################################################

if(cell=="4Tcm"&&model=="tempHet"){
  kMin<-1
  bMin<-1	
  
  kB<-5
  bB<-2

  sourceList<-c(0,5,10.6,10,15,20,25,30)/100
}

################################################################
# USER MODIFICATION AREA ENDS HERE
################################################################

################################################################
#used by main.R and bootstrap.R codes (basically actual simulation codes)
bCount<-1
modelConfigs<-matrix(-1, nrow =(kB-kMin+1)*(bB-bMin+1)*length(sourceList), ncol = 3)
for(sL in 1:length(sourceList)){
  for(kBM in kMin:kB){
    for(bBM in bMin:bB){
      modelConfigs[bCount,1] <- kBM
      modelConfigs[bCount,2] <- bBM
      modelConfigs[bCount,3] <- sourceList[sL]
      bCount<-bCount+1
    }
  }
}
################################################################
#used by analysis.R
bCount<-1
boxConfigs<-matrix(-1, nrow =(kB-kMin+1)*(bB-bMin+1), ncol = 2)
for(kBM in kMin:kB){
  for(bBM in bMin:bB){
    boxConfigs[bCount,1] <- kBM
    boxConfigs[bCount,2] <- bBM
    bCount<-bCount+1
  }
}
