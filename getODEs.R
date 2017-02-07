
########################################################################################################################################################################################################################
# Import equations, functions being used, initial conditions, and the two parameters 
# knocked out by forcing the ki67+/- ratio to be constant overall for the whole experiment.
# warnings turned off here because of annoying complaining to do with lines not being terminated properly in csv's or some shit
#edit oct 19 removed equations up/down import, now in getODEsolver.R for c++ function.

# Note: Now that ODEs are solved in C++ form, the actual equations are imported in getODEsolver.R instead of here, which just imports the names of the functions (boxes)
# and the expressions for the initial conditions (in .csv not .cpp form)

#import function names
options(warn=-1)
functions<-read.csv(paste(fileI,"fns.csv",sep=""), header = FALSE,quote = "", stringsAsFactors=F)[[1]]
options(warn=0)
# NOTE: don't want to parse here: fns just represents a vector of strings which are the names of the boxes we're solving for.
# these are subsequently used to call certain columns in dataframes, thus they must remain a string and *not* an expression, like the 'IC' or 'condition' things imported below.

#import initial conditions for upslope
options(warn=-1)
IC<-read.csv(paste(fileI,"IC.csv",sep=""), header = FALSE,quote = "", stringsAsFactors=F)[[1]]
options(warn=0)
IC<-parse(text=IC)

#import pair of equations which enforce Ki+ and Ki- populations being constant over time. These are solved for alpha and delta in each population
options(warn=-1)  
conditions<-read.csv(paste(fileI,"conditions.csv",sep=""), header = FALSE,sep = "\n", stringsAsFactors=F)[[1]]
options(warn=0)
conditions<-parse(text=conditions)


########################################################################################################################################################################################################################
# organize intermediate boxes into appropriate global function groupings
bmkmF<-c("BmKm")
bmkpF<-c()
bpkmF<-c()
bpkpF<-c()
for(i in 1:length(functions)){
  
  if(substr(functions[i],1,4)=="BmKp"){
    bmkpF<-c(bmkpF,functions[i])
    #bmkpF<-c(bmkpF,i)
  } else if(substr(functions[i],1,4)=="BpKm"){
    bpkmF<-c(bpkmF,functions[i])
    #bpkmF<-c(bpkmF,i)
  } else if(substr(functions[i],1,4)=="BpKp"){
    bpkpF<-c(bpkpF,functions[i])
    #bpkpF<-c(bpkpF,i)
  }
}

functionsChar<-functions #save string versions of functions for labelling purposes (used in ndSolve function in getFunctions package)
functions<-parse(text=functions)

fracTol<-10^-5
########################################################################################################################################################################################################################
#to avoid repeatedly calling if statements, we just do a single if-statement to define how to add the populations
#together for a given box configuration - essentially this is a problem because of lists of length 1 being jerks for the addColumns function (I think? can come back and fix probably)
if(kBoxMax>1&&bBoxMax>1){
  
  computeFinalFracs<-function(population){
    bmkm<-population[,bmkmF]
    bmkp<-addColumns(population[,bmkpF])
    bpkm<-addColumns(population[,bpkmF])
    bpkp<-addColumns(population[,bpkpF])
    
    kMinusFrac<-bpkm/(bmkm+bpkm)
    kPosFrac<-bpkp/(bmkp+bpkp)
    
    kMinusFrac[abs(kMinusFrac)<fracTol]<-0.
    kPosFrac[abs(kPosFrac)<fracTol]<-0.
    
    kMinusFrac[kMinusFrac>1]<-1.
    kPosFrac[kPosFrac>1]<-1.
    
    result<-cbind(kMinusFrac,kPosFrac)
    
#    if(any(is.nan(Tform(unlist(result))))) {browser()}
    return(result)
  }
  
  ##################################################
} else if(kBoxMax==1&&bBoxMax>1){
  
  computeFinalFracs<-function(population){
    bmkm<-population[,bmkmF]
    bmkp<-population[,bmkpF]
    bpkm<-addColumns(population[,bpkmF])
    bpkp<-addColumns(population[,bpkpF])
    
    kMinusFrac<-bpkm/(bmkm+bpkm)
    kPosFrac<-bpkp/(bmkp+bpkp)
    
    kMinusFrac[abs(kMinusFrac)<fracTol]<-0.
    kPosFrac[abs(kPosFrac)<fracTol]<-0.
    
    kMinusFrac[kMinusFrac>1]<-1.
    kPosFrac[kPosFrac>1]<-1.   
    
    result<-cbind(kMinusFrac,kPosFrac)
    
 #   if(any(is.nan(Tform(unlist(result))))) {browser()}
    return(result)
  } 
  
  ##################################################
} else if(kBoxMax>1&&bBoxMax==1) {
  computeFinalFracs<-function(population){
    bmkm<-population[,bmkmF]
    bmkp<-addColumns(population[,bmkpF])
    bpkm<-population[,bpkmF]
    bpkp<-addColumns(population[,bpkpF])
    
    kMinusFrac<-bpkm/(bmkm+bpkm)
    kPosFrac<-bpkp/(bmkp+bpkp)
    
    kMinusFrac[abs(kMinusFrac)<fracTol]<-0.
    kPosFrac[abs(kPosFrac)<fracTol]<-0.
    
    kMinusFrac[kMinusFrac>1]<-1.
    kPosFrac[kPosFrac>1]<-1.
      
    result<-cbind(kMinusFrac,kPosFrac)
    
#    if(any(is.nan(Tform(unlist(result))))) {browser()}
    return(result)
  }
  ##################################################
} else if(kBoxMax==1&&bBoxMax==1){
  
  computeFinalFracs<-function(population){
    bmkm<-population[,bmkmF]
    bmkp<-population[,bmkpF]
    bpkm<-population[,bpkmF]
    bpkp<-population[,bpkpF]
    
    kMinusFrac<-bpkm/(bmkm+bpkm)
    kPosFrac<-bpkp/(bmkp+bpkp)
    
    kMinusFrac[abs(kMinusFrac)<fracTol]<-0.
    kPosFrac[abs(kPosFrac)<fracTol]<-0.
    
    kMinusFrac[kMinusFrac>1]<-1.
    kPosFrac[kPosFrac>1]<-1.
       
    result<-cbind(kMinusFrac,kPosFrac)
    
 #   if(any(is.nan(Tform(unlist(result))))) {browser()}
    return(result)
  }
  
}
rm(i)
print('getODEs internal')
########################################################################################################################################################################################################################
#ENDFILE
########################################################################################################################################################################################################################
########################################################################################################################################################################################################################
