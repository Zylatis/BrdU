############################################################################################################
# getODEsolver.R
# contains code to import C++-ready ODEs, defines the ndsolve() function as well as the 
# getFunctionSingle/Doublepop() functions.
# IMPORTANT NOTE: this code contains more than just function definitions - it actually compiles the C++ code
# when run therefore this CANNOT BE RUN IN ISOLATION as it will try import files and compile the solver
# without knowing which cell type, which box configuration, etc, etc.
############################################################################################################
# Packages
library(odeintr) #C++ ODE solver package
library(xlsx)
############################################################################################################

# get c++ files. These are prepared in a series of MM files (but can be done with string manipulation in R if one so desires...)
cFnFileTemplate<-paste("CFunctions/",source.switch,"/B",toString(bBoxMax),"K",toString(kBoxMax),model,sep="")
options(warn=-1)
odeSystemUp <- paste(readLines(paste(cFnFileTemplate,"CppUp.txt",sep="")), collapse=" ")
odeSystemDown <- paste(readLines(paste(cFnFileTemplate,"CppDown.txt",sep="")), collapse=" ")
options(warn=0)

# Check that these files are imported correctly. We have to print to a file here because all of this is generally inside a parallel loop in main.R, so printing to the terminal is problematic.
if(nchar(odeSystemDown)==0||nchar(odeSystemUp)==0){
  write.table("Error in getODEsolver.R - failed to import C++ ODE files correctly", file=paste("ERROR",modelPrint,".txt", sep=""), sep="\t", row.names=F)
}
############################################################################################################

# By default the C++ solver has no additional headers or definitions. Here we can tell it to get some libraries and define some custom functions we need.
cppHeaderCode<-"
#include <math.h>
"

cppDef<-"double HeavisideTheta(double x){
  double val=-1;
  if(x<0){
   val=0;
  } else if(x>=0){
   val=1;
  }
  return(val);
}"

########################################################################################################################################################################################################################

# Cumbersome if statement time!
# Currently I don't know how to make the 'set_params' function more modular, so i've just split tempHet and kinHet with if() statements.
# Fear not, these are only run once per model (not once per SSR-calculation) so it's not costly.
# Each type of het has a different set of parameters which is why they need to be compiled differently.
# The kinHet case will be annotatd but the tempHet case is essentially identical but with different parameters

if(model=="kinHet"||model=="kinHetExtended"||model=="kinHetExtended2"){

  # The C++ solver is just a C++ function that takes in some parameter values and initial conditions and outputs the function values.
  # Thus, when we compile it we have to tell the compiler which internal bits are to be part of the function call i.e. what are parameters
  # whose values we will feed in. This is done in internalHeader
  internalHeader<-c("b","eps","sig","phi","N0","nu","a","del") 
  
  # This simply calls the compiler to make the C++ object files (stored ...somewhere temporarily). Specifying the sys_dim is important, failure to do so will result in 
  # it trying to guess and for D>3 it usually fucks it up.
  # We compile the up and downslope solver separately.

  compile_sys("systemUp", odeSystemUp, internalHeader,sys_dim=length(functions))
  
  # Only the downslope requires the Heaviside definition - during the upslope there is not concept of 'source switch' stuff.
  compile_sys("systemDown", odeSystemDown, internalHeader,sys_dim=length(functions),headers=cppHeaderCode,globals=cppDef)
  
  # Define our own internal NDSolve() function in R using the C++ code. This is the funciton that is called repeatedly.
  # It takes in tau, which is a specific value of t (so if you want to do a calculation over a vector of times you need sapply)
  # and 'parameters' which MUST BE a dataframe. 
  # Finally, IC is a NUMERIC VECTOR of initial conditions. The order of the initial conditions must match the order of the functions and ODEs (which is done in MM)
  ndSolveNewUp <- function(tau,parameters,IC){
    
    # each compliled system is given it's own set of functions.
    # so SYSTEM_set_params() sets the parameters for the current instance to be solved
    systemUp_set_params(b = with(as.list(parameters),eval(b)),
                     eps = with(as.list(parameters),eval(eps)),
                     sig = with(as.list(parameters),eval(sig)),
                     phi = with(as.list(parameters),eval(phi)),
                     N0 = with(as.list(parameters),eval(N0)), 
                     nu = with(as.list(parameters),eval(nu)),
                     a = with(as.list(parameters),eval(a)),
                     del = with(as.list(parameters),eval(del))
                     )
    # SYSTEM_at() propogates the system from t = 0 where it uses the initial condition values in 'IC" to time tau, returning the results as a dataframe
    resultRaw<-systemUp_at(IC, tau)
    # drop time column from dataframe, not needed really.
    result<-resultRaw[,-1] 
    # We keep track of the boxes by assigning the appropraite headers to the solution dataframe. (functionsChar just contains a vector of strings of the funciton names i.e. c("BmKM","BmKp1"...))
    colnames(result)<-functionsChar
    return(result)
  }
    
  # This is essentially identical except that it corresponds to the down-slope ODEs only. epsilon is set to zero (downslope condition) in higher functions, but in principle could be set to zero here.
  # Furthermore, a meta-function could be made to make both up and downslope systems at once, but it's hard because of the way odeintr defines these weird functions based on the SYSTEM pragma...
  ndSolveNewDown <- function(tau,parameters,IC){

    systemDown_set_params(b = with(as.list(parameters),eval(b)),
                        eps = with(as.list(parameters),eval(eps)),
                        nu = with(as.list(parameters),eval(nu)),
                        phi = with(as.list(parameters),eval(phi)),
                        sig = with(as.list(parameters),eval(sig)),
                        N0 = with(as.list(parameters),eval(N0)), 
                        a = with(as.list(parameters),eval(a)),
                        del = with(as.list(parameters),eval(del))
    )
    resultRaw<-systemDown_at(IC, tau)
    result<-resultRaw[,-1] 
    colnames(result)<-functionsChar
    return(result)
  }
  

}

############################################################################################################

if(model=="tempHet"){
  #compile solver - takes time!
  internalHeader<-c("b","eps","kappa","sig","phi","N0","nu","a","del") #unsure if the ordering here must match ordering in set_params, for now assume it does so we hard code - only 2 models so no that bad
  compile_sys("systemUp", odeSystemUp, internalHeader,sys_dim=length(functions))
  compile_sys("systemDown", odeSystemDown, internalHeader,sys_dim=length(functions),headers=cppHeaderCode,globals=cppDef)
  
  ndSolveNewUp <- function(tau,parameters,IC){
    #set pars
    systemUp_set_params(b = with(as.list(parameters),eval(b)),
                        eps = with(as.list(parameters),eval(eps)),
                        kappa = with(as.list(parameters),eval(kappa)),
                        sig = with(as.list(parameters),eval(sig)),
                        phi = with(as.list(parameters),eval(phi)),
                        N0 = with(as.list(parameters),eval(N0)), 
                        nu = with(as.list(parameters),eval(nu)),
                        a = with(as.list(parameters),eval(a)),
                        del = with(as.list(parameters),eval(del))
    )
    resultRaw<-systemUp_at(IC, tau)
    result<-resultRaw[,-1] #drop time column from dataframe, not needed really.
    colnames(result)<-functionsChar
    return(result)
  }
  
  ndSolveNewDown <- function(tau,parameters,IC){
    #set pars
    systemDown_set_params(b = with(as.list(parameters),eval(b)),
                          eps = with(as.list(parameters),eval(eps)),
                          kappa = with(as.list(parameters),eval(kappa)),
                          nu = with(as.list(parameters),eval(nu)),
                          phi = with(as.list(parameters),eval(phi)),
                          sig = with(as.list(parameters),eval(sig)),
                          N0 = with(as.list(parameters),eval(N0)), 
                          a = with(as.list(parameters),eval(a)),
                          del = with(as.list(parameters),eval(del))
    )
    resultRaw<-systemDown_at(IC, tau)
    result<-resultRaw[,-1] #drop time column from dataframe, not needed really.
    colnames(result)<-functionsChar
    return(result)
  }
  
}


#################################################################################################################################
# Define functions to be called in other areas (i.e. no 'running stuff' below here)

# Function to calculate all up and down-slopes in the case of kinetic hetereogeneity
calcFunctionDoublePop<-function(parsA,parsB,tUp,tDownList,downStartPos){
  
  # Calculate numeric values of upslope initial conditions based on current parameter set. The functional form for the IC_A and IC_B are identical
  # but the input parameters are different.
  tempA<-rep(NA,length(IC))
  tempB<-tempA
  for(i in 1:length(IC)){
    tempA[[i]] <- with(as.list(parsA),eval(IC[[i]]))
    tempB[[i]] <- with(as.list(parsB),eval(IC[[i]]))
  }
  
  ICA<-tempA
  ICB<-tempB
  remove(tempA,tempB)
  names(ICA)=functions
  names(ICB)=functions 
  
  # WARNING: DISGUSTING CODING HACK
  # So...the problem here is that for some reason the odeintr solver is fine having repeated timepoints fed to it SO LONG AS they're not the first ones. 
  # i.e. if you feed it t = 0 twice (so that we get two output values to match to the upslope which has two data points at each time point) it will freak out and die
  # In principle, what i've done here is stupid but changes introduce errors so i've left it for now.
  # What should actually be done is only the unique time points are calculated and then the multiplicities (to compare to data and make residuals) are dealt with after
  # I had a go at this but it required declaring more nasty global variable stuff or feeding even more parameters to these functions so I left it
  
  # WARNING MKII: This assumes that the upslope has two mice per timepoint. If you feed in another upslope this needs to change or be made modular as described above!
  tUp[1]<-10^-5
  tUp[2]<-1.1*tUp[1]
  
  # solve ODEs with given parameter set and time vector
  populationUpA<-ndSolveNewUp(tUp,parsA,ICA)
  populationUpB<-ndSolveNewUp(tUp,parsB,ICB)
  
  # Calculate initial conditions for the downslope. Recall that the odeintr solver is an IVP solver so the downslopes start (as far as the solver is concerned) at 't=0' with the initial 
  # conditions corresponding to the values of the upslope at the start of the downslope.
  
  tempListA <- vector("list", nDown)
  tempListB<-tempListA
  for(i in 1:nDown){
    tempListA[[i]] <- populationUpA[downStartPos[i],]
    tempListB[[i]] <- populationUpB[downStartPos[i],]
  }
  
  downICAlist <-  do.call(rbind, tempListA)
  downICBlist <-  do.call(rbind, tempListB)
  
  # do sep for pop A and B so can use list and rbind once only - rbind in for loop a very bad thing!
  downPopulationA<-list()
  downPopulationB<-list()
  
  #  set brdu uptake to zero for downslope calculation
  epsPos<-match("eps",header)
  parsA[epsPos]<-0
  parsB[epsPos]<-0
  
  #temp - shift initial time points here to stop new ODE solver crapping itself. Doesn't like multiple points with the same exact value, especially if zero
  for(i in 1:nDown){
    tDownList[[i]][2]<-tDownList[[i]][2]+0.0001
    tDownList[[i]][3]<-tDownList[[i]][3]+0.0002
  }
  
  for(i in 1:nDown){
    taui<-tDownList[[i]]
    # have to reshift times for downslope here because by default the ode integrator starts at zero i.e. IVP not BVP, sorta.
    tau<-taui-min(taui)
    downPopulationA[[i]]<-ndSolveNewDown(tau,parsA,unlist(downICAlist[i,]))
    downPopulationB[[i]]<-ndSolveNewDown(tau,parsB,unlist(downICBlist[i,]))
  }  
  
  # combine the 'count' values (total N0=1 here so not actual counts) of the twosub-populations
  populationUp<-populationUpA+populationUpB
  pooledDownPopulations<-list()
  for(i in 1:nDown){
    pooledDownPopulations[[i]]<-downPopulationA[[i]]+downPopulationB[[i]]
  }
  #can only be one upslope curve so we duplicate this curve for the number of upslope datasets we have!
  output<-c(rep(list(populationUp),nUp),pooledDownPopulations)
  # Return vector of lists of dataframes, yay!
  return(output)
}

#################################################################################################################################
# Same as above but for single population (tempHet)
# Coding: todo - make single function to take in pars and output both up and downslope counts, then can be used by temp and kinHet?
calcFunctionSinglePop<-function(pars,tUp,tDownList,downStartPos){
  
  temp<-rep(NA,length(IC))
  #browser()
  for(i in 1:length(IC)){
    temp[[i]] <- with(as.list(pars),eval(IC[[i]]))
  }
  
  IC<-temp
  remove(temp)
  names(IC)=functions
  
  tUp[1]<-10^-5
  tUp[2]<-1.1*tUp[1]

  #solve ODEs with given parameter set and time vector
  populationUp<-ndSolveNewUp(tUp,pars,IC)#ndSolve(tauUp1,parsA,equationsUp,ICA)
  
  tempList <- vector("list", nDown)
  for(i in 1:nDown){
    tempList[[i]] <- populationUp[downStartPos[i],]
  }
  
  downIClist <-  do.call(rbind, tempList)

  downPopulation<-list()
  
  #  set brdu uptake to zero for downslope calculation
  epsPos<-match("eps",header)
  pars[epsPos]<-0

  #temp - shift initial time points here to stop new ODE solver crapping itself. Doesn't like multiple points with the same exact value, especially if zero
  for(i in 1:nDown){
    tDownList[[i]][2]<-tDownList[[i]][2]+0.0001
    tDownList[[i]][3]<-tDownList[[i]][3]+0.0002
  }
  
  for(i in 1:nDown){
    taui<-tDownList[[i]]
    # have to reshift times for downslope here because by default the ode integrator starts at zero i.e. IVP not BVP, sorta.
    tau<-taui-min(taui)
    downPopulation[[i]]<-ndSolveNewDown(tau,pars,unlist(downIClist[i,]))
  }  
  
  populationList<-NULL
  return(c(rep(list(populationUp),nUp),downPopulation))
}
print('getODESolver internal')
#################################################################################################################################
