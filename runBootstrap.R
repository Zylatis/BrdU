rm(list = ls())
############################################################################################

# Code to perform N # of bootstrap replicates using sample with replacement of residuals

# Requires best fit to exist already and have had Analysis.R run to compute the optimal
# box configuration for each source term, which is then used as the best fit around which
# the replicates are constructed.

# Outputs raw bootstrap replicate data, i.e. N # of best fit parameters (in a single file)

############################################################################################
# BEGIN USER EDIT AREA

upSets<-c('Expt185','Expt203')
downSets<-c('Expt184','Expt196','Expt197')

# Number of bootstrap replicates to cover
iters<-1500

# END USER EDIT AREA
########################################################################################

args <- commandArgs(trailingOnly = TRUE) # get command line arguments as detailed below
cell<-args[1] # Choose cell
source.switch<-args[2]  # choose source switch type (immediate, delayStep, expSwitch)
sigV<-as.numeric(args[3]) # Choose source term (as a fraction! i.e. 40% is 0.4)
model<-args[4]
availCores<-as.numeric(args[5])

fitUp<-1
zero.alphaB<-0
nLogicalCores<-min(nrow(modelConfigs),availCores)
#######################################################################################
#slope select (cowabunga dude!)
justDown<-0
source('workDir.R')

boot.exp.folder<-paste(folderEMain,"bootstrap/",sep="")
if(!dir.exists(boot.exp.folder)){
  print("No boot folder dummmy!")
  dir.create(boot.exp.folder)
}

#######################################################################################
library(iterators)
library(parallel)
library(doParallel)
library(foreach)
library(Rcpp)
library(xlsx)

# This code assumes that Analysis.R has been run at least once (probably with have.CIs = FALSE) in order to generate the BestCompare.xlsx
# This .xlsx file is needed to extract what the best fit kBoxMax and bBoxMax are for each sigma, which this code needs to know (to import the correct files)
# Make sure the BestCompare.xlsx has entries corresponding to the source terms you want to fit! 
summary.file<-paste(folderEMain,"summaries/",cell,model,source.switch,"BestCompare.xlsx",sep="")
summary.data<-read.xlsx(summary.file,sheetIndex=1)


# find best-fit box configuration for given source term
best.data<-summary.data[round(summary.data$sig,5)==sigV,]

#select box config and sigma list
bBoxMax<-best.data["bBox"][1,]
kBoxMax<-best.data["kBox"][1,]


impFile<-paste("Data/",cell,"BRDU.xls",sep="")

impFile<-paste("Data/",cell,"BRDU.xls",sep="")
if(model=='kinHetExtended'||model=='kinHetExtended2'){
  folderE<- paste(folderEMain,model,"-",100*eval(sigV),"percSource/",sep="")
} else {
  folderE<- paste(folderEMain,model,100*eval(sigV),"percSource/",sep="")
}

bestFitFile<-paste(folderE,"/B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"bestfit.xlsx",sep="")
bestFitPars<-read.xlsx(bestFitFile,sheetIndex=1)
bestVals<-bestFitPars[,2]
bestNames<-as.vector(bestFitPars[,1])
bestFitPars<-bestVals
names(bestFitPars)<-bestNames

nuV<- -2 #setting nuV=-2 sets it to be free and fitted, all other values taken to be fixed inputs. if source.switch=immediate this doesn't matter, nu doesn't appear

#import file for functions etc from Mathemagica
fileI<-paste("Functions/",source.switch,"/B",bBoxMax,"K",kBoxMax,model,source.switch,sep="")
E<-exp(1)
#choose which SSR function to use (kinetic or temporal heterogeneity. In principle could make totally elastic in terms of abstracting to arbitrary numbers of populations, but am lazy for now)
fn<-paste("computeSSR",model,sep="")
fn<-parse(fn,text=fn)
#get data. Duh.
source("getData.R", local = TRUE)
phiGlobal<-mean(ki67Data) # global best fit ki67+%
#get ODEs...double duh
source("getODEs.R", local = TRUE)
#you get the point
source("getFunctions.R", local = TRUE)
source("getODEsolver.R", local = TRUE)
computeSSR<-eval(fn)
source("getSettings.R", local = TRUE)
# invisible(rm(i))
headerSave<-header
usefulPars<-bestFitPars[header]
############################################################################################################
#check best fit pars give same LogL
SSR<-computeSSR(usefulPars,tauUp1,tauDownList,UpDataList,DownDataList)
if(abs(bestFitPars[["lnL"]]-SSR)>10^-4){
  print("You dun gawn fucked up! Old SSR!=current SSR - check!")
} else {
  # ############################################################################################################
  invTform<-function(x){
    (sin(x))^2
  }
  
  sample.rows<-function(frame,number){
    frame[sample(nrow(frame),number,replace=TRUE),] # 30/nov/15 noticed that sampling wasn't with replacement here!!!
  }
  
  if(model=="kinHet"||model=="kinHetExtended"||model=='kinHetExtended2'){
    pars<-calcParsKinHet(usefulPars)
    parsA<-pars[[1]]
    parsB<-pars[[2]]
    functionVals<-calcFunctionDoublePop(parsA,parsB,tauUp1,tauDownList,downStartUpPos)
    
  } else if(model=="tempHet"){
    pars<-calcParsTempHet(usefulPars)
    
    functionVals<-calcFunctionSinglePop(pars,tauUp1,tauDownList,downStartUpPos)
  }  
  
  #best fit Km and Kp fractions
  fracs<-lapply(functionVals,computeFinalFracs)
  fracsSSR<-fracs
  for(kk in 1:nUp){
    fracs[[kk]][1:2,]<-0 
  }
  dataList<-c(UpDataList,DownDataList)
  bestFitResiduals<-calcResiduals(fracs,dataList[1:nUp],dataList[nUp+1:nDown]) #data-function
  #   # ############################################################################################################
  do.boot<-function(){
    newData<-list() # contains list of 2-column dataframes containing residuals for Km and Kp curves respectively
    for(jj in 1:length(fracs)){   
      curve.resids<-data.frame(bestFitResiduals[[2*jj-1]],bestFitResiduals[[2*jj]])
      choice<-sample.rows(curve.resids,nrow(curve.resids)) #random choice of mice but keeping kM and kP data together as rows!
      # residuals are defined on the transformed scale, but the best fit is not! need to transform best fit curve THEN add residual back
      tempKm<-Tform(fracs[[jj]][,1])+choice[,1] #bestfit+residual
      tempKp<-Tform(fracs[[jj]][,2])+choice[,2]
      # need to transform the curves back - our routine needs to take in 'untransformed' data
      temp<-cbind(invTform(tempKm),invTform(tempKp))
      # finally, save new replicate data to be fed into the ordinary routine (remember though that we've redefined phiGlobal)
      colnames(temp)<-c("kMinusFrac","kPosFrac")
      newData[[jj]]<-temp
    }
    #try using smaller number of iterations in bootstrap
    doFit<-findFit(tauUp1,tauDownList,newData[1:nUp],dataList[nUp+1:nDown])
    names(doFit)<-header

    source("formatOutput.R", local = TRUE)

    format<-paste(model,"Format",sep="")
    format<-eval(parse(format,text=format))
    format(dataList[1:nUp],dataList[nUp+1:nDown],0)# 0 = don't output plots, takes up space and only needed for diagnostics really.

    return(doFit)
  }
  #   # ############################################################################################################
  
  cl <- makeCluster(nLogicalCores)
  registerDoParallel(cl)

  
  #export file name (well, part of it - is pasted together with other things, see formatOutput.R)
  folderE<-paste(boot.exp.folder,"B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"/",sep="")
  #   # ############################################################################################################
  
  if(!file.exists(folderE)){
    dir.create(file.path(folderE))
  }
  print("start boot strap")
  ############################################################################################################
  boot.results<- foreach(kk = 1:iters,.export=c("cell","impFile","upSets","downSets","justDown","fitUp","zero.alphaB","workDir")) %dopar% {
    
    fileE<-paste(folderE,"B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"Replicate",toString(kk),sep="")
    BF.file<-paste(fileE,"bestfit.xlsx",sep="")
    if(!file.exists(BF.file)){
      nuV<- -2 
      library(GenSA)
      #import file for functions etc from Mathemagica
      fileI<-paste("Functions/",source.switch,"/B",bBoxMax,"K",kBoxMax,model,source.switch,sep="")
      E<-exp(1)
      fn<-paste("computeSSR",model,sep="")
      fn<-parse(fn,text=fn)
      source("getData.R", local = TRUE)
      #resample phiGlobal here - should be done with paired mice but just a quick code change here as a first pass.
      phiGlobal<-mean(sample(ki67Data,length(ki67Data),replace=TRUE))
      source("getODEs.R", local = TRUE)
      source("getFunctions.R", local = TRUE)
      source("getODEsolver.R", local = TRUE)
      computeSSR<-eval(fn)
      #for some reason had to change settings file rather than assign (globally) new maxit/temp, annoying but whatever works! Formally gives more freedom to change things but weird still...
      source("getBootSettings.R", local = TRUE)
      
      # noted on 19/5 that thing about error handling below is bollocks
      #custom error handling. seems sometimes the choice of resdiuals can lead to 0 or infinity-order parameters being solved by nleqslv
      #this fucks shit up so we hack a while loop to keep choosing residuals until we get some that don't give an error!
      #20/11/15 doesn't seem to work properly...still getting 'knockedParsA/B' errors, needs to be in actual nleqslv function so left for now :/
      do.boot()
    } else if(file.exists(BF.file)){
      #do nothing - skip this replicate
    }
    print("Done.")
    ############################################################################################################
  } # end foreach
} #end if statement checking SSR against best fit vals
