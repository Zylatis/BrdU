rm(list = ls())
############################################################################################

# Main code to loop over box configurations and perform best fit calculations

# Outputs best fit files in experimentt/cell/source-hetmodel/ directory


############################################################################################
# BEGIN USER EDIT AREA

upSets<-c('Expt185','Expt203')
downSets<-c('Expt184','Expt196','Expt197')

# END USER EDIT AREA
############################################################################################

fitUp<-1
zero.alphaB<-0


# Currently taking command line arguments
args <- commandArgs(trailingOnly = TRUE)
cell<-args[1]
source.switch<-args[2]
model<-args[3]
availCores<-as.numeric(args[4])
source('workDir.R')
#slope select (cowabunga dude!)
justDown<-0

print(cell)
print(model)
print(source.switch)

library(iterators)
library(parallel)
library(doParallel)
library(foreach)
library(Rcpp)

#choose max box configuration
source("getModelConfigs.R", local = TRUE) #this now gets possible source configurations now

############################################################################################

print(nrow(modelConfigs))
nCores<-min(nrow(modelConfigs),availCores)
print(nCores)
cl <- makeCluster(nCores)
registerDoParallel(cl)
############################################################################################
# # #START PARALLEL LOOP HERE   ############################################################################################
out<-foreach(ModelIterator=1:nrow(modelConfigs), .combine=rbind,.export=c("upSets","downSets","fitUp",'zero.alphaB')) %dopar%  { # ,.errorhandling = 'pass' 
  E<-exp(1)
  kBoxMax<-modelConfigs[ModelIterator,1]
  bBoxMax<-modelConfigs[ModelIterator,2]
  sigV<-modelConfigs[ModelIterator,3]
  
  #######################################################################################
  #other preliminaries
  #######################################################################################
  #data file
  impFile<-paste("Data/",cell,"BRDU.xls",sep="")
  if(model=='kinHetExtended'||model=='kinHetExtended2'){
    folderE<- paste(folderEMain,model,"-",100*eval(sigV),"percSource/",sep="")
  } else {
  folderE<- paste(folderEMain,model,100*eval(sigV),"percSource/",sep="")
  }
  print(folderE)

  # Export file management bits (makes folders etc)
  modelPrint<-paste("B",toString(bBoxMax),"K",toString(kBoxMax)," ", model," ",source.switch," sig",100*sigV,sep="")
  if(!file.exists(folderE)){
    dir.create(file.path(folderE))
  }
  fileE<-paste(folderE,"B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,sep="")
  #check if this configuration already computed NOTE: this says nothing of whether or not it was computed with the same internal solver settings!
  ######################################################################################
  bestFileE<-paste(fileE,"bestfit.xlsx",sep="")

  # Check if best fit.xlsx file already exists
  if(!file.exists(bestFileE)){
    # Note: nu has been kept free for almost all current fits, so this nuV freedom is somewhat vestigial (ignore!)
    nuV<- -2 #setting nuV=-2 sets it to be free and fitted, all other values taken to be fixed inputs. if source.switch=immediate this doesn't matter, nu doesn't appear
    #export file name (well, part of it - is pasted together with other things, see kinFormat.R and tempFormat.R)

    #import file for functions etc from Mathemagica
    fileI<-paste("Functions/",source.switch,"/B",bBoxMax,"K",kBoxMax,model,source.switch,sep="")
    
    #choose which SSR function to use (kinetic or temporal heterogeneity. In principle could make totally elastic in terms of abstracting to arbitrary numbers of populations, but am lazy for now)
    fn<-paste("computeSSR",model,sep="")
    fn<-parse(fn,text=fn)
    #get data. Duh.
    source("getData.R", local = TRUE)

    #as of 30/nov/2015 getDataMultidown now only imports the set of Ki67%'s, the mean is calcualted here to allow it to be more easily bootstrapped in the bootstrap code.
    #phiglobal must be declared before getSettings.R is sourced!
    phiGlobal<-mean(ki67Data)
    
    #get ODEs...double duh
    source("getODEs.R", local = TRUE)
    
    #you get the point
    source("getFunctions.R", local = TRUE)

    source("getODEsolver.R", local = TRUE)

    # Define the SSR function we wish to use (tempHet or kinHet, see above definition of 'fn')
    computeSSR<-eval(fn)
    source("getSettings.R", local = TRUE)
    headerSave<-header
  ######################################################################################
    #run fit
#     options(error = recover)
#     options(warn = 2)

    # Function findFit defined in getFunctionsMultiDown.R
    # We simply pass it the upslope times, the downslope times, and the up and downslope data (see getDataMultiDown.R for data structure details)
    doFit<-findFit(tauUp1,tauDownList,UpDataList,DownDataList)
    names(doFit)<-header

    #format results and output images and .xls files
    source("formatOutput.R", local = TRUE)
    format<-paste(model,"Format",sep="")
    format<-eval(parse(format,text=format))
    format(UpDataList,DownDataList,1)#1 correspods to 'yes' for outputting plots as well as best fit data. want plots for best fit, sometimes (usually not) for bootstraps

    tracker<-cbind(modelPrint,"done")
    write.table("done", file=paste("progress/done",modelPrint, sep=""), sep="\t", row.names=F)
######################################################################################
#end file existence if -statement
  } else {
    tracker<-cbind(modelPrint,"skipped")
   
   }
# ######################################################################################
   tracker

}
stopCluster(cl)

############################################################################################
#END PARALLEL LOOP HERE
############################################################################################
