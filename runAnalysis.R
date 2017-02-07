rm(list = ls())
############################################################################################
# Analysis file for outputs of fitting procedures.

# Output is the BestCompare.xlsx files in /summaries/ as well as other files in /bulk/

############################################################################################
# BEGIN USER EDIT AREA

target.source.list<-c(7) #in percent, not fraction
upSets<-c('Expt185','Expt203')
downSets<-c('Expt184','Expt196','Expt197')

# END USER EDIT AREA
############################################################################################

library(foreach)
library(xlsx)
library(gdata)
library(lattice)
library(grid)
library(gridExtra)
library(ggplot2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
model<-args[1]
cell<-args[2]
source.switch<-args[3]
modelLabel<-model
source('workDir.R')

mainDir<-outputFolder
summary.folder<-paste(mainDir,cell,"/summaries/",sep="")
if(!dir.exists(summary.folder)){
  print("Summaries folder does not exist, creating:")
  dir.create(summary.folder)
  dir.create(paste(summary.folder,"bulk/",sep=""))
}

expDir<-paste(mainDir,cell,"/summaries/",sep="")
paste.vectors<-function(name.vector){
  c(paste(name.vector,".lower",sep=""),paste(name.vector,".upper",sep=""))
}

############################################################################################
#choose max box configuration
source("getModelConfigs.R", local = TRUE)
############################################################################################
############################################################################################
#Loop over source/source.switch types
sourceList<-target.source.list
nModels<-length(sourceList)
#combine source term with source.switch function to get all possible models under consideration
source.config<-matrix(-1, nrow =nModels, ncol = 2)

for(sourceC in 1:length(sourceList)){
  source.config[sourceC,1] <- sourceList[sourceC]
  source.config[sourceC,2] <- source.switch	
}


############################################################################################
#Keep track of best fit toplogy across models
bestPerSig<-NULL

############################################################################################
#Start model loop
for(modelIterator in 1:nrow(source.config)){
  sig.size<-source.config[modelIterator,1]

  sigV<-eval(parse(text=sig.size))/100. # need to divide by 100 here because our sigma Lists as defined by hand above in this file are in % and need to be fractions.
  # As a comparison: see the definition of source.frac in analysisHeaders.R - the sigma there is pulled directly from summary.xlsx and so is of the form 0.05 rather than 5% (for example)
  # Can't have this file pull from the summary files, this file *makes* the summary files!
  source.switch<-source.config[modelIterator,2]
  source("analysisHeaders.R")
  
  impFile<-paste("Data/",cell,"BRDU.xls",sep="")
  if(model=='kinHetExtended'||model=='kinHetExtended2'){
    # folderE<- paste(folderEMain,model,"-",100*eval(sigV),"percSource/",sep="")
    labelCommon<-paste(model,'-',sig.size,"percSource",sep="")
  } else {
    labelCommon<-paste(model,sig.size,"percSource",sep="")
    # folderE<- paste(folderEMain,model,100*eval(sigV),"percSource/",sep="")
  }

  # labelCommon<-paste(model,sig.size,"percSource",sep="")
  fileCommon<-paste(mainDir,cell,"/",labelCommon,"/",sep="")
  ############################################################################################
  output=NULL
  for(boxModelIterator in 1:nrow(boxConfigs)){
    
    k<-boxConfigs[boxModelIterator,1]
    b<-boxConfigs[boxModelIterator,2]
    file<-paste(fileCommon,"B",toString(b),"K",toString(k),model,source.switch,"bestfit.xlsx",sep="")
    fit<-read.xlsx(file, sheetIndex=1)
    header<-unlist(lapply(fit[[1]],as.character)) #only needs to be defined once but meh, if statement just as costly
    vals<-as.data.frame(rbind(fit[[2]]))
    
    header<-c(header,"bBox","kBox")
    vals<-cbind(vals,b,k) #add box configs here
    
    output<-rbind(output,vals)
  }#end for loop
  names(output)<-header
  
  #function defined in analysisHeaders.R: takes in raw excel data and formats data by making global/population averaged values
  compiled.pars<-get.parameters(header,output,model,sigV)
  output<-compiled.pars[[1]]
  header<-compiled.pars[[2]]

  output<-output[order(output$lnL),]
  
  targetOutputName<-"bestPerSig"
  targetOutput<-eval(parse(text=targetOutputName)) 
  
  assign(targetOutputName,rbind(targetOutput,output[1,]))
  
  write.xlsx(output,paste(expDir,"bulk/",cell,labelCommon,source.switch,"SUMMARY.xlsx",sep=""))
  
  ############################################################################################
  #end model loop
}

targetOutputName<-"bestPerSig"
targetOutput<-eval(parse(text=targetOutputName)) 

write.xlsx(targetOutput,paste(expDir,cell,model,source.switch,"BestCompare.xlsx",sep=""))

print('Done.')
