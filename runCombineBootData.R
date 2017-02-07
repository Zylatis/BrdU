rm(list = ls())


library(foreach)
library(xlsx)
library(gdata)
library(lattice)
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)
library(doParallel)
################################################################################################################################################################
# This code assumes that Analysis.R has been run in order to generate the BestCompare.xlsx
# This .xlsx file is needed to extract what the best fit kBoxMax and bBoxMax are for each sigma, which this code needs to know (to import the correct files)
# Alternatively, one could make it so this file simply picks a sigma then looks to see what bootstrap files exist, but the above method already exists for the
# actual bootstrap code (as it must, to generate the bootstrapped fit files) and so we just recycle that idea here

################################################################################################################################################################

args <- commandArgs(trailingOnly = TRUE)
model<-args[1]
cell<-args[2]
availCores<-args[5]

source.switch<-"delayStep"
n.cores<-min(nrow(modelConfigs),availCores) #detectCores(all.tests = FALSE, logical = TRUE)-1
nReplicates<- 1500 #must be manually set!
print(model)

upSets<-c('Expt185','Expt203')
downSets<-c('Expt184','Expt196','Expt197')
source('workDir.R')


ki67File<-paste("Data/",cell,"AllKi67Pos.xls",sep="")
ki67Data<-as.numeric(read.xls(ki67File, sheet="Sheet1",header=F)[[1]])
phiGlobal<-mean(ki67Data)
# analysisHeader.R requires a source.switch and model to be defined, so only contains functions that are used in situations where those variables exist
# i.e. it contains some function definitions, but actually also contains code which is RUN on source, and will throw errors if certain variables are
# missing.
source("analysisHeaders.R")

summary.file<-paste(outputFolder,cell,"/summaries/",cell,model,source.switch,"BestCompare.xlsx",sep="")
summary.data<-read.xlsx(summary.file,sheetIndex=1)

final.set<-foreach(i = 1:nrow(summary.data),.combine='rbind') %do% { #this must remain in serial if the next (nested) loop is dopar
################################################################################################################################################################
  # Here we are actually looping over whatever configurations are in the BestCompare.xlsx file - no freedom to pick and choose specific source terms!
  # If you want to limit which source terms you examine, you must either hack this code or re-run Analysis.R to generate a new BestCompare.xlsx file
  # This is done with good reason: several files make use of a loop over sigmas and bBoxes etc, having them all have their own toggles to pick
  # and mix is a recipe for disaster! Better to keep everything strongly coupled to one set of toggles...
  bBoxMax<-summary.data["bBox"][i,]
  kBoxMax<-summary.data["kBox"][i,]
  sigV<-summary.data["sig"][i,]

  folder<-paste(outputFolder,cell,"/bootstrap/B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"/",sep="")

  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  boot.vals<-foreach(i = 1:nReplicates,.combine='rbind') %dopar% {
    output<-NULL
    
    library(xlsx)
    file<-paste(folder,"B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"Replicate",toString(i),"bestfit.xlsx",sep="")
    vals<-read.xlsx(file,sheetIndex=1)
    header<-unlist(lapply(vals[[1]],as.character)) #only needs to be defined once but meh, if statement just as costly
    vals<-as.data.frame(rbind(vals[[2]]))
    output<-rbind(output,vals)
    #function defined in analysisHeaders.R: takes in raw excel data and formats data by making global/population averaged values

    compiled.pars<-get.parameters(header,output,model,sigV)
    output<-compiled.pars[[1]] # header not required here
  }
    stopCluster(cl)
    
    export.file.allboot<-paste(outputFolder,cell,"/bootstrap/","B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,"sig",toString(sigV*100),"AllBootVals.xlsx",sep="")
    write.xlsx(boot.vals,export.file.allboot)
     
################################################################################################################################################################

    
make.CI<-function(col.vector){
  col.vector<-col.vector[,1]
  lower<-quantile(col.vector,0.025)
  upper<-quantile(col.vector,0.975)
  data.frame(lower=lower[[1]],upper=upper[[1]])
}

#this allows us to just define a vector of parameters and loop over it to make the CI data
set.CI<-function(parameter){#parameter needs to be string
  parameter.CI<-make.CI(boot.vals[parameter])
  label<-paste(parameter,".CI",sep="")
  assign(label,parameter.CI,envir=.GlobalEnv)
}



# see analysisHeaders.R for definition of all/global pars
all.pars<-c(global.pars,other.pars,implicit.plot.pars)
CI.headers<-unlist(lapply(all.pars,paste.vectors))
 
compiled.CI.data<-as.data.frame(t(unlist(lapply(all.pars,set.CI))))
names(compiled.CI.data)<-CI.headers

fixed.data<-data.frame(sig=sigV,model=model,bBoxMax=bBoxMax,kBoxMax=kBoxMax,source.switch=source.switch)

result<-cbind(fixed.data,compiled.CI.data)
################################################################################################################################################################
} #end foreach

CI.export.file<-paste(outputFolder,cell,"/bootstrap/",model,source.switch,"ParameterCIs.xlsx",sep="")
write.xlsx(final.set,CI.export.file,row.names=FALSE)
################################################################################################################################################################
################################################################################################################################################################
