#######################################################################################
# Code to plot the Brdu+ frac in Ki67+- for a single box-config. Requires CIs to exist
# i.e. bootstraps must have been run for this configuration.

# Produces plot outputs on both the transformed and un-transformed scale.

# Has optional function calls (at end of file) to perform consistency checks for this
# set of parametres, i.e. make sure that the ki67+ frac across each level of ki67 (each
# ki67 box) is constant with time.
#######################################################################################
args <- commandArgs(trailingOnly = TRUE)
cell<-args[1]
sigV <-as.numeric(args[2])
model<- args[3]


source.switch<-"delayStep"
fitUp<-1
zero.alphaB<-0
#######################################################################################
#slope select (cowabunga dude!)
justDown<-0

upSets<-c('Expt185','Expt203')
downSets<-c('Expt184','Expt196','Expt197')

source('workDir.R')
#######################################################################################
library(iterators)
library(parallel)
library(doParallel)
library(foreach)
library(Rcpp)
library(xlsx)

summary.file<-paste(folderEMain,"summaries/",cell,model,source.switch,"BestCompare.xlsx",sep="")
summary.data<-read.xlsx(summary.file,sheetIndex=1)
cur.focus <- summary.data[round(summary.data$sig,5) == sigV, ]

bBoxMax<-cur.focus['bBox'][[1]]
kBoxMax<-cur.focus['kBox'][[1]]

#data file
impFile<-paste("Data/",cell,"BRDU.xls",sep="")
if(model=='kinHetExtended'||model=='kinHetExtended2'){
fit.folder<-paste(outputFolder,"/",cell,"/",model,"-",100*eval(sigV),"percSource",sep="")
folderE<-paste(outputFolder,"/",cell,"/Examine/",model,"-",100*eval(sigV),"percSource/",sep="")

} else {
  fit.folder<-paste(outputFolder,"/",cell,"/",model,100*eval(sigV),"percSource",sep="")
  folderE<-paste(outputFolder,"/",cell,"/Examine/",model,100*eval(sigV),"percSource/",sep="")
  
}
modelPrint<-paste("B",toString(bBoxMax),"K",toString(kBoxMax)," ", model," ",source.switch," sig",100*sigV,sep="")

fit.file<-paste(fit.folder,"/B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,sep="")
#check if this configuration already computed NOTE: this says nothing of whether or not it was computed with the same internal solver settings!
######################################################################################
fit.file.pars<-paste(fit.file,"bestfit.xlsx",sep="")
nuV<- -2 #setting nuV=-2 sets it to be free and fitted, all other values taken to be fixed inputs. if source.switch=immediate this doesn't matter, nu doesn't appear
#export file name (well, part of it - is pasted together with other things, see kinFormat.R and tempFormat.R)

#import file for functions etc from Mathemagica
fileI<-paste("Functions/",source.switch,"/B",bBoxMax,"K",kBoxMax,model,source.switch,sep="")
E<-exp(1)
#choose which SSR function to use (kinetic or temporal heterogeneity. In principle could make totally elastic in terms of abstracting to arbitrary numbers of populations, but am lazy for now)
fn<-paste("computeSSR",model,sep="")
fn<-parse(fn,text=fn)
#get data. Duh.
source("getData.R", local = TRUE)
phiGlobal<-mean(ki67Data)
#get ODEs...double duh
source("getODEs.R", local = TRUE)
#you get the point
#source("getFunctionsMultiDownN.R", local = TRUE)
source("getFunctions.R", local = TRUE)
source("getODEsolver.R", local = TRUE)
computeSSR<-eval(fn)
source("getSettings.R", local = TRUE)

headerSave<-header
######################################################################################
f1<-paste(outputFolder,"/",cell,"/Examine/",sep="")
if(!dir.exists(f1)){
  dir.create(f1)
}
if(!dir.exists(folderE)){
  dir.create(folderE)
}
fileE<-paste(folderE,"/B",toString(bBoxMax),"K",toString(kBoxMax),model,source.switch,sep="")


fit<-read.xlsx(fit.file.pars, sheetIndex=1)
imported.header<-unlist(lapply(fit[[1]],as.character)) #only needs to be defined once but meh, if statement just as costly
imported.pars<-as.data.frame(rbind(fit[[2]]))

names(imported.pars)<-imported.header
fit.vector<-NULL
for(i in 1:length(header)){
  fit.vector[i]<-with(imported.pars,eval(parse(text=header[[i]])))
}


twice.log.likelihood<-computeSSR(fit.vector,tauUp1,tauDownList,UpDataList,DownDataList)
doFit<-fit.vector
names(doFit)<-header



source("singleSetCI.R",local=TRUE)
plot.with.CI(fit.vector,UpDataList,DownDataList)

source("formatOutput.R", local = TRUE)

format<-paste(model,"Format",sep="")
format<-eval(parse(format,text=format))
format(UpDataList,DownDataList,1)

header<-headerSave

#only works if kinHet
#check to see that the nleqslv routine isn't shitting its pants trying to get the numeric solutions to alpha_i and delta_i
#test.nleqslv(fit.vector)

# testConstancy(fit.vector,tauUp1,tauDownList)


