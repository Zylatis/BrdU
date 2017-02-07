
workDir<-'/home/graeme/Dropbox/BrdUFiles/'
#workDir<-'C:/Users/Graeme/Dropbox/BrdUFiles/'
setwd(workDir)
current.data.strings <- c(upSets,downSets)
outputFolder <- paste(workDir,sep="")

for(i in 1:length(current.data.strings)){
  outputFolder<-paste(outputFolder,current.data.strings[i],sep="")
}
outputFolder<-paste(outputFolder,'/',sep="")
folderEMain<-paste(outputFolder,cell,"/",sep="")

if(!dir.exists(outputFolder)){
  print("No output folder dummmy!")
  dir.create(outputFolder)
}

if(!dir.exists(folderEMain)){
  dir.create(folderEMain)
  dir.create(paste(folderEMain,'summaries',sep=""))
  dir.create(paste(folderEMain,'summaries/bulk/',sep=""))
}
