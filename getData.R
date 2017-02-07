############################################################################################################
# The code below is run by both the main best fit and bootstrap fitting routines to import the prepared
# data for fitting.
############################################################################################################

require(gdata)
############################################################################################################
#get mean global ki67 pos percentage
ki67File<-paste("Data/",cell,"AllKi67Pos.xls",sep="")
ki67Data<-as.numeric(read.xls(ki67File, sheet="Sheet1",header=F)[[1]])
############################################################################################################

nDown<-length(downSets)
nUp<-length(upSets)
#almost certainly easier way to do this but meh!
kmL<-"KmRatio"
kpL<-"KpRatio"

appendi<-function(string,i){
  paste(string,toString(i),sep="")
}

assignData<-function(slope,dSet){
  kmT<-paste(kmL,slope,sep="")
  kpT<-paste(kpL,slope,sep="")
  tauL<-paste("tau",slope,sep="")

  finalVar<-paste(slope,"Data",sep="")
  counter<-1
  for(i in 1:length(dSet)){

    mVals<-read.xls(impFile, sheet=paste(kmT,dSet[i], sep=""),header=F)
    pVals<-read.xls(impFile, sheet=paste(kpT,dSet[i], sep=""),header=F)
    out<-cbind(mVals[,2],pVals[,2])
    
    assign(appendi(finalVar,counter),out, envir = .GlobalEnv)
    assign(appendi(tauL,counter),mVals[,1], envir = .GlobalEnv)
    counter<-counter+1
  }
}
# ############################################################################################################
#definitions here: internal definitions of down and upslope now just go from 1:nDown. E.g., if we select downSets<-c(2,3)
#these will be mapped to DownData1 and DownData2 - internal counters only keep track of data IN USE!! Labels of which curves are which are checked at the door
assignData("Up",upSets)
assignData("Down",downSets)

UpDataList<-vector("list", nUp)
DownDataList<-vector("list", nDown)

#Up/DownDataLists contain all data in use according to the numbering scheme described in the comment above
for(i in 1:nUp){
  UpDataList[[i]]<-eval(parse(text=appendi("UpData",i)))
}

for(i in 1:nDown){
  DownDataList[[i]]<-eval(parse(text=appendi("DownData",i)))
}

############################################################################################################

join<-function(stringList){
  temp<-NULL
  for(i in 1:length(stringList)){
    temp<-paste(temp,stringList[i],sep="")
  }
  return(parse(text=temp))
}

############################################################################################################
#currently working on the assumption that the downslope start times are included in the upslope data case too
#this means we can save computer time and just solve the upslope at certain points and use these as initial conditions for the downslope,
#without having dedicatied solver routines to computer some special downslope time point

tauL<-"tauDown"
popL<-c("A","B")
tauDownList<-list()
downStartUpPos<-list()#stores position in upTimes of downStartTimes to be used in calculation of ICs
downStartTimes<-list()

labelsKm<-list()
labelsKp<-list()

for(i in 1:nUp){
  labelsKm[[i]]<-paste("upKm",toString(upSets[i]),sep="")
  labelsKp[[i]]<-paste("upKp",toString(upSets[i]),sep="")
}


for(i in 1:nDown){
  tau<-eval(join(c(tauL,i)))
  tauDownList[[i]]<-(tau)  
  downStartUpPos[[i]]<-match(min(tau),tauUp1) #IMPORTANT!!! - here we assume that both upslopes have the same exact time points (basically so we know where to look for matching conditions for downslope initial conditions)
  downStartTimes[[i]]<-min(tau)
  #need to shift by nUp because it's a global list including up and downslope data
  labelsKm[[i+nUp]]<-paste("downKm",toString(downSets[i]),sep="")
  labelsKp[[i+nUp]]<-paste("downKp",toString(downSets[i]),sep="")
}

downStartUpPos<-unlist(downStartUpPos)
downStartTimes<-unlist(downStartTimes)
############################################################################################################
#keep track of number of points

nrows<-list()
counter<-1
for(i in 1:nUp){
  nrows[[counter]]<-nrow(UpDataList[[i]])
  nrows[[counter+1]]<-nrows[[counter]]
  counter<-counter+2
}

for(i in 1:nDown){
  nrows[[counter]]<-nrow(DownDataList[[i]])
  nrows[[counter+1]]<-nrows[[counter]]
  counter<-counter+2
}
nrows<-unlist(nrows)
#old def of nrows just for reference
#nrows<-c(rep(length(tauUp1),2),rep(length(tauDown2),2),rep(length(tauDown3),2))

#############################################################################################
#horribly coded attempt to keep track of unique time points to avoid lots of duplicate NDSolving
tauUpUnique1<-unique(tauUp1)
tauDownUniqueList<-lapply(tauDownList,unique)

rm(i,counter) # be sure to delete common iterators so we are aware of nasty nested loops that re-use iterators!
print('getData internal')
