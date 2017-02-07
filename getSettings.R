############################################################################################################
# define basic header and GenSA search parameters for core parameters. Parameters are then added in the
# submodel part below. Currently sigma is not fitted, but its (constant) value can be changed manually
# in runBestFits.R


kinHetHeader<-c("b","eps","phiA") 

kinHetExtendedHeader<-kinHetHeader
kinHetExtended2Header<-kinHetHeader

tempHetHeader<-c("b","eps","kappa") 
noHetHeader<-c("b","eps") 

#choose from above pre-made headers based on model in use
headerN<-paste(model,"Header",sep="")
header<-eval(parse(headerN,text=headerN))

if(model=="kinHet"||model=="kinHetExtended"||model=="kinHetExtended2"){
  modelLabel<-"Kinetic heterogeneity"
  #basic starting set of parameters
  start<-c(2,   0.75,  phiGlobal+0.1)
  lower<-c(0.1, 0.25,   phiGlobal+10^-2) 
  upper<-c(3,   1,     1-10^-2)
  
  if(zero.alphaB==1){ 
    #if zero.alphaB=1, alpha_B = 0 so phiB is used in it's place when imposing conditions - it is not a free parameter.
  }
  
  if(zero.alphaB==0){
    header<-c(header,"phiB")
    start<-c(start,phiGlobal*0.1)
    lower<-c(lower,0.01)
    upper<-c(upper,phiGlobal-10^-2)
  }
  
}

if(model=="tempHet"){
  start<-c(2,0.75,0.1)
  lower<-c(0.5,0.5,10^-5)
  upper<-c(5,1,5)
  modelLabel<-"Temporal heterogeneity"
}

if(model=="noHet"){
  start<-c(2,0.75)
  lower<-c(0.5,0.5)
  upper<-c(3,1)
}

#by default we assume that nu is fixed, changed below and this fact is utilized in kinFormat. Shitty way of doing things but...
nuVstat<-"fixed"
############################################################################################################
# choose sub-model (i.e. free or fixed sigma and nu parameters, as governed by the search bounds)
if(eval(sigV)!=0&&(model=="kinHet"||model=="kinHetExtended"||model=="kinHetExtended2")){
  header<-c(header,"zeta") 
  start<-c(start,0.5)
  lower<-c(lower,10^-5) 
  upper<-c(upper,1-10^-5)
}

# only need zeta and source switch if source is nonzero!
if(eval(sigV)==0&&(model=="kinHet"||model=="kinHetExtended"||model=="kinHetExtended2")){
  zeta<-0
  nuV<-expression(-1) #-1 means not used, no changes to header means nu not seen by GenSA function
}

if(source.switch=="immediate"){
  nuV<-expression(-1) #-1 means not used, no changes to header means nu not seen by GenSA function
}


# NOTE: The reason that the delayStep and expSwitch cases are kept separate here is just the way that nu is defined (notice the range is different)
# In principle this should actually be made as in the paper so that the range of nu can be the same (1/nu vs nu etc) but this has not been done (changes introduce errors!)

# to force nu to a fixed value, set nuV = val in main.R (-2 will make it be free and fitted)
if(source.switch=="delayStep"&&eval(nuV)== -2){ 
  nuVstat<- "free"
  nuV<-expression(nu) 
  header<-c(header,"nu")
  # now nu is being seen by the GenSA and inputting into the SSR function. Need to add search range
  start<-c(start,0.5)
  lower<-c(lower,0)
  upper<-c(upper,5)
}

# to force nu to a fixed value, set nuV = val in main.R (-2 will make it be free and fitted)
if(source.switch=="expSwitch"&&eval(nuV)== -2){ 
  nuVstat<- "free"
  nuV<-expression(nu) 
  header<-c(header,"nu")
  
  start<-c(start,0.5)
  lower<-c(lower,0)
  upper<-c(upper,5)
}

############################################################################################################
# Define starting temperature and maximum number of iterations for GenSA package
temp<-10^8

if(model=="kinHet"||model=="kinHetExtended"||model=="kinHetExtended2"){
  maxit<-3*10^3
} 

if(model=="tempHet"){
  maxit<-10^4
}

if(source.switch == 'expSwitch'){
  maxit<-3*10^3
}

# Finally set parameter search limits
start<-as.numeric(start)
lower<-as.numeric(lower)
upper<-as.numeric(upper)
# Toggles allows us to keep track of different settings used for each fit. For each fit, along with the bestfit.xlsx
# we export an xlsx of all the numeric settings used in that fit.
toggles<-as.data.frame(rbind(start,lower,upper,temp,maxit))
names(toggles)<-header

solTolerance<-10^-4

# This code allows us to just fit the downslopes if we want. 
# We must always calcualte the upslope to get the ICs for down, but we can exclude it from the SSR if we want.
# The SSR is calculated by summing over a list of residuals. By starting this sum at 3 instead of 1 we exclude
# the Km and Kp upslope contributions from being included in the SSR. See computeSSRkinHet() in getFunctionsMultiDown.R.
if(fitUp==1){
  ssrStart<-1 #include upslope in SSR calc
} else if(fitUp==0){
  ssrStart<-2*nUp+1 #start at first downslope
}
print('settings internal')
