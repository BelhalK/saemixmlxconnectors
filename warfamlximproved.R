source('R/aaa_generics.R')
source('R/compute_LL.R')
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 

source('R/main.R')
source('R/main_estep.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 
source('R/SaemixObject.R') 
source('R/zzz.R') 

library("mlxR")
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix")

# #convert line endings '\r' to '\t'
# warfa_data <- read.table("data/warfarin_data.txt", header=T)
# write.table(warfa_data,"data/warfa.txt",sep="\t",row.names=FALSE)


################################################################ MonolixProject ####################################################################################################################################
#data
data = list(dataFile= "data/warfarin_data_monolix.txt",
    headerTypes =c("id", "time", "amount", "observation", "obsid", "ignore", "ignore","ignore"),
    observationTypes = list(dv = "continuous"))
warfarin.saemix <- read.table(data$dataFile, header=T)
# warfarin.saemix <- warfarin.saemix[which(warfarin.saemix$dvid==1),]

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time","amt"),name.response=c("dv"), name.X="time")

#model
model.file <- 'mlxProjects/warfarinmlx/model/oral1_1cpt_kaVCl.txt'
saemix.model<-saemixModel(model=model.file,description="warfarin",type="structural", monolix=TRUE,
  ,psi0=matrix(c(1,1,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))
#saemix Run
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(100,50),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,nb.chains=1)
warfa<-saemix(model=saemix.model,data=saemix.data,options_warfa)

################################################################ MonolixProject ####################################################################################################################################




################################################################ SAEMIX ####################################################################################################################################
warfarin.saemix <- read.table("data/warfa.txt", header=T)
saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")
dim(warfarin.saemix)

# warfarin.saemix <- read.table(data$dataFile, header=T)
# warfarin.saemix <- warfarin.saemix[which(warfarin.saemix$dvid==1),]
# saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("time","amt"),name.response=c("dv"), name.X="time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k<-Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}


saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,1,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE), monolix=FALSE)

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(100,50),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,monolix=FALSE)
warfa<-saemix(saemix.model,saemix.data,options_warfa)

################################################################ SAEMIX ####################################################################################################################################





# data = list(dataFile= "data/warfarin_data_mlx.txt",
#     headerTypes =c("id", "time", "amount", "observation", "ignore", "ignore", "ignore","ignore"),
#     observationTypes = list(dv = "continuous"))
# warfarin.saemix <- read.table(data$dataFile, header=T)
# saemix.data<-saemixData(name.data=warfarin.saemix, data.mlx = data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("dv"), name.X="time")
