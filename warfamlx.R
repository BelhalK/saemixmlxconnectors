source('oldR/aaa_generics.R')
source('oldR/compute_LL.R')
source('oldR/func_aux.R') 
source('oldR/func_distcond.R') 
source('oldR/func_FIM.R')
source('oldR/func_plots.R') 
source('oldR/func_simulations.R') 

source('oldR/main.R')
source('oldR/main_estep.R')
source('oldR/main_initialiseMainAlgo.R') 
source('oldR/main_mstep.R') 
source('oldR/SaemixData.R')
source('oldR/SaemixModel.R') 
source('oldR/SaemixRes.R') 
source('oldR/SaemixObject.R') 
source('oldR/zzz.R') 

library("mlxR")
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix")


# if no Monolix project, you can start by creating one from R directly
data = list(dataFile= "data/warfarin_data_mlx.txt",
headerTypes =c("id", "time", "amount", "observation", "ignore", "ignore", "ignore","ignore"),
observationTypes = list(dv = "continuous"))
modelFile <- 'models/pkmodel.txt'

#create mlx project
newProject(modelFile = modelFile, data = data)
#save mlx project
saveProject("/path/to/file.mlxtran")


################################################################ SAEMIX ####################################################################################################################################
project.file <- "mlxProjects/warfarinmlx/warfarinPK_project.mlxtran"  #your mlx project
loadProject(project.file)

warfa_data <- readDatamlx(project = project.file) 
# OR READ IT DIRECTLY FROM THE .txt FILE
warfarin.saemix <- read.table("data/warfarin_data.txt", header=T)

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k<-Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}


saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,1,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

##RUNS
K1 = 200
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,monolix=FALSE)
warfa<-saemix(saemix.model,saemix.data,options_warfa)

################################################################ SAEMIX ####################################################################################################################################


################################################################ MonolixProject ####################################################################################################################################
project.file <- "mlxProjects/warfarinmlx/warfarinPK_project.mlxtran"
loadProject(project.file)
warfa_data <- readDatamlx(project = project.file)
treat <- warfa_data$treatment[,c(1,3)]
warfarin.saemix <- merge(treat ,warfa_data$y_1,by="id")
warfarin.saemix <- warfarin.saemix[order(warfarin.saemix$id),]


saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  ypred <- 2  #dummy model
  return(ypred)
}

#SINCE monolix=TRUE model1cpt is a dummy model
#The image of the structural model will be directly computed from Monolix
#The saemixModel is useful for all the param values, only its 'model' argument is not used

saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,1,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,nb.chains=1,monolix=TRUE)
warfa<-saemix(model=saemix.model,data=saemix.data,options_warfa)



################################################################ MonolixProject ####################################################################################################################################





