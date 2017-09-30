#----------------------------------------------------------------------------------------#
# K Lloyd 2016_09_16
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# All models applied to the same data. Models are:	GP, GPS1, GPS2, GPS3, AFT, Cox PH, 
# 													Coxnet, GBM and RSF
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Boxplot to compare concordance index results of different models
#---------------------------------------------------------------------------------------#

##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##
# args = commandArgs(trailingOnly=TRUE)

library(fields)
library(gbm)
library(glmnet)
library(ipred)
library(MASS)
library(Matrix)
library(nlme)
library(NORMT3)
library(pdist)
library(randomForestSRC)
library(rgl)
library(rms)
library(survcomp)
library(survival)
library(zoo)


##-------------------------------------------------------------------------------------##
##--------------------------------- Source Functions ----------------------------------##
##-------------------------------------------------------------------------------------##
source('../toSource/add.alpha.R')
source('../toSource/AdjustTrainingSurvivalMeanVariance.R')
source('../toSource/ApplyAFT.R')
source('../toSource/ApplyCoxph.R')
source('../toSource/ApplyGBM.R')
source('../toSource/ApplyCoxnet.R')
source('../toSource/ApplyGP.R')
source('../toSource/ApplyGPDiffInf.R')
source('../toSource/ApplyGPLearnCens.R')
source('../toSource/ApplyGPLaplace.R')
source('../toSource/ApplyRF.R')
source('../toSource/ApplyRFSurvival.R')
source('../toSource/CalculateMetrics.R')
source('../toSource/CensorData.R')
source('../toSource/CovFunc.R')
source('../toSource/DataExtraction.R')
source('../toSource/GenerateData.R')
source('../toSource/ImputeMissingData.R')
source('../toSource/LogPriorX.R')
source('../toSource/MakeSyntheticData.R')
source('../toSource/MeanFunc.R')
source('../toSource/NormaliseExpressionData.R')
source('../toSource/PlotHyperparamAnyKernel.R')
source('../toSource/PlotKaplanMeier.R')
source('../toSource/PlotPredictedMeasured.R')
source('../toSource/PlotTrainingTargetsLearned.R')
source('../toSource/PreLearnHyperparam.R')
source('../toSource/PrintOptions.R')
source('../toSource/PrintOptions2.R')
source('../toSource/RemoveCensored.R')


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
outerFolder 			<- 'Runs'
folderName 				<- paste0(outerFolder,'/',unid)
dir.create(file.path(getwd(),outerFolder),showWarnings=FALSE)

nReps 					<- 30

modelsList 				<- c('GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV','GPSurvInfUnif','GPSurvInfMed','AFT','Coxph','Coxnet','GBM','RF','RFSurvival')


##-------------------------------------------------------------------------------------##
##----------------------------- Plot and Save Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
printOptions 			<- FALSE
printPlots 				<- TRUE
savePlots 				<- TRUE
saveFiles 				<- TRUE
plotRuns 				<- TRUE
printResults 			<- TRUE
plotSaveOptions 		<- list('printOptions'=printOptions,'printPlots'=printPlots,'savePlots'=savePlots,'saveFiles'=saveFiles,'plotRuns'=plotRuns,
								'printResults'=printResults,'folderName'=folderName,'outerFolder'=outerFolder,'nReps'=nReps)


##-------------------------------------------------------------------------------------##
##---------------------------------- Data Parameters ----------------------------------##
##-------------------------------------------------------------------------------------##
dataSource 				<- 'Generate' 
dimension 				<- 5 
extraDimensions 		<- 0 
proportionTest 			<- NA
nTraining 				<- 400
nTest 					<- 50 

covFuncFormGen 			<- 'SqExp'
maternParamGen 			<- 3
meanFuncFormGen 		<- 'Zero'
if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE

if(useARD){
	logHypGenerate 		<- list('noise'=log(0.01),'func'=log(0.5),'length'=log(c(1.1,0.9,1.2,0.7)),'mean'=c(rep(0,dimension),0))
} else {
	logHypGenerate 		<- list('noise'=log(0.10),'func'=log(0.7),'length'=log(1.1),'mean'=c(rep(0,dimension),0))
}

gridMinimum 			<- 0
gridMaximum 			<- 8
censoringType 			<- 'NormalLoopSample'
censoringSD 			<- 50
censoringMean 			<- 0
censoredProportion 		<- 0.75
nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)

dataOptionsStructure 	<- list('dataSource'=dataSource,'logHypGenerate'=logHypGenerate,'covFuncFormGen'=covFuncFormGen,'meanFuncFormGen'=meanFuncFormGen,
								'maternParamGen'=maternParamGen,'dimension'=dimension,'nTraining'=nTraining,'nTest'=nTest,'gridMinimum'=gridMinimum,
								'gridMaximum'=gridMaximum,'censoringSD'=censoringSD,'censoringMean'=censoringMean,'censoringType'=censoringType,'nCensored'=nCensored,
								'useARD'=useARD,'extraDimensions'=extraDimensions,'folderName'=folderName,'proportionTest'=proportionTest,
								'censoredProportion'=censoredProportion)


##-------------------------------------------------------------------------------------##
##---------------------------- Gaussian Process Parameters ----------------------------##
##-------------------------------------------------------------------------------------##
modelType 				<- 'survival' 			# 'non-survival' or 'survival'
covFuncForm 			<- 'SqExp'
meanFuncForm 			<- 'Zero'
burnIn 					<- FALSE
maternParam 			<- 3
tolerance 				<- 0.0001*nTraining*censoredProportion
toleranceLaplace 		<- 0.0001*nTraining*censoredProportion
hypChangeTol 			<- 10^-10
maxCount 				<- 500
maxit 					<- 2000
maxitPreLearn 			<- 1
maxitSurvival 			<- 10
maxitLaplaceHyp			<- 2
maxitLaplaceFHat 		<- 100
optimType 				<- 'Nelder-Mead' 					# 'CG' = conjugate gradient optimisation, 'Nelder-Mead' = numerical gradients
noiseCorr 				<- FALSE
imposePriors 			<- TRUE
if(covFuncForm!='ARD'){
	logHypStart 	<- list('noise'=log(0.2)+rnorm(1,mean=0,sd=0.01),'func'=log(0.8)+rnorm(1,mean=0,sd=0.05),'length'=log(1.5)+rnorm(1,mean=0,sd=0.05),'mean'=c(rep(0,dimension),0))
} else {
	logHypStart 		<- list('noise'=log(0.2),'func'=log(0.8),'length'=rep(log(0.9),dimension),'mean'=c(rep(0,dimension),0))
}
parameterStructure 		<- list('meanFuncForm'=meanFuncForm,'covFuncForm'=covFuncForm,'maternParam'=maternParam,'maxit'=maxit,'maxitPreLearn'=maxitPreLearn,
								'maxitSurvival'=maxitSurvival,'optimType'=optimType,'logHypStart'=logHypStart,'modelType'=modelType,'unid'=unid,'tolerance'=tolerance,
								'toleranceLaplace'=toleranceLaplace,'maxCount'=maxCount,'hypChangeTol'=hypChangeTol,'burnIn'=burnIn,'maxitLaplaceHyp'=maxitLaplaceHyp,
								'maxitLaplaceFHat'=maxitLaplaceFHat,'noiseCorr'=noiseCorr,'imposePriors'=imposePriors)

if(printOptions) PrintOptions(dataOptionsStructure,parameterStructure,modelsList)
if(saveFiles){
	dir.create(file.path(paste0(getwd(),"/",outerFolder),unid),showWarnings=FALSE)
	if(dataSource=='Generate'){
		dir.create(file.path(paste0(getwd(),"/",outerFolder,'/',unid),'data'),showWarnings=FALSE)
		write.table(exp(unlist(logHypStart))[1:3],paste0(folderName,'/data/hypStart.csv'),row.names=FALSE,col.names=FALSE,sep=',')
		write.table(exp(unlist(logHypGenerate))[1:3],paste0(folderName,'/data/hypGenerate.csv'),row.names=FALSE,col.names=FALSE,sep=',')
	}
	sink(paste0(getwd(),'/',outerFolder,'/',unid,'/',unid,'RunOptions.txt'))
	PrintOptions2(dataOptionsStructure,parameterStructure)
	sink()
}


##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureGPNonSurvNoCens 	<- list()
outputStructureGPSurvNoCorr 	<- list()
outputStructureGPSurvCorrV		<- list()
outputStructureGPSurvCorrL 		<- list()
outputStructureGPSurvInfUnif 	<- list()
outputStructureGPSurvInfMed 	<- list()
outputStructureGPSurvLearnCens 	<- list()
outputStructureAFT				<- list()
outputStructureCoxph 			<- list()
outputStructureGBM				<- list()
outputStructureCoxnet 			<- list()
outputStructureRF				<- list()
outputStructureRFSurvival 		<- list()

c.index.GPNonSurvNoCens 		<- rep(NA,nReps)
c.index.GPSurvNoCorr 			<- rep(NA,nReps)
c.index.GPSurvCorrV				<- rep(NA,nReps)
c.index.GPSurvCorrL 			<- rep(NA,nReps)
c.index.GPSurvInfUnif 			<- rep(NA,nReps)
c.index.GPSurvInfMed 			<- rep(NA,nReps)
c.index.GPSurvLearnCens 		<- rep(NA,nReps)
c.index.AFT						<- rep(NA,nReps)
c.index.Coxph 					<- rep(NA,nReps)
c.index.GBM						<- rep(NA,nReps)
c.index.Coxnet 					<- rep(NA,nReps)
c.index.RF						<- rep(NA,nReps)
c.index.RFSurvival 				<- rep(NA,nReps)

rmse.GPNonSurvNoCens 			<- rep(NA,nReps)
rmse.GPSurvNoCorr 				<- rep(NA,nReps)
rmse.GPSurvCorrV				<- rep(NA,nReps)
rmse.GPSurvCorrL 				<- rep(NA,nReps)
rmse.GPSurvInfUnif 				<- rep(NA,nReps)
rmse.GPSurvInfMed 				<- rep(NA,nReps)
rmse.GPSurvLearnCens 			<- rep(NA,nReps)
rmse.AFT						<- rep(NA,nReps)
rmse.Coxph 						<- rep(NA,nReps)
rmse.GBM						<- rep(NA,nReps)
rmse.Coxnet 					<- rep(NA,nReps)
rmse.RF							<- rep(NA,nReps)
rmse.RFSurvival 				<- rep(NA,nReps)

brier.GPNonSurvNoCens 			<- rep(NA,nReps)
brier.GPSurvNoCorr 				<- rep(NA,nReps)
brier.GPSurvCorrV				<- rep(NA,nReps)
brier.GPSurvCorrL 				<- rep(NA,nReps)
brier.GPSurvInfUnif 			<- rep(NA,nReps)
brier.GPSurvInfMed 				<- rep(NA,nReps)
brier.GPSurvLearnCens 			<- rep(NA,nReps)
brier.AFT						<- rep(NA,nReps)
brier.Coxph 					<- rep(NA,nReps)
brier.GBM						<- rep(NA,nReps)
brier.Coxnet 					<- rep(NA,nReps)
brier.RF						<- rep(NA,nReps)
brier.RFSurvival 				<- rep(NA,nReps)


##-------------------------------------------------------------------------------------##
##----------------------------------- Generate Data -----------------------------------##
##-------------------------------------------------------------------------------------##
trainingTestStructureForNReps 			<- GenerateData(dataOptionsStructure,outerFolder,nReps)
for(i in 1:nReps){
	trainingTestStructureForNReps[[i]] 	<- NormaliseExpressionData(trainingTestStructureForNReps[[i]],normaliseFlag=TRUE,winsoriseFlag=FALSE) # normaliseFlag TRUE -> genes normalised to (mean=0,sd=1), winsoriseFlag TRUE -> outside (0.05,0.95) quantiles clipped to here
	cat('Data set for run',i,'normalised',fill=TRUE)
}


##-------------------------------------------------------------------------------------##
##----------------------------- Run GPNonSurvNoCens Model -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- 'None'
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'non-survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	outputStructureGPNonSurvNoCens[[i]] <- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$c.index)!=0,outputStructureGPNonSurvNoCens[[i]]$c.index,NA)
	rmse.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$rmse)!=0,outputStructureGPNonSurvNoCens[[i]]$rmse,NA)
	brier.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$brier)!=0,outputStructureGPNonSurvNoCens[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPS1 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvNoCorr[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvNoCorr[i] 			<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$c.index)!=0,outputStructureGPSurvNoCorr[[i]]$c.index,NA)
	rmse.GPSurvNoCorr[i] 				<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$rmse)!=0,outputStructureGPSurvNoCorr[[i]]$rmse,NA)
	brier.GPSurvNoCorr[i] 				<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$brier)!=0,outputStructureGPSurvNoCorr[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPS2 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'noiseCorrLearned'
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvCorrL[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$c.index)!=0,outputStructureGPSurvCorrL[[i]]$c.index,NA)
	rmse.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$rmse)!=0,outputStructureGPSurvCorrL[[i]]$rmse,NA)
	brier.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$brier)!=0,outputStructureGPSurvCorrL[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run GPS3 Model -----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvCorrV[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$c.index)!=0,outputStructureGPSurvCorrV[[i]]$c.index,NA)
	rmse.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$rmse)!=0,outputStructureGPSurvCorrV[[i]]$rmse,NA)
	brier.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$brier)!=0,outputStructureGPSurvCorrV[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##------------------------------- Run GPR1 InfUnif Model ------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'none'
parameterStructure$modelType 			<- 'survival'
parameterStructure$inferenceType 		<- 'uniform'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInfUnif[[i]] 	<- ApplyGPDiffInf(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInfUnif[i] 			<- ifelse(length(outputStructureGPSurvInfUnif[[i]]$c.index)!=0,outputStructureGPSurvInfUnif[[i]]$c.index,NA)
	rmse.GPSurvInfUnif[i] 				<- ifelse(length(outputStructureGPSurvInfUnif[[i]]$rmse)!=0,outputStructureGPSurvInfUnif[[i]]$rmse,NA)
	brier.GPSurvInfUnif[i] 				<- ifelse(length(outputStructureGPSurvInfUnif[[i]]$brier)!=0,outputStructureGPSurvInfUnif[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##------------------------------- Run GPR2 InfMed Model -------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'none'
parameterStructure$modelType 			<- 'survival'
parameterStructure$inferenceType 		<- 'median'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInfMed[[i]] 	<- ApplyGPDiffInf(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInfMed[i] 			<- ifelse(length(outputStructureGPSurvInfMed[[i]]$c.index)!=0,outputStructureGPSurvInfMed[[i]]$c.index,NA)
	rmse.GPSurvInfMed[i] 				<- ifelse(length(outputStructureGPSurvInfMed[[i]]$rmse)!=0,outputStructureGPSurvInfMed[[i]]$rmse,NA)
	brier.GPSurvInfMed[i] 				<- ifelse(length(outputStructureGPSurvInfMed[[i]]$brier)!=0,outputStructureGPSurvInfMed[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##------------------------------ Run GPR3 LearnCens Model -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- 'none'
parameterStructure$modelType 			<- 'survival'
parameterStructure$inferenceType 		<- 'learn'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvLearnCens[[i]] <- ApplyGPLearnCens(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvLearnCens[i] 			<- ifelse(length(outputStructureGPSurvLearnCens[[i]]$c.index)!=0,outputStructureGPSurvLearnCens[[i]]$c.index,NA)
	rmse.GPSurvLearnCens[i] 			<- ifelse(length(outputStructureGPSurvLearnCens[[i]]$rmse)!=0,outputStructureGPSurvLearnCens[[i]]$rmse,NA)
	brier.GPSurvLearnCens[i] 			<- ifelse(length(outputStructureGPSurvLearnCens[[i]]$brier)!=0,outputStructureGPSurvLearnCens[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##----------------------------------- Run AFT Model -----------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureAFT[[i]] 			<- ApplyAFT(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.AFT[i] 						<- ifelse(length(outputStructureAFT[[i]]$c.index)!=0,outputStructureAFT[[i]]$c.index,NA)
	rmse.AFT[i] 						<- ifelse(length(outputStructureAFT[[i]]$rmse)!=0,outputStructureAFT[[i]]$rmse,NA)
	brier.AFT[i] 						<- ifelse(length(outputStructureAFT[[i]]$brier)!=0,outputStructureAFT[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run Coxph Model ----------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureCoxph[[i]] 			<- ApplyCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Coxph[i] 					<- ifelse(length(outputStructureCoxph[[i]]$c.index)!=0,outputStructureCoxph[[i]]$c.index,NA)
	rmse.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$rmse)!=0,outputStructureCoxph[[i]]$rmse,NA)
	brier.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$brier)!=0,outputStructureCoxph[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##----------------------------------- Run GBM Model -----------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGBM[[i]] 			<- ApplyGBM(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GBM[i] 						<- ifelse(length(outputStructureGBM[[i]]$c.index)!=0,outputStructureGBM[[i]]$c.index,NA)
	rmse.GBM[i] 						<- ifelse(length(outputStructureGBM[[i]]$rmse)!=0,outputStructureGBM[[i]]$rmse,NA)
	brier.GBM[i] 						<- ifelse(length(outputStructureGBM[[i]]$brier)!=0,outputStructureGBM[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##---------------------------------- Run Coxnet Model ---------------------------------##
##-------------------------------------------------------------------------------------##
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureCoxnet[[i]] 			<- ApplyCoxnet(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Coxnet[i] 					<- ifelse(length(outputStructureCoxnet[[i]]$c.index)!=0,outputStructureCoxnet[[i]]$c.index,NA)
	rmse.Coxnet[i] 						<- ifelse(length(outputStructureCoxnet[[i]]$rmse)!=0,outputStructureCoxnet[[i]]$rmse,NA)
	brier.Coxnet[i] 					<- ifelse(length(outputStructureCoxnet[[i]]$brier)!=0,outputStructureCoxnet[[i]]$brier,NA)
}


#-------------------------------------------------------------------------------------##
#------------------------------------ Run RF Model -----------------------------------##
#-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- 'None'
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'non-survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	outputStructureRF[[i]] 				<- ApplyRF(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.RF[i] 						<- ifelse(length(outputStructureRF[[i]]$c.index)!=0,outputStructureRF[[i]]$c.index,NA)
	rmse.RF[i] 							<- ifelse(length(outputStructureRF[[i]]$rmse)!=0,outputStructureRF[[i]]$rmse,NA)
	brier.RF[i] 						<- ifelse(length(outputStructureRF[[i]]$brier)!=0,outputStructureRF[[i]]$brier,NA)
}


#-------------------------------------------------------------------------------------##
#------------------------------- Run RFSurvival Model --------------------------------##
#-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureRFSurvival[[i]] 		<- ApplyRFSurvival(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.RFSurvival[i] 				<- ifelse(length(outputStructureRFSurvival[[i]]$c.index)!=0,outputStructureRFSurvival[[i]]$c.index,NA)
	rmse.RFSurvival[i] 					<- ifelse(length(outputStructureRFSurvival[[i]]$rmse)!=0,outputStructureRFSurvival[[i]]$rmse,NA)
	brier.RFSurvival[i] 				<- ifelse(length(outputStructureRFSurvival[[i]]$brier)!=0,outputStructureRFSurvival[[i]]$brier,NA)
}


##-------------------------------------------------------------------------------------##
##------------------------------ Print and Save Results -------------------------------##
##-------------------------------------------------------------------------------------##
							##-------------------------------##
							##---- Print C Index & RMSE -----##
							##-------------------------------##
if(printResults){
	cat('---------------------------------------',fill=TRUE)
	cat('Concordance Indices:',fill=TRUE)
	cat('GPNonSurvNoCens mean c index =',paste0(round(c.index.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS1 mean c index =',paste0(round(c.index.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPS2 mean c index =',paste0(round(c.index.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('GPS3 mean c index =',paste0(round(c.index.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('AFT mean c index =',paste0(round(c.index.AFT,4),collapse=', '),fill=TRUE)
	cat('Coxph mean c index =',paste0(round(c.index.Coxph,4),collapse=', '),fill=TRUE)
	cat('GBM mean c index =',paste0(round(c.index.GBM,4),collapse=', '),fill=TRUE)
	cat('Coxnet mean c index =',paste0(round(c.index.Coxnet,4),collapse=', '),fill=TRUE)
	cat('RF mean c index =',paste0(round(c.index.RF,4),collapse=', '),fill=TRUE)
	cat('RFSurvival mean c index =',paste0(round(c.index.RFSurvival,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
	cat('RMSE:',fill=TRUE)
	cat('GP mean rmse =',paste0(round(rmse.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS1 mean rmse =',paste0(round(rmse.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPS2 mean rmse =',paste0(round(rmse.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('GPS3 mean rmse =',paste0(round(rmse.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
}

							##-------------------------------##
							##--- Plot Kaplan-Meier Plots ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GPNonSurv','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPNonSurvNoCens[[i]]$funcMeanPred,outputStructureGPNonSurvNoCens[[i]]$trainingTestStructure$testTargets,outputStructureGPNonSurvNoCens[[i]]$trainingTestStructure$testEvents,'GPNonSurv')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvNoCorr[[i]]$funcMeanPred,outputStructureGPSurvNoCorr[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvNoCorr[[i]]$trainingTestStructure$testEvents,'GPSurvNoCorr')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvCorrL[[i]]$funcMeanPred,outputStructureGPSurvCorrL[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvCorrL[[i]]$trainingTestStructure$testEvents,'GPSurvCorrL')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvCorrV[[i]]$funcMeanPred,outputStructureGPSurvCorrV[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvCorrV[[i]]$trainingTestStructure$testEvents,'GPSurvCorrV')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvInfUnif[[i]]$funcMeanPred,outputStructureGPSurvInfUnif[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvInfUnif[[i]]$trainingTestStructure$testEvents,'GPSurvInfUnif')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvInfMed[[i]]$funcMeanPred,outputStructureGPSurvInfMed[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvInfMed[[i]]$trainingTestStructure$testEvents,'GPSurvInfMed')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfLearn','/',unid,'GPSurvInfLearn','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvLearnCens[[i]]$funcMeanPred,outputStructureGPSurvLearnCens[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvLearnCens[[i]]$trainingTestStructure$testEvents,'GPSurvNoCorr')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','AFT','/',unid,'AFT','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureAFT[[i]]$testPredictions,outputStructureAFT[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureAFT[[i]]$trainingTestStructure$testEvents,'AFT')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','Coxph','/',unid,'Coxph','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureCoxph[[i]]$testPredictions,outputStructureCoxph[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureCoxph[[i]]$trainingTestStructure$testEvents,'Coxph')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GBM','/',unid,'GBM','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGBM[[i]]$testPredictions,outputStructureGBM[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGBM[[i]]$trainingTestStructure$testEvents,'GMB')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','Coxnet','/',unid,'Coxnet','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureCoxnet[[i]]$testPredictions,outputStructureCoxnet[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureCoxnet[[i]]$trainingTestStructure$testEvents,'Coxnet')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','RF','/',unid,'RF','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureRF[[i]]$testPredictions,outputStructureRF[[i]]$trainingTestStructure$testTargets,outputStructureRF[[i]]$trainingTestStructure$testEvents,'RF')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','RFSurvival','/',unid,'RFSurvival','PlotMeasuredKM.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureRFSurvival[[i]]$testPredictions,outputStructureRFSurvival[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureRFSurvival[[i]]$trainingTestStructure$testEvents,'RFSurvival')
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							##--- Plot Measured/Predicted ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GPNonSurv','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPNonSurvNoCens[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPSurvNoCorr[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPSurvCorrL[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPSurvCorrV[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPSurvInfUnif[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPSurvInfMed[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfLearn','/',unid,'GPSurvInfLearn','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGPSurvLearnCens[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','AFT','/',unid,'AFT','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureAFT[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','Coxph','/',unid,'Coxph','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureCoxph[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GBM','/',unid,'GBM','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureGBM[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','Coxnet','/',unid,'Coxnet','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureCoxnet[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','RF','/',unid,'RF','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureRF[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','RFSurvival','/',unid,'RFSurvival','PlotMeasuredPredicted.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotPredictedMeasured(outputStructureRFSurvival[[i]])
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							## Plot Training Targets Learned ##
							##-------------------------------##							
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotTrainingTargetsLearned.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotTrainingTargetsLearned(outputStructureGPSurvNoCorr[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotTrainingTargetsLearned.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotTrainingTargetsLearned(outputStructureGPSurvCorrL[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotTrainingTargetsLearned.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotTrainingTargetsLearned(outputStructureGPSurvCorrV[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotTrainingTargetsLearned.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotTrainingTargetsLearned(outputStructureGPSurvInfUnif[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotTrainingTargetsLearned.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotTrainingTargetsLearned(outputStructureGPSurvInfMed[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfLearn','/',unid,'GPSurvInfLearn','PlotTrainingTargetsLearned.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotTrainingTargetsLearned(outputStructureGPSurvLearnCens[[i]])
		}
	dev.off(pdf.output)
}

							##-------------------------------##
							##--- Plot GP Hyperparameters ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GPNonSurv','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPNonSurvNoCens[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvNoCorr[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvCorrL[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvCorrV[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfUnif','/',unid,'GPSurvInfUnif','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvInfUnif[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfMed','/',unid,'GPSurvInfMed','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvInfMed[[i]])
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvInfLearn','/',unid,'GPSurvInfLearn','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvLearnCens[[i]])
		}
	dev.off(pdf.output)
}


							##--------------------------------##
							##--- Plot Concordance Indices ---##
							##--------------------------------##
modelNames 		<- c('AFT','Coxph','Coxnet','GBM','RFSurvival','GPNonSurvNoCens','GPSurvInfMed','GPSurvInfUnif','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV') # Reordered version of modelsList
modelNamesPlot 	<- c('AFT','Cox PH','Coxnet','GBM','RSF','GP','GPR1','GPR2','GPS1','GPS2','GPS3')
toInvert 		<- c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
c.index.mean 	<- numeric()
c.index.mat 	<- matrix(0,nrow=nReps,ncol=length(modelNames))
brier.mean 		<- numeric()
brier.mat 		<- matrix(0,nrow=nReps,ncol=length(modelNames))
rmse.mean 		<- numeric()
rmse.mat 		<- matrix(0,nrow=nReps,ncol=length(modelNames))
for(i in 1:length(modelNames)){
	c.index.mat[,i] 	<- get(paste0('c.index.',modelNames[i]))
	if(toInvert[i]) c.index.mat[,i] <- 1-c.index.mat[,i]
	c.index.mean[i] 	<- mean(c.index.mat[,i],na.rm=TRUE)

	brier.mat[,i] 		<- get(paste0('brier.',modelNames[i]))
	brier.mean[i] 		<- mean(brier.mat[,i],na.rm=TRUE)

	rmse.mat[,i] 		<- get(paste0('rmse.',modelNames[i]))
	rmse.mean[i] 		<- mean(rmse.mat[,i],na.rm=TRUE)
}
cols 			<- c(rep('green4',5),rep('chartreuse4',3),rep('chartreuse3',3))

if(savePlots){
	pdf(file=paste0(getwd(),"/",outerFolder,"/",unid,'/',unid,'PlotCIndexAllModelsBoxplot.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(c.index.mat,col=add.alpha(cols,0.8),xaxt='n',xlab='',ylab='Concordance Index',ylim=c(0.4,1),las=2)
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}

if(savePlots){
	pdf(file=paste0(getwd(),"/",outerFolder,"/",unid,'/',unid,'PlotBrierAllModelsBoxplot.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(brier.mat,col=add.alpha('chartreuse3',0.8),xaxt='n',xlab='',ylab='Integrated Brier Score')
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}

if(savePlots){
	pdf(file=paste0(getwd(),"/",outerFolder,"/",unid,'/',unid,'PlotRMSEAllModelsBoxplot.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(rmse.mat,col=add.alpha('chartreuse3',0.8),xaxt='n',xlab='',ylab='RMSE')
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}


##-------------------------------------------------------------------------------------##
##------------------------------------ Save Output ------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureAll <- list('outputStructureGPNonSurvNoCens'=outputStructureGPNonSurvNoCens,
							'outputStructureGPSurvNoCorr'=outputStructureGPSurvNoCorr,
							'outputStructureGPSurvCorrV'=outputStructureGPSurvCorrV,
							'outputStructureGPSurvCorrL'=outputStructureGPSurvCorrL,
							'outputStructureGPSurvLaplace'=outputStructureGPSurvLaplace,
							'outputStructureGPSurvInfUnif'=outputStructureGPSurvInfUnif,
							'outputStructureGPSurvInfMed'=outputStructureGPSurvInfMed,
							'outputStructureGPSurvLearnCens'=outputStructureGPSurvLearnCens,
							'outputStructureAFT'=outputStructureAFT,
							'outputStructureCoxph'=outputStructureCoxph,
							'outputStructureCoxnet'=outputStructureCoxnet,
							'outputStructureGBM'=outputStructureGBM,
							'outputStructureRF'=outputStructureRF,
							'outputStructureRFSurvival'=outputStructureRFSurvival)

save(list='outputStructureAll',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureAll','_','Workspace.RData'))