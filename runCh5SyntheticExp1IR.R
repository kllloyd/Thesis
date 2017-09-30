#----------------------------------------------------------------------------------------#
# K Lloyd 2017_08_09
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# Models investigating feature selection
# All models applied to the same data. Models are:	GP, GPS3SqExp, GPS3ARD, GPS3IARD1, GPS3IARD2, GPS3SqExpRSFS,
# 													GPSurvBICForward, GPSurvBICBackward,
# 													Cox PH, Coxnet, StepCox, FilterCoxph1, FilterCoxph2,
# 													RF and RSF
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Compare c index results of different models
#---------------------------------------------------------------------------------------#

##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##
# args = commandArgs(trailingOnly=TRUE)

library(devtools)
library(fields)
library(gbm)
library(glmnet)
library(impute)
library(ipred)
library(MASS)
library(Matrix)
library(mice)
library(nlme)
library(NORMT3)
library(pdist)
library(randomForestSRC)
# library(rgl)
# library(rms)
library(survcomp)
library(survival)
library(zoo)


##-------------------------------------------------------------------------------------##
##--------------------------------- Source Functions ----------------------------------##
##-------------------------------------------------------------------------------------##
source('../toSource/add.alpha.R')
source('../toSource/AddFeature.R')
source('../toSource/AdjustTrainingSurvivalMeanVariance.R')
source('../toSource/ApplyAFT.R')
source('../toSource/ApplyCoxph.R')
source('../toSource/ApplyFilterCoxph.R')
source('../toSource/ApplyGBM.R')
source('../toSource/ApplyCoxnet.R')
source('../toSource/ApplyGP.R')
source('../toSource/ApplyGPDiffInf.R')
source('../toSource/ApplyGPLaplace.R')
source('../toSource/ApplyGPSurvBIC.R')
source('../toSource/ApplyRF.R')
source('../toSource/ApplyRFSurvival.R')
source('../toSource/ApplyStepCoxph.R')
source('../toSource/CalculateBIC.R')
source('../toSource/CalculateEnsembleResults.R')
source('../toSource/CalculateMetrics.R')
source('../toSource/CensorData.R')
source('../toSource/ChooseSubsetsRandomly.R')
source('../toSource/CovFunc.R')
source('../toSource/DataExtraction.R')
source('../toSource/GenerateData.R')
source('../toSource/ImputeMissingData.R')
source('../toSource/LogPriorX.R')
source('../toSource/MakeEnsemblePlots.R')
source('../toSource/MakeSyntheticData.R')
source('../toSource/MeanFunc.R')
source('../toSource/ModifyFeatureListsInformedARD2.R')
source('../toSource/NormaliseExpressionData.R')
source('../toSource/PlotHyperparamAnyKernel.R')
source('../toSource/PlotKaplanMeier.R')
source('../toSource/PostProcessingHyperparameterBoxplots.R')
source('../toSource/PreLearnHyperparam.R')
source('../toSource/PrintOptions.R')
source('../toSource/PrintOptions2.R')
source('../toSource/RemoveCensored.R')
source('../toSource/RemoveFeature.R')
source('../toSource/RunModelCombinationsBackward.R')
source('../toSource/RunModelCombinationsForward.R')
source('../toSource/SubsetIndices.R')

##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
# k 						<- as.numeric(args[1])
# unid 					<- paste0(format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S'),'k',k)
unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
outerFolder 			<- 'Runs'
folderName 				<- paste0(outerFolder,'/',unid)
dir.create(file.path(getwd(),outerFolder),showWarnings=FALSE)

nReps 					<- 50

modelsList 				<- c('GPNonSurvNoCens','GPSurvSqExp','GPSurvARD','GPSurvInformedARD2','GPSurvInformedARDV3','GPSurvInformedARDV3','GPSurvInformedARDV4','GPSurvInfMed','GPSurvInfUnif','GPSurvBICForward','GPSurvBICBackward','RF','RFSurvival','Coxph','Coxnet','StepCoxph')


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
dimension 				<- 20
extraDimensions 		<- 5 
proportionTest 			<- NA
nTraining 				<- 500
nTest 					<- 100 

covFuncFormGen 			<- 'ARD'
extraParamGen 			<- 3
meanFuncFormGen 		<- 'Linear'
if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE

if(useARD){
	# 20 D
	# groups: 1,2,3,4,5,6 ; 4,5,6,7,8,9,10,11,12 ; 11,12,13,14,15, ; 16,17,18,19,20
	logHypGenerate 	<- list('noise'=log(0.05),'func'=log(0.7),'length'=log(c(0.9,1.1,1.3,1.6,1.65,1.7,1.9,1.95,2.0,2.2,2.6,2.7,3.5,3.7,4)),'mean'=c(rep(0,dimension),0))
} else {
	logHypGenerate 		<- list('noise'=log(0.01),'func'=log(0.5),'length'=log(1.1),'mean'=c(rep(0,dimension),0))
}

gridMinimum 			<- 0
gridMaximum 			<- 8
censoringType 			<- 'NormalLoopSample'
censoringSD 			<- 50
censoringMean 			<- 0
censoredProportion 		<- 0.70
nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)

dataOptionsStructure 	<- list('dataSource'=dataSource,'logHypGenerate'=logHypGenerate,'covFuncFormGen'=covFuncFormGen,'meanFuncFormGen'=meanFuncFormGen,
								'extraParamGen'=extraParamGen,'dimension'=dimension,'nTraining'=nTraining,'nTest'=nTest,'gridMinimum'=gridMinimum,
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
extraParam 				<- 3
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
	logHypStart 		<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(2),'mean'=c(rep(0,dimension),0))
} else {
	logHypStart 		<- list('noise'=log(0.2),'func'=log(0.8),'length'=rep(log(2),dimension),'mean'=c(rep(0,dimension),0))

}
parameterStructure 		<- list('meanFuncForm'=meanFuncForm,'covFuncForm'=covFuncForm,'extraParam'=extraParam,'maxit'=maxit,'maxitPreLearn'=maxitPreLearn,
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
outputStructureGPNonSurvNoCens 		<- list()
outputStructureGPSurvSqExp 			<- list()
outputStructureGPSurvARD			<- list()
outputStructureGPSurvInformedARD1 	<- list()
outputStructureGPSurvInformedARD2 	<- list()
outputStructureGPSurvBICForward 	<- list()
outputStructureGPSurvBICBackward 	<- list()
outputStructureRF					<- list()
outputStructureRFSurvival 			<- list()
outputStructureCoxph 				<- list()
outputStructureFilterCoxphU 		<- list()
outputStructureFilterCoxphC 		<- list()
outputStructureCoxnet 				<- list()
outputStructureStepCoxph 			<- list()
c.index.GPNonSurvNoCens 			<- rep(NA,nReps)
c.index.GPSurvSqExp 				<- rep(NA,nReps)
c.index.GPSurvARD					<- rep(NA,nReps)
c.index.GPSurvInformedARD1 			<- rep(NA,nReps)
c.index.GPSurvInformedARD2 			<- rep(NA,nReps)
c.index.GPSurvBICForward 			<- rep(NA,nReps)
c.index.GPSurvBICBackward 			<- rep(NA,nReps)
c.index.RF							<- rep(NA,nReps)
c.index.RFSurvival 					<- rep(NA,nReps)
c.index.Coxph 						<- rep(NA,nReps)
c.index.FilterCoxphU 				<- rep(NA,nReps)
c.index.FilterCoxphC 				<- rep(NA,nReps)
c.index.Coxnet 						<- rep(NA,nReps)
c.index.StepCoxph 					<- rep(NA,nReps)
rmse.GPNonSurvNoCens 				<- rep(NA,nReps)
rmse.GPSurvSqExp 					<- rep(NA,nReps)
rmse.GPSurvARD						<- rep(NA,nReps)
rmse.GPSurvInformedARD1 			<- rep(NA,nReps)
rmse.GPSurvInformedARD2 			<- rep(NA,nReps)
rmse.GPSurvBICForward 				<- rep(NA,nReps)
rmse.GPSurvBICBackward 				<- rep(NA,nReps)
rmse.RF								<- rep(NA,nReps)
rmse.RFSurvival 					<- rep(NA,nReps)
rmse.Coxph 							<- rep(NA,nReps)
rmse.FilterCoxphU 					<- rep(NA,nReps)
rmse.FilterCoxphC 					<- rep(NA,nReps)
rmse.Coxnet 						<- rep(NA,nReps)
rmse.StepCoxph 						<- rep(NA,nReps)
brier.GPNonSurvNoCens 				<- rep(NA,nReps)
brier.GPSurvSqExp 					<- rep(NA,nReps)
brier.GPSurvARD						<- rep(NA,nReps)
brier.GPSurvInformedARD1 			<- rep(NA,nReps)
brier.GPSurvInformedARD2 			<- rep(NA,nReps)
brier.GPSurvBICForward 				<- rep(NA,nReps)
brier.GPSurvBICBackward 			<- rep(NA,nReps)
brier.RF							<- rep(NA,nReps)
brier.RFSurvival 					<- rep(NA,nReps)
brier.Coxph 						<- rep(NA,nReps)
brier.FilterCoxphU 					<- rep(NA,nReps)
brier.FilterCoxphC 					<- rep(NA,nReps)
brier.Coxnet 						<- rep(NA,nReps)
brier.StepCoxph 					<- rep(NA,nReps)

##-------------------------------------------------------------------------------------##
##----------------------------------- Generate Data -----------------------------------##
##-------------------------------------------------------------------------------------##
trainingTestStructureForNReps 			<- GenerateData(dataOptionsStructure,outerFolder,nReps)
for(i in 1:nReps){
	trainingTestStructureForNReps[[i]] 	<- NormaliseExpressionData(trainingTestStructureForNReps[[i]],normaliseFlag=TRUE,winsoriseFlag=FALSE) # normaliseFlag TRUE -> genes normalised to (mean=0,sd=1), winsoriseFlag TRUE -> outside (0.05,0.95) quantiles clipped to here
	cat('Data set for run',i,'normalised',fill=TRUE)
	colnames(trainingTestStructureForNReps[[i]]$trainingData) 	<- paste0('x',1:dataOptionsStructure$dimension)
	colnames(trainingTestStructureForNReps[[i]]$testData) 		<- paste0('x',1:dataOptionsStructure$dimension)
}
save('trainingTestStructureForNReps',file=paste0('Runs','/','trainingTestStructureForNReps','_','Workspace_',dimension,'D.RData'))


##-------------------------------------------------------------------------------------##
##----------------------------- Run GPNonSurvNoCens Model -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- 'None'
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'non-survival'
parameterStructure$covFuncForm 			<- 'SqExp'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dimension),0))
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	outputStructureGPNonSurvNoCens[[i]] <- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$c.index)!=0,outputStructureGPNonSurvNoCens[[i]]$c.index,NA)
	rmse.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$rmse)!=0,outputStructureGPNonSurvNoCens[[i]]$rmse,NA)
	brier.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$brier)!=0,outputStructureGPNonSurvNoCens[[i]]$brier,NA)
}
save(list='outputStructureGPNonSurvNoCens',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPNonSurvNoCens','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##------------------------------- Run GPS3 SqExp Model --------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$covFuncForm 			<- 'SqExp'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dimension),0))
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
parameterStructure$tolerance 			<- tolerance
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvSqExp[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvSqExp[i] 				<- ifelse(length(outputStructureGPSurvSqExp[[i]]$c.index)!=0,outputStructureGPSurvSqExp[[i]]$c.index,NA)
	rmse.GPSurvSqExp[i] 				<- ifelse(length(outputStructureGPSurvSqExp[[i]]$rmse)!=0,outputStructureGPSurvSqExp[[i]]$rmse,NA)
	brier.GPSurvSqExp[i] 				<- ifelse(length(outputStructureGPSurvSqExp[[i]]$brier)!=0,outputStructureGPSurvSqExp[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExp'))
save(list='outputStructureGPSurvSqExp',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExp','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##-------------------------------- Run GPS3 ARD Model ---------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$covFuncForm 			<- 'ARD'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=rep(log(0.9),dimension),'mean'=c(rep(0,dimension),0))
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
parameterStructure$tolerance 			<- tolerance
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvARD[[i]] 		<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvARD[i] 				<- ifelse(length(outputStructureGPSurvARD[[i]]$c.index)!=0,outputStructureGPSurvARD[[i]]$c.index,NA)
	rmse.GPSurvARD[i] 					<- ifelse(length(outputStructureGPSurvARD[[i]]$rmse)!=0,outputStructureGPSurvARD[[i]]$rmse,NA)
	brier.GPSurvARD[i] 					<- ifelse(length(outputStructureGPSurvARD[[i]]$brier)!=0,outputStructureGPSurvARD[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvARD'))
save(list='outputStructureGPSurvARD',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvARD','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##---------------------------- Run GPS3 InformedARD1 Model ----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$covFuncForm 			<- 'InformedARD'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(c(1.1,2,3,40)),'mean'=c(rep(0,dimension),0))
parameterStructure$extraParam 			<- list('list1'=c('x1','x2','x3','x4','x5','x6'),
												'list2'=c('x4','x5','x6','x7','x8','x9','x10','x11','x12'),
												'list3'=c('x11','x12','x13','x14','x15'),
												'list4'=c('x16','x17','x18','x19','x20'))	
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
parameterStructure$tolerance 			<- tolerance
for(i in 1:nReps){
	trainingTestStructure 						<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInformedARD1[[i]] 		<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInformedARD1[i] 				<- ifelse(length(outputStructureGPSurvInformedARD1[[i]]$c.index)!=0,outputStructureGPSurvInformedARD1[[i]]$c.index,NA)
	rmse.GPSurvInformedARD1[i] 					<- ifelse(length(outputStructureGPSurvInformedARD1[[i]]$rmse)!=0,outputStructureGPSurvInformedARD1[[i]]$rmse,NA)
	brier.GPSurvInformedARD1[i] 				<- ifelse(length(outputStructureGPSurvInformedARD1[[i]]$brier)!=0,outputStructureGPSurvInformedARD1[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvInformedARD1'))
save(list='outputStructureGPSurvInformedARD1',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvInformedARD1','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##--------------------------- Run GPS3 InformedARD2 Model ----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$covFuncForm 			<- 'InformedARD'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(c(1.1,2,3,40)),'mean'=c(rep(0,dimension),0))
parameterStructure$extraParam 			<- list('list1'=c('x1','x2','x3','x4','x5','x6'),
												'list2'=c('x4','x5','x6','x7','x8','x9','x10','x11','x12'),
												'list3'=c('x11','x12','x13','x14','x15'),
												'list4'=c('x16','x17','x18','x19','x20'))
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
parameterStructure$tolerance 			<- tolerance
parameterStructure 						<- ModifyFeatureListsInformedARD2(dataOptionsStructure,parameterStructure,trainingTestStructureForNReps[[1]])
for(i in 1:nReps){
	trainingTestStructure 					<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInformedARD2[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInformedARD2[i] 			<- ifelse(length(outputStructureGPSurvInformedARD2[[i]]$c.index)!=0,outputStructureGPSurvInformedARD2[[i]]$c.index,NA)
	rmse.GPSurvInformedARD2[i] 				<- ifelse(length(outputStructureGPSurvInformedARD2[[i]]$rmse)!=0,outputStructureGPSurvInformedARD2[[i]]$rmse,NA)
	brier.GPSurvInformedARD2[i] 			<- ifelse(length(outputStructureGPSurvInformedARD2[[i]]$brier)!=0,outputStructureGPSurvInformedARD2[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvInformedARD2'))
save(list='outputStructureGPSurvInformedARD2',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvInformedARD2','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##--------------------------- Run GPS3 BIC Forward Selection --------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 			<- censoringType
parameterStructure$covFuncForm 				<- 'SqExp'
parameterStructure$logHypStart 				<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dimension),0))
parameterStructure$noiseCorr 				<- 'noiseCorrVec'
parameterStructure$modelType 				<- 'survival'
parameterStructure$tolerance 				<- tolerance
parameterStructure$direction 				<- 'forward'
allParameterStructures 						<- list('dataOptionsStructure'=dataOptionsStructure,'parameterStructure'=parameterStructure,'plotSaveOptions'=plotSaveOptions)
for(i in 1:nReps){
	trainingTestStructure 					<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvBICForward[[i]] 	<- ApplyGPSurvBIC(trainingTestStructure,allParameterStructures)
	c.index.GPSurvBICForward[i] 			<- ifelse(length(outputStructureGPSurvBICForward[[i]]$c.index)!=0,outputStructureGPSurvBICForward[[i]]$c.index,NA)
	rmse.GPSurvBICForward[i] 				<- ifelse(length(outputStructureGPSurvBICForward[[i]]$rmse)!=0,outputStructureGPSurvBICForward[[i]]$rmse,NA)
	brier.GPSurvBICForward[i] 				<- ifelse(length(outputStructureGPSurvBICForward[[i]]$brier)!=0,outputStructureGPSurvBICForward[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvBICForward'))
save(list='outputStructureGPSurvBICForward',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvBICForward','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##------------------------- Run GPS3 BIC Backward Elimination -------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 			<- censoringType
parameterStructure$covFuncForm 				<- 'SqExp'
parameterStructure$logHypStart 				<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dimension),0))
parameterStructure$noiseCorr 				<- 'noiseCorrVec'
parameterStructure$modelType 				<- 'survival'
parameterStructure$tolerance 				<- tolerance
parameterStructure$direction 				<- 'backward'
allParameterStructures 						<- list('dataOptionsStructure'=dataOptionsStructure,'parameterStructure'=parameterStructure,'plotSaveOptions'=plotSaveOptions)
for(i in 1:nReps){
	trainingTestStructure 					<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvBICBackward[[i]] 	<- ApplyGPSurvBIC(trainingTestStructure,allParameterStructures)
	c.index.GPSurvBICBackward[i] 			<- ifelse(length(outputStructureGPSurvBICBackward[[i]]$c.index)!=0,outputStructureGPSurvBICBackward[[i]]$c.index,NA)
	rmse.GPSurvBICBackward[i] 				<- ifelse(length(outputStructureGPSurvBICBackward[[i]]$rmse)!=0,outputStructureGPSurvBICBackward[[i]]$rmse,NA)
	brier.GPSurvBICBackward[i] 				<- ifelse(length(outputStructureGPSurvBICBackward[[i]]$brier)!=0,outputStructureGPSurvBICBackward[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvBICBackward'))
save(list='outputStructureGPSurvBICBackward',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvBICBackward','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##----------- Run GPS3 SqExp Model Applied to Bootstrapped Feature Subsets ------------##
##-------------------------------------------------------------------------------------##
nBootstraps 							<- 200 
subsetDimension 						<- c(3) 		# e.g. c(3,4)
toCount 								<- 'features'
parameterStructure$covFuncForm 			<- 'SqExp'
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'

outputStructureGPSurvSqExpRSFS 			<- rep(list(list()),nReps)
c.index.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
rmse.GPSurvSqExpRSFS 					<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
brier.GPSurvSqExpRSFS 					<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
bic.GPSurvSqExpRSFS 					<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
ensembleStructure 						<- list()
subsetIndices 							<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
c.index.GPSurvSqExpRSFS.ensemble 		<- numeric()
rmse.GPSurvSqExpRSFS.ensemble 			<- numeric()
brier.GPSurvSqExpRSFS.ensemble 			<- numeric()
c.index.GPSurvSqExpRSFS.ensembleUnif 	<- numeric()
rmse.GPSurvSqExpRSFS.ensembleUnif 		<- numeric()
brier.GPSurvSqExpRSFS.ensembleUnif 		<- numeric()

for(i in 1:nReps){
	subsetIndices[i,] <- SubsetIndices(dimension,subsetDimension,nBootstraps)
	for(j in 1:nBootstraps){
		trainingTestStructure 						<- ChooseSubsetsRandomly(trainingTestStructureForNReps[[i]],subsetDimension,subsetIndices[i,j])
		dataOptionsStructure$dimension 				<- trainingTestStructure$dimension
		parameterStructure$logHypStart 				<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(2),'mean'=c(rep(0,dataOptionsStructure$dimension),0))
		outputStructureGPSurvSqExpRSFS[[i]][[j]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
		c.index.GPSurvSqExpRSFS[i,j] 				<- ifelse(length(outputStructureGPSurvSqExpRSFS[[i]][[j]]$c.index)!=0,outputStructureGPSurvSqExpRSFS[[i]][[j]]$c.index,NA)
		rmse.GPSurvSqExpRSFS[i,j] 					<- ifelse(length(outputStructureGPSurvSqExpRSFS[[i]][[j]]$rmse)!=0,outputStructureGPSurvSqExpRSFS[[i]][[j]]$rmse,NA)
		brier.GPSurvSqExpRSFS[i,j] 					<- ifelse(length(outputStructureGPSurvSqExpRSFS[[i]][[j]]$brier)!=0,outputStructureGPSurvSqExpRSFS[[i]][[j]]$brier,NA)
		bic.GPSurvSqExpRSFS[i,j] 					<- CalculateBIC(outputStructureGPSurvSqExpRSFS[[i]][[j]],toCount)
	}
	ensembleStructure[[i]] <- CalculateEnsembleResults(outputStructureGPSurvSqExpRSFS[[i]],trainingTestStructureForNReps[[i]],bic.GPSurvSqExpRSFS[i,],c.index.GPSurvSqExpRSFS[i,])
	c.index.GPSurvSqExpRSFS.ensemble[i] 			<- ensembleStructure[[i]]$ensembleMetrics$c.index
	rmse.GPSurvSqExpRSFS.ensemble[i] 				<- ensembleStructure[[i]]$ensembleMetrics$rmse
	brier.GPSurvSqExpRSFS.ensemble[i] 				<- ensembleStructure[[i]]$ensembleMetrics$brier
	c.index.GPSurvSqExpRSFS.ensembleUnif[i] 		<- ensembleStructure[[i]]$ensembleMetricsUnifWeight$c.index
	rmse.GPSurvSqExpRSFS.ensembleUnif[i] 			<- ensembleStructure[[i]]$ensembleMetricsUnifWeight$rmse
	brier.GPSurvSqExpRSFS.ensembleUnif[i] 			<- ensembleStructure[[i]]$ensembleMetricsUnifWeight$brier
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExpRSFS'))
save(list='outputStructureGPSurvSqExpRSFS',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExpRSFS','_','Workspace.RData'))
save(list='ensembleStructure',file=paste0('Runs','/',unid,'/',unid,'_','ensembleStructure','_','Workspace.RData'))
	

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
save(list='outputStructureRF',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureRF','_','Workspace.RData'))


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
save(list='outputStructureRFSurvival',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureRFSurvival','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##---------------------------------- Run Coxph Model ----------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureCoxph[[i]] 			<- ApplyCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Coxph[i] 					<- ifelse(length(outputStructureCoxph[[i]]$c.index)!=0,outputStructureCoxph[[i]]$c.index,NA)
	rmse.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$rmse)!=0,outputStructureCoxph[[i]]$rmse,NA)
	brier.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$brier)!=0,outputStructureCoxph[[i]]$brier,NA)
}
save(list='outputStructureCoxph',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureCoxph','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##------------------------------ Run Filter Coxph Model 1 -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
parameterStructure$filterMethod 		<- 'univariate'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureFilterCoxphU[[i]] 	<- ApplyFilterCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.FilterCoxphU[i] 			<- ifelse(length(outputStructureFilterCoxphU[[i]]$c.index)!=0,outputStructureFilterCoxphU[[i]]$c.index,NA)
	rmse.FilterCoxphU[i] 				<- ifelse(length(outputStructureFilterCoxphU[[i]]$rmse)!=0,outputStructureFilterCoxphU[[i]]$rmse,NA)
	brier.FilterCoxphU[i] 				<- ifelse(length(outputStructureFilterCoxphU[[i]]$brier)!=0,outputStructureFilterCoxphU[[i]]$brier,NA)
}
save(list='outputStructureFilterCoxphU',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureFilterCoxphU','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##------------------------------ Run Filter Coxph Model 2 -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
parameterStructure$filterMethod 		<- 'cox'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureFilterCoxphC[[i]] 	<- ApplyFilterCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.FilterCoxphC[i] 			<- ifelse(length(outputStructureFilterCoxphC[[i]]$c.index)!=0,outputStructureFilterCoxphC[[i]]$c.index,NA)
	rmse.FilterCoxphC[i] 				<- ifelse(length(outputStructureFilterCoxphC[[i]]$rmse)!=0,outputStructureFilterCoxphC[[i]]$rmse,NA)
	brier.FilterCoxphC[i] 				<- ifelse(length(outputStructureFilterCoxphC[[i]]$brier)!=0,outputStructureFilterCoxphC[[i]]$brier,NA)
}
save(list='outputStructureFilterCoxphC',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureFilterCoxphC','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##---------------------------------- Run Coxnet Model ---------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureCoxnet[[i]] 			<- ApplyCoxnet(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.Coxnet[i] 					<- ifelse(length(outputStructureCoxnet[[i]]$c.index)!=0,outputStructureCoxnet[[i]]$c.index,NA)
	rmse.Coxnet[i] 						<- ifelse(length(outputStructureCoxnet[[i]]$rmse)!=0,outputStructureCoxnet[[i]]$rmse,NA)
	brier.Coxnet[i] 					<- ifelse(length(outputStructureCoxnet[[i]]$brier)!=0,outputStructureCoxnet[[i]]$brier,NA)
}
save(list='outputStructureCoxnet',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureCoxnet','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##----------------------------- Run Stepwise Coxph Model ------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureStepCoxph[[i]] 		<- ApplyStepCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.StepCoxph[i] 				<- ifelse(length(outputStructureStepCoxph[[i]]$c.index)!=0,outputStructureStepCoxph[[i]]$c.index,NA)
	rmse.StepCoxph[i] 					<- ifelse(length(outputStructureStepCoxph[[i]]$rmse)!=0,outputStructureStepCoxph[[i]]$rmse,NA)
	brier.StepCoxph[i] 					<- ifelse(length(outputStructureStepCoxph[[i]]$brier)!=0,outputStructureStepCoxph[[i]]$brier,NA)
}
save(list='outputStructureStepCoxph',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureStepCoxph','_','Workspace.RData'))



##-------------------------------------------------------------------------------------##
##------------------------------ Print and Save Results -------------------------------##
##-------------------------------------------------------------------------------------##
							##--------------------------------##
							##--- Plot Concordance Indices ---##
							##--------------------------------##
modelNames 		<- c('RF','RFSurvival','Coxph','Coxnet','StepCoxph','GPNonSurvNoCens','GPSurvSqExp','GPSurvARD','GPSurvBICForward','GPSurvBICBackward','GPSurvInformedARD','GPSurvInformedARD2','GPSurvSqExpRSFS.ensemble') # Reordered version of modelsList
modelNamesPlot 	<- c('RF','RFSurvival','Coxph','Coxnet','StepCoxph','GP','GPS3SqExp','GPS3ARD','GPS3BICForward','GPS3BICBackward','GPS3IARD1','GPS3IARD2','GPS3SqExpRSFS')
toInvert 		<- c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
cols 			<- c(rep('green4',5),rep('chartreuse4',5),rep('chartreuse3',3))
c.index.mean 	<- numeric()
c.index.mat 	<- matrix(0,nrow=nReps,ncol=length(modelNames))
for(i in 1:length(modelNames)){
	c.index.mat[,i] 	<- get(paste0('c.index.',modelNames[i]))
	if(toInvert[i]) c.index.mat[,i] <- 1-c.index.mat[,i]
	c.index.mean[i] 	<- mean(c.index.mat[,i],na.rm=TRUE)
}

if(plotSaveOptions$savePlots){
	pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotCIndexAllModelsBoxplot.pdf'),width=8,height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(c.index.mat,col=add.alpha(cols,0.8),xaxt='n',xlab='',ylab='Concordance Index',ylim=c(0.4,1))
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}

							##--------------------------------##
							##----------- Plot BIC -----------##
							##--------------------------------##
bic.GPNonSurvNoCens 	<- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPNonSurvNoCens[[x]],'hyperparameters'))
bic.GPSurvARD 			<- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvARD[[x]],'hyperparameters'))
bic.GPSurvBICBackward 	<- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvBICBackward[[x]],'hyperparameters'))
bic.GPSurvBICForward 	<- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvBICForward[[x]],'hyperparameters'))
bic.GPSurvInformedARD 	<- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvInformedARD[[x]],'hyperparameters'))
bic.GPSurvInformedARDV2 <- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvInformedARDV2[[x]],'hyperparameters'))
bic.GPSurvInformedARDV3 <- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvInformedARDV3[[x]],'hyperparameters'))
bic.GPSurvInformedARDV4 <- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvInformedARDV4[[x]],'hyperparameters'))
bic.GPSurvSqExp 		<- sapply(1:nReps,function(x) CalculateBIC(outputStructureGPSurvSqExp[[x]],'hyperparameters'))

modelNames 				<- c('GPSurvSqExp','GPSurvARD','GPSurvBICForward','GPSurvBICBackward','GPSurvInformedARD','GPSurvInformedARDV2') # Reordered version of modelsList
modelNamesPlot 			<- c('GPS3SqExp','GPS3ARD','GPS3BICForward','GPS3BICBackward','GPS3IARD1','GPS3IARD2')
bic.mat 				<- cbind(bic.GPSurvSqExp,bic.GPSurvARD,bic.GPSurvBICForward,bic.GPSurvBICBackward,bic.GPSurvInformedARD,bic.GPSurvInformedARDV2)

cols 					<- c(rep('chartreuse4',4),rep('chartreuse3',4))
if(plotSaveOptions$savePlots){
	pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotBICAllModelsBoxplot.pdf'),width=8,height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(bic.mat,col=add.alpha(cols,0.8),xaxt='n',xlab='',ylab='BIC')
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}

							##--------------------------------##
							##----- Plot Hyperparameters -----##
							##--------------------------------##
nBootstraps <- length(outputStructureGPSurvSqExpRSFSFull[[1]])
lengthHypRSFS <- matrix(rep(NA,nReps*nBootstraps),nrow=nReps)
for(i in 1:nReps){
	lengthHypRSFS[i,] <- sapply(1:nBootstraps,function(x) outputStructureGPSurvSqExpRSFSFull[[i]][[x]]$logHypChosen$length)
}
summary(exp(c(lengthHypRSFS)))

if(plotSaveOptions$savePlots){
	pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotHyperparam_1IR_20D.pdf'),width=8,height=6)
	pdf.output <- dev.cur()
		PostProcessingHyperparameterBoxplots(modelNames,nReps,FALSE,TRUE)
	dev.off(pdf.output)
}