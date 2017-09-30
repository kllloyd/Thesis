#----------------------------------------------------------------------------------------#
# K Lloyd 2017_08_09
#----------------------------------------------------------------------------------------#
# Applying models to data extracted from package CuratedOvarianData.
# Data set from Tothill et al., 2008. Gene expression and clinical measurements with survival data.
# All models applied to the same data. Models are:	GP, GPS3SqExp, GPS3ARD, GPS3IARD1, GPS3IARD2, GPS3SqExpRSFS,
# 													GPSurvBICForward, GPSurvBICBackward,
# 													Cox PH, Coxnet, StepCox, FilterCoxph1, FilterCoxph2,
# 													RF and RSF
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Compare c index results of different models
# To be run with different combinations of clinical and molecular features
# 16 combinations, change using k
#---------------------------------------------------------------------------------------#

##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##
# args = commandArgs(trailingOnly=TRUE)

library(devtools)
library(curatedOvarianData)
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
source('../toSource/GetDataSetFeatures.R')
source('../toSource/GetFeatureLists.R')
source('../toSource/ImputeMissingData.R')
source('../toSource/LogPriorX.R')
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
source('../toSource/RandomGeneSelectionCuratedOvarian.R')
source('../toSource/RemoveCensored.R')
source('../toSource/RemoveFeature.R')
source('../toSource/RenameGeneFeaturesAccordingToDataSetCuratedOvarian.R')
source('../toSource/RunModelCombinationsBackward.R')
source('../toSource/RunModelCombinationsForward.R')
source('../toSource/SetParametersBIC.R')
source('../toSource/SetParametersCh5RealTothill.R')
source('../toSource/SubsetIndices.R')


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
nReps 			<- 30

dataSetChoice 	<- 'OCGS+SRGS+Clin+Rand'
# dataSetChoice 	<- args
nRand 			<-  50
switch(dataSetChoice,
	'OCGS' 					={geneSubsetFlag 		<- 'TaqMan' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 97},
	'OCGS+Clinical'			={geneSubsetFlag 		<- 'TaqMan' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 100},
	'SRGS'					={geneSubsetFlag 		<- 'SRGS1'
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 84},
	'SRGS+Clinical'			={geneSubsetFlag 		<- 'SRGS1'
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 87},
	'OCGS+Rand' 			={geneSubsetFlag 		<- 'TaqMan+Rand' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 97+nRand},
	'OCGS+Clin+Rand' 		={geneSubsetFlag 		<- 'TaqMan+Rand' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 97+3+nRand},
	'SRGS+Rand' 			={geneSubsetFlag 		<- 'SRGS1+Rand' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 84+nRand},
	'SRGS+Clin+Rand' 		={geneSubsetFlag 		<- 'SRGS1+Rand' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 84+3+nRand},
	'OCGS+SRGS' 			={geneSubsetFlag 		<- 'TaqMan+SRGS1' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 173},
	'OCGS+SRGS+Clin' 		={geneSubsetFlag 		<- 'TaqMan+SRGS1' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 173+3},
	'OCGS+SRGS+Rand' 		={geneSubsetFlag 		<- 'TaqMan+SRGS1+Rand' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 173+nRand},
	'OCGS+SRGS+Clin+Rand' 	={geneSubsetFlag 		<- 'TaqMan+SRGS1+Rand' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 173+3+nRand},
	'Random' 				={geneSubsetFlag 		<- 'Rand' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- nRand},
	'OCGSsublists' 			={geneSubsetFlag 		<- 'TaqMan' 
							  clinicalSubsetFlag 	<- 'None'
							  dimension 			<- 97},
	'Clinical' 				={geneSubsetFlag 		<- 'None' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 3},
	'Clinicalsublists' 		={geneSubsetFlag 		<- 'None' 
							  clinicalSubsetFlag 	<- 'Four'
							  dimension 			<- 3})

# k 						<- as.numeric(args[1])
# unid 					<- paste0(format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S'),'k',k)
unid 					<- format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S')
outerFolder 			<- 'Runs'
folderName 				<- paste0(outerFolder,'/',unid)
dir.create(file.path(getwd(),outerFolder),showWarnings=FALSE)
allParameterStructures 	<- SetParametersCh5RealTothill(geneSubsetFlag=geneSubsetFlag,clinicalSubsetFlag=clinicalSubsetFlag,nRand=nRand,dimension=dimension,nReps=nReps,unid=unid)
plotSaveOptions 		<- allParameterStructures$plotSaveOptions
dataOptionsStructure 	<- allParameterStructures$dataOptionsStructure
parameterStructure 		<- allParameterStructures$parameterStructure

censoringType 			<- dataOptionsStructure$censoringType
tolerance 				<- parameterStructure$tolerance


##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureGPNonSurvNoCens 		<- list()
outputStructureGPSurvSqExp 			<- list()
outputStructureGPSurvARD			<- list()
outputStructureGPSurvInformedARD1 	<- list()
outputStructureGPSurvInformedARD2 	<- list()
outputStructureGPSurvInfMed 		<- list()
outputStructureGPSurvInfUnif 		<- list()
outputStructureGPSurvBICForwards 	<- list()
outputStructureGPSurvBICBackwards 	<- list()
outputStructureRF					<- list()
outputStructureRFSurvival 			<- list()
outputStructureCoxph 				<- list()
outputStructureCoxnet 				<- list()
outputStructureStepCoxph 			<- list()
outputStructureFilterCoxphC 		<- list()
outputStructureFilterCoxphU 		<- list()
c.index.GPNonSurvNoCens 			<- rep(NA,nReps)
c.index.GPSurvSqExp 				<- rep(NA,nReps)
c.index.GPSurvARD					<- rep(NA,nReps)
c.index.GPSurvInformedARD1 			<- rep(NA,nReps)
c.index.GPSurvInformedARD2 			<- rep(NA,nReps)
c.index.GPSurvInfMed 				<- rep(NA,nReps)
c.index.GPSurvInfUnif 				<- rep(NA,nReps)
c.index.GPSurvBICForwards 			<- rep(NA,nReps)
c.index.GPSurvBICBackwards 			<- rep(NA,nReps)
c.index.RF							<- rep(NA,nReps)
c.index.RFSurvival 					<- rep(NA,nReps)
c.index.Coxph 						<- rep(NA,nReps)
c.index.Coxnet 						<- rep(NA,nReps)
c.index.StepCoxph 					<- rep(NA,nReps)
c.index.FilterCoxphC 				<- rep(NA,nReps)
c.index.FilterCoxphU 				<- rep(NA,nReps)
trainingTestStructureForNReps 		<- list()


##-------------------------------------------------------------------------------------##
##----------------------------------- Generate Data -----------------------------------##
##-------------------------------------------------------------------------------------##
trainingTestStructureForNReps 			<- GenerateData(dataOptionsStructure,outerFolder,nReps)
for(i in 1:nReps){
	trainingTestStructureForNReps[[i]] 	<- NormaliseExpressionData(trainingTestStructureForNReps[[i]],normaliseFlag=TRUE,winsoriseFlag=FALSE) # normaliseFlag TRUE -> genes normalised to (mean=0,sd=1), winsoriseFlag TRUE -> outside (0.05,0.95) quantiles clipped to here
}
save(list='trainingTestStructureForNReps',file=paste0('Runs','/',unid,'/',unid,'_','trainingTestStructureForNReps','_','Workspace.RData'))

##-------------------------------------------------------------------------------------##
##----------------------------- Run GPNonSurvNoCens Model -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- 'None'
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'non-survival'
parameterStructure$covFuncForm 			<- 'SqExp'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dataOptionsStructure$dimension),0))
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	trainingTestStructure 				<- RemoveCensored(trainingTestStructure,dataOptionsStructure)
	outputStructureGPNonSurvNoCens[[i]] <- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPNonSurvNoCens[i] 			<- ifelse(length(outputStructureGPNonSurvNoCens[[i]]$c.index)!=0,outputStructureGPNonSurvNoCens[[i]]$c.index,NA)
}
save(list='outputStructureGPNonSurvNoCens',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPNonSurvNoCens','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##------------------------------- Run GPS3 SqExp Model --------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$covFuncForm 			<- 'SqExp'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(0.9),'mean'=c(rep(0,dataOptionsStructure$dimension),0))
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
parameterStructure$tolerance 			<- tolerance
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvSqExp[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvSqExp[i] 				<- ifelse(length(outputStructureGPSurvSqExp[[i]]$c.index)!=0,outputStructureGPSurvSqExp[[i]]$c.index,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExp'))
save(list='outputStructureGPSurvSqExp',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExp','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##-------------------------------- Run GPS3 ARD Model ---------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$covFuncForm 			<- 'ARD'
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=rep(log(0.9),dataOptionsStructure$dimension),'mean'=c(rep(0,dataOptionsStructure$dimension),0))
parameterStructure$noiseCorr 			<- 'noiseCorrVec'
parameterStructure$modelType 			<- 'survival'
parameterStructure$tolerance 			<- tolerance
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvARD[[i]] 		<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvARD[i] 				<- ifelse(length(outputStructureGPSurvARD[[i]]$c.index)!=0,outputStructureGPSurvARD[[i]]$c.index,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvARD'))
save(list='outputStructureGPSurvARD',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvARD','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##---------------------------- Run GPS3 InformedARD1 Model -----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 			<- censoringType
parameterStructure$covFuncForm 				<- 'InformedARD'
featureListStructure 						<- GetFeatureLists(dataSetChoice,'Inbuilt',dataOptionsStructure$randomGeneList,dataOptionsStructure$chosenClinicalFeatures,parameterStructure$covFuncForm,dataOptionsStructure$dimension)
parameterStructure$extraParam 				<- featureListStructure$extraParam
parameterStructure$logHypStart 				<- featureListStructure$logHypStart
parameterStructure$noiseCorr 				<- 'noiseCorrVec'
parameterStructure$modelType 				<- 'survival'
parameterStructure$tolerance 				<- tolerance
for(i in 1:nReps){
	trainingTestStructure 					<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInformedARD1[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInformedARD1[i] 			<- ifelse(length(outputStructureGPSurvInformedARD1[[i]]$c.index)!=0,outputStructureGPSurvInformedARD1[[i]]$c.index,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvInformedARD1'))
save(list='outputStructureGPSurvInformedARD1',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvInformedARD1','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##--------------------------- Run GPS3 InformedARDV2 Model ----------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 			<- censoringType
parameterStructure$covFuncForm 				<- 'InformedARD'
featureListStructure 						<- GetFeatureLists(dataSetChoice,'Inbuilt',dataOptionsStructure$randomGeneList,dataOptionsStructure$chosenClinicalFeatures,parameterStructure$covFuncForm,dataOptionsStructure$dimension)
parameterStructure$extraParam 				<- featureListStructure$extraParam
parameterStructure$logHypStart 				<- featureListStructure$logHypStart
parameterStructure$noiseCorr 				<- 'noiseCorrVec'
parameterStructure$modelType 				<- 'survival'
parameterStructure$tolerance 				<- tolerance
parameterStructure 							<- ModifyFeatureListsInformedARD2(dataOptionsStructure,parameterStructure,trainingTestStructureForNReps)
for(i in 1:nReps){
	trainingTestStructure 					<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvInformedARD2[[i]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvInformedARD2[i] 			<- ifelse(length(outputStructureGPSurvInformedARD2[[i]]$c.index)!=0,outputStructureGPSurvInformedARD2[[i]]$c.index,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvInformedARD2'))
save(list='outputStructureGPSurvInformedARD2',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvInformedARD2','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##----------------- Run GPS3 SqExp Model with RSFS Feature Selection ------------------##
##-------------------------------------------------------------------------------------##
nBootstraps 						<- 600
subsetDimension 					<- c(5)
toCount 							<- 'features'
dataOptionsStructure$censoringType 	<- censoringType
parameterStructure$covFuncForm 		<- 'SqExp'
parameterStructure$noiseCorr 		<- 'noiseCorrVec'
parameterStructure$modelType 		<- 'survival'
	
outputStructureGPSurvSqExpRSFS 		<- rep(list(list()),nReps)
c.index.GPSurvSqExpRSFS 			<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
bic.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
subsetIndices 						<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
ensembleStructure 					<- list()
for(i in 1:nReps){
	subsetIndices <- SubsetIndices(dimension,subsetDimension,nBootstraps)
	for(j in 1:nBootstraps){
		trainingTestStructure 						<- ChooseSubsetsRandomly(trainingTestStructureForNReps[[i]],subsetDimension,subsetIndices[i,j])
		dataOptionsStructure$dimension 				<- trainingTestStructure$dimension
		parameterStructure$logHypStart 				<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(2),'mean'=c(rep(0,dataOptionsStructure$dimension),0))
		outputStructureGPSurvSqExpRSFS[[i]][[j]] 	<- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
		c.index.GPSurvSqExpRSFS[i,j] 				<- ifelse(length(outputStructureGPSurvSqExpRSFS[[i]][[j]]$c.index)!=0,outputStructureGPSurvSqExpRSFS[[i]][[j]]$c.index,NA)
		bic.GPSurvSqExpRSFS[i,j] 					<- CalculateBIC(outputStructureGPSurvSqExpRSFS[[i]][[j]],toCount)
	}
	ensembleStructure[[i]] 				<- CalculateEnsembleResults(outputStructureGPSurvSqExpRSFS[[i]],trainingTestStructureForNReps[[i]],bic.GPSurvSqExpRSFS[i,],c.index.GPSurvSqExpRSFS[i,])
	c.index.GPSurvSqExpRSFS.ensemble[i] <- ensembleStructure[[i]]$ensembleMetrics$c.index
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExpRSFS'))
save(list='outputStructureGPSurvSqExpRSFS',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExpRSFS','_','Workspace.RData'))
save(list='ensembleStructure',file=paste0('Runs','/',unid,'/',unid,'_','ensembleStructureGPSurvSqExpRSFS','_','Workspace.RData'))


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
}
file.rename(from=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/GPSurvCorrV'),to=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/GPSurvBICForward'))
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
}
file.rename(from=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/GPSurvCorrV'),to=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/GPSurvBICBackward'))
save(list='outputStructureGPSurvBICBackward',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvBICBackward','_','Workspace.RData'))


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
}
save(list='outputStructureCoxph',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureCoxph','_','Workspace.RData'))


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
}
save(list='outputStructureStepCoxph',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureStepCoxph','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##----------------------------- Run Filter Coxph 1 Model ------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
parameterStructure$filterMethod 		<- 'univariate'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureFilterCoxphU[[i]] 	<- ApplyFilterCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.FilterCoxphU[i] 			<- ifelse(length(outputStructureFilterCoxphU[[i]]$c.index)!=0,outputStructureFilterCoxphU[[i]]$c.index,NA)
}
save(list='outputStructureFilterCoxphU',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureFilterCoxphU','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##----------------------------- Run Filter Coxph 2 Model ------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 		<- censoringType
parameterStructure$noiseCorr 			<- FALSE
parameterStructure$modelType 			<- 'survival'
parameterStructure$filterMethod 		<- 'cox'
for(i in 1:nReps){
	trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
	outputStructureFilterCoxphC[[i]] 	<- ApplyFilterCoxph(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.FilterCoxphC[i] 			<- ifelse(length(outputStructureFilterCoxphC[[i]]$c.index)!=0,outputStructureFilterCoxphC[[i]]$c.index,NA)
}
save(list='outputStructureFilterCoxphC',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureFilterCoxphC','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##------------------------------ Print and Save Results -------------------------------##
##-------------------------------------------------------------------------------------##
							##--------------------------------##
							##--- Plot Concordance Indices ---##
							##--------------------------------##
modelNames 		<- c('RF','RFSurvival','Coxph','FilterCoxphU','FilterCoxphC','Coxnet','StepCoxph','GPNonSurvNoCens','GPSurvSqExp','GPSurvARD','GPSurvBICForward','GPSurvBICBackward','GPSurvInformedARD1','GPSurvInformedARD2','GPSurvSqExpRSFS.ensemble') # Reordered version of modelsList
modelNamesPlot 	<- c('RF','RFSurvival','Coxph','FilterCoxph1','FilterCoxph2','Coxnet','StepCoxph','GP','GPS3SqExp','GPS3ARD','GPS3BICForward','GPS3BICBackward','GPS3IARD1','GPS3IARD2','GPS3SqExpRSFS')
toInvert 		<- c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
c.index.mean 	<- numeric()
c.index.mat 	<- matrix(0,nrow=nReps,ncol=length(modelNames))
for(i in 1:length(modelNames)){
	c.index.mat[,i] 	<- get(paste0('c.index.',modelNames[i]))
	if(toInvert[i]) c.index.mat[,i] <- 1-c.index.mat[,i]
	c.index.mean[i] 	<- mean(c.index.mat[,i],na.rm=TRUE)
}

cols 			<- c(rep('green4',7),rep('chartreuse4',5),rep('chartreuse3',3))

if(allParameterStructures$plotSaveOptions$savePlots){
	pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotCIndexAllModelsBoxplot.pdf'),width=8,height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(8,1))
		boxplot(c.index.mat,col=add.alpha(cols,0.8),xaxt='n',xlab='',ylab='Concordance Index',ylim=c(0.4,1))
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}