#----------------------------------------------------------------------------------------#
# K Lloyd 2016_09_16
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# All models applied to the same data. Models are:	GP, GPS1, GPS2 and GPS3
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Making 1D plots to show fitted curve and training targets before and after learning
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
source('../toSource/PlotTrainingTargetsPrePostLearning.R')
source('../toSource/PlotKaplanMeier.R')
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

nReps 					<- 100

modelsList 				<- c('GPNonSurvNoCens','GPSurvNoCorr','GPSurvCorrL','GPSurvCorrV')


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
dimension 				<- 1
proportionTest 			<- NA
nTraining 				<- 100
nTest 					<- 50 

logHypGenerate 			<- list('noise'=log(0.1),'func'=log(0.7),'length'=log(1.5),'mean'=c(rep(0,dimension),0))

covFuncFormGen 			<- 'SqExp'
maternParamGen 			<- 3
meanFuncFormGen 		<- 'Linear'
gridMinimum 			<- 0
gridMaximum 			<- 8
censoringType 			<- 'NormalLoopSample'
censoringSD 			<- 50
censoringMean 			<- 0
censoredProportion 		<- 0.40
nCensored 				<- ceiling((nTraining+nTest)*censoredProportion)
if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE
extraDimensions 		<- 0 

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
tolerance 				<- 0.00001*nTraining*censoredProportion
toleranceLaplace 		<- 0.00001*nTraining*censoredProportion
hypChangeTol 			<- 10^-10
maxCount 				<- 500
maxit 					<- 2000
maxitPreLearn 			<- 1000
maxitSurvival 			<- 10
maxitLaplaceHyp			<- 10
maxitLaplaceFHat 		<- 100
optimType 				<- 'Nelder-Mead' 					# 'CG' = conjugate gradient optimisation, 'Nelder-Mead' = numerical gradients
noiseCorr 				<- FALSE
imposePriors 			<- TRUE
logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(1.7),'mean'=c(rep(0,dimension),0))
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

c.index.GPNonSurvNoCens 		<- rep(NA,nReps)
c.index.GPSurvNoCorr 			<- rep(NA,nReps)
c.index.GPSurvCorrV				<- rep(NA,nReps)
c.index.GPSurvCorrL 			<- rep(NA,nReps)

rmse.GPNonSurvNoCens 			<- rep(NA,nReps)
rmse.GPSurvNoCorr 				<- rep(NA,nReps)
rmse.GPSurvCorrV				<- rep(NA,nReps)
rmse.GPSurvCorrL 				<- rep(NA,nReps)

brier.GPNonSurvNoCens 			<- rep(NA,nReps)
brier.GPSurvNoCorr 				<- rep(NA,nReps)
brier.GPSurvCorrV				<- rep(NA,nReps)
brier.GPSurvCorrL 				<- rep(NA,nReps)


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
##------------------------------ Print and Save Results -------------------------------##
##-------------------------------------------------------------------------------------##
							##-------------------------------##
							##---- Print C Index & RMSE -----##
							##-------------------------------##
if(printResults){
	cat('---------------------------------------',fill=TRUE)
	cat('Concordance Indices:',fill=TRUE)
	cat('GPNonSurvNoCens mean c index =',paste0(round(c.index.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPSurvNoCorr mean c index =',paste0(round(c.index.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPSurvCorrV mean c index =',paste0(round(c.index.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('GPSurvCorrL mean c index =',paste0(round(c.index.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
	cat('RMSE:',fill=TRUE)
	cat('GPNonSurvNoCens mean rmse =',paste0(round(rmse.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPSurvNoCorr mean rmse =',paste0(round(rmse.GPSurvNoCorr,4),collapse=', '),fill=TRUE)
	cat('GPSurvCorrV mean rmse =',paste0(round(rmse.GPSurvCorrV,4),collapse=', '),fill=TRUE)
	cat('GPSurvCorrL mean rmse =',paste0(round(rmse.GPSurvCorrL,4),collapse=', '),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
}

							##-------------------------------##
							##--- Plot Kaplan-Meier Plots ---##
							##-------------------------------##
if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GPNonSurv','PlotKaplanMeier.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPNonSurvNoCens[[i]]$funcMeanPred,outputStructureGPNonSurvNoCens[[i]]$trainingTestStructure$testTargets,outputStructureGPNonSurvNoCens[[i]]$trainingTestStructure$testEvents,'GPNonSurv')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','PlotKaplanMeier.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvNoCorr[[i]]$funcMeanPred,outputStructureGPSurvNoCorr[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvNoCorr[[i]]$trainingTestStructure$testEvents,'GPSurvNoCorr')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','PlotKaplanMeier.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvCorrL[[i]]$funcMeanPred,outputStructureGPSurvCorrL[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvCorrL[[i]]$trainingTestStructure$testEvents,'GPSurvCorrL')
		}
	dev.off(pdf.output)

	pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','PlotKaplanMeier.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotKaplanMeier(outputStructureGPSurvCorrV[[i]]$funcMeanPred,outputStructureGPSurvCorrV[[i]]$trainingTestStructure$testTargetsPreCensoring,outputStructureGPSurvCorrV[[i]]$trainingTestStructure$testEvents,'GPSurvCorrV')
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
}

							##--------------------------------##
							##--- Plot 1D training targets ---##
							##--------------------------------##
if(savePlots){
	for(i in 1:nReps){
		pdf(paste0('Runs/',unid,'/','GPSurvNoCorr','/',unid,'GPSurvNoCorr','Run',i,'PlotTrainingTargets1D.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
		PlotTrainingTargetsPrePostLearning(outputStructureGPSurvNoCorr[[i]])
		dev.off(pdf.output)
	}

	for(i in 1:nReps){
		pdf(paste0('Runs/',unid,'/','GPSurvCorrL','/',unid,'GPSurvCorrL','Run',i,'PlotTrainingTargets1D.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
		PlotTrainingTargetsPrePostLearning(outputStructureGPSurvCorrL[[i]])
		dev.off(pdf.output)
	}

	for(i in 1:nReps){
		pdf(paste0('Runs/',unid,'/','GPSurvCorrV','/',unid,'GPSurvCorrV','Run',i,'PlotTrainingTargets1D.pdf'),width=10,height=8,onefile=TRUE)
		pdf.output <- dev.cur()
		PlotTrainingTargetsPrePostLearning(outputStructureGPSurvCorrV[[i]])
		dev.off(pdf.output)
	}
}

##-------------------------------------------------------------------------------------##
##------------------------------------ Save Output ------------------------------------##
##-------------------------------------------------------------------------------------##
outputStructureAll <- list('outputStructureGPNonSurvNoCens'=outputStructureGPNonSurvNoCens,
							'outputStructureGPSurvNoCorr'=outputStructureGPSurvNoCorr,
							'outputStructureGPSurvCorrV'=outputStructureGPSurvCorrV,
							'outputStructureGPSurvCorrL'=outputStructureGPSurvCorrL)

save(list='outputStructureAll',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureAll','_','Workspace.RData'))