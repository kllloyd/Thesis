#----------------------------------------------------------------------------------------#
# K Lloyd 2017_08_09
#----------------------------------------------------------------------------------------#
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# Models investigating feature selection
# All models applied to the same data. Models are:	GP, GPS3SqExp and GPS3SqExpRSFS
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Compare c index results of different models
# To be run with different numbers of data dimensions, informative and non-informative
# 32 combinations, change using k
#---------------------------------------------------------------------------------------#

##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##
# args = commandArgs(trailingOnly=TRUE)

library(akima)
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
source('../toSource/AdjustTrainingSurvivalMeanVariance.R')
source('../toSource/ApplyAFT.R')
source('../toSource/ApplyCoxph.R')
source('../toSource/ApplyGBM.R')
source('../toSource/ApplyCoxnet.R')
source('../toSource/ApplyGP.R')
source('../toSource/ApplyGPLaplace.R')
source('../toSource/ApplyGPDiffInf.R')
source('../toSource/ApplyRF.R')
source('../toSource/ApplyRFSurvival.R')
source('../toSource/ApplyStepCoxph.R')
source('../toSource/CalculateMetrics.R')
source('../toSource/CensorData.R')
source('../toSource/CovFunc.R')
source('../toSource/DataExtraction.R')
source('../toSource/GenerateData.R')
source('../toSource/ImputeMissingData.R')
source('../toSource/LogPriorX.R')
source('../toSource/MakeSyntheticData.R')
source('../toSource/MeanFunc.R')
source('../toSource/ModifyFeatureListsInformedARDV2.R')
source('../toSource/NormaliseExpressionData.R')
source('../toSource/PlotHyperparamAnyKernel.R')
source('../toSource/PlotKaplanMeier.R')
source('../toSource/PreLearnHyperparam.R')
source('../toSource/PrintOptions.R')
source('../toSource/PrintOptions2.R')
source('../toSource/RemoveCensored.R')
source('../toSource/CalculateBIC.R')
source('../toSource/SubsetIndices.R')
source('../toSource/ChooseSubsetsRandomly.R')
source('../toSource/CalculateEnsembleResults.R')
source('../toSource/scatter_fill.R')

# TO BE RUN FOR k <- 1 to k <- length(dimensionToRun)
# k 						<- as.numeric(args[1])
k 						<- 1

##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
unid 					<- paste0(format(Sys.time(),format='y%Ym%md%dh%Hm%Ms%S'),'k',k)
outerFolder 			<- 'Runs'
folderName 				<- paste0(outerFolder,'/',unid)
dir.create(file.path(getwd(),outerFolder),showWarnings=FALSE)

nReps 					<- 50
nBootstrapsToRun		<- c(10,10,100,100,100,150,150,150,150,200,200,200,200,300,300,300,300,500,500,500,500,750,750,750,750,750,1000,1000,1000,1000,1000,1000)
nBootstraps 			<- nBootstrapsToRun[k]

modelsList 				<- c('GPNonSurvNoCens','GPSurvSqExp','GPSurvSqExpRSFS')


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
dimensionToRun 			<- c(5,5,10,10,10,15,15,15,15,20,20,20,20,30,30,30,30,50,50,50,50,75,75,75,75,75,100,100,100,100,100,100)
dimension 				<- dimensionToRun[k]
extraDimensionsToRun 	<- c(2,4,2,5,8,3,7,10,13,3,10,15,18,3,10,19,27,3,15,30,45,3,25,40,55,70,3,20,40,60,75,90)
extraDimensions 		<- extraDimensionsToRun[k]
proportionTest 			<- NA
nTraining 				<- 500
nTest 					<- 100 

covFuncFormGen 			<- 'ARD'
extraParamGen 			<- 3
meanFuncFormGen 		<- 'Linear'
if(covFuncFormGen=='ARD') useARD <- TRUE else useARD <- FALSE

if(useARD){
	logHypGenerate 		<- list('noise'=log(0.01),'func'=log(0.5),'length'=log(rgamma(dimension-extraDimensions,shape=2,scale=3)),'mean'=c(rep(0,dimension),0))
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
covFuncForm 			<- 'ARD'
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
imposePriors 			<- FALSE
if(!useARD){
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
c.index.GPNonSurvNoCens 			<- rep(NA,nReps)
rmse.GPNonSurvNoCens 				<- rep(NA,nReps)
brier.GPNonSurvNoCens 				<- rep(NA,nReps)
outputStructureGPSurvSqExp 			<- list()
c.index.GPSurvSqExp 				<- rep(NA,nReps)
rmse.GPSurvSqExp 					<- rep(NA,nReps)
brier.GPSurvSqExp 					<- rep(NA,nReps)
outputStructureGPSurvSqExpRSFS 		<- rep(list(list()),nReps)
c.index.GPSurvSqExpRSFS 			<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
rmse.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
brier.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
trainingTestStructureForNReps 		<- list()


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

if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPNonSurv','/',unid,'GP','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPNonSurvNoCens[[i]])
		}
	dev.off(pdf.output)
}


##-------------------------------------------------------------------------------------##
##------------------------------- Run GPS3 SqExp Model  -------------------------------##
##-------------------------------------------------------------------------------------##
dataOptionsStructure$censoringType 	<- censoringType
parameterStructure$covFuncForm 		<- 'SqExp'
parameterStructure$logHypStart 		<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(2),'mean'=c(rep(0,dimension),0))
parameterStructure$noiseCorr 		<- 'noiseCorrVec'
parameterStructure$modelType 		<- 'survival'

for(i in 1:nReps){
	trainingTestStructure 			<- trainingTestStructureForNReps[[i]]
	outputStructureGPSurvSqExp[[i]] <- ApplyGP(trainingTestStructure,dataOptionsStructure,parameterStructure,plotSaveOptions)
	c.index.GPSurvSqExp[i] 			<- ifelse(length(outputStructureGPSurvSqExp[[i]]$c.index)!=0,outputStructureGPSurvSqExp[[i]]$c.index,NA)
	rmse.GPSurvSqExp[i] 			<- ifelse(length(outputStructureGPSurvSqExp[[i]]$rmse)!=0,outputStructureGPSurvSqExp[[i]]$rmse,NA)
	brier.GPSurvSqExp[i] 			<- ifelse(length(outputStructureGPSurvSqExp[[i]]$brier)!=0,outputStructureGPSurvSqExp[[i]]$brier,NA)
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExp'))
save(list='outputStructureGPSurvSqExp',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExp','_','Workspace.RData'))

if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPSurvSqExp','/',unid,'GPS3SqExp','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in c(1:nReps)){
			PlotHyperparamAnyKernel(outputStructureGPSurvSqExp[[i]])
		}
	dev.off(pdf.output)
}


##-------------------------------------------------------------------------------------##
##----------- Run GPS3 SqExp Model Applied to Bootstrapped Feature Subsets ------------##
##-------------------------------------------------------------------------------------##
subsetDimension 					<- min(round(dimension/3),10) 		# e.g. c(3,4)
toCount 							<- 'features'
bic.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
ensembleStructure 					<- list()
subsetIndices 						<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
c.index.GPSurvSqExpRSFS.ensemble 	<- numeric()
rmse.GPSurvSqExpRSFS.ensemble 		<- numeric()
brier.GPSurvSqExpRSFS.ensemble 		<- numeric()

parameterStructure$covFuncForm 		<- 'SqExp'
parameterStructure$noiseCorr 		<- 'noiseCorrVec'
parameterStructure$modelType 		<- 'survival'

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
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExpRSFS'))
save(list='outputStructureGPSurvSqExpRSFS',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExpRSFS','_','Workspace.RData'))
save(list='ensembleStructure',file=paste0('Runs','/',unid,'/',unid,'_','ensembleStructureGPSurvSqExpRSFS','_','Workspace.RData'))

if(savePlots){
	pdf(paste0('Runs/',unid,'/','GPSurvSqExpRSFS','/',unid,'GPS3SqExpRSFS','PlotHyperparam.pdf'),width=10,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
		for(i in 1:nReps){
			for(j in c(1:nBootstraps)){
				PlotHyperparamAnyKernel(outputStructureGPSurvSqExpRSFS[[i]][[j]])
			}
		}
	dev.off(pdf.output)
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
	cat('GPNonSurvNoCens c index =',paste0(round(c.index.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS3SqExp c index =',paste0(round(c.index.GPSurvSqExp,4),collapse=', '),fill=TRUE)
	cat('GPS3SqExpRSFS c index =',sapply(1:nReps,function(x) paste0(round(ensembleStructure[[x]]$ensembleMetrics$c.index,4),collapse=', ')),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
	cat('RMSE:',fill=TRUE)
	cat('GPNonSurvNoCens rmse =',paste0(round(rmse.GPNonSurvNoCens,4),collapse=', '),fill=TRUE)
	cat('GPS3SqExp rmse =',paste0(round(rmse.GPSurvSqExp,4),collapse=', '),fill=TRUE)
	cat('GPS3SqExpRSFS rmse =',sapply(1:nReps,function(x) paste0(round(ensembleStructure[[x]]$ensembleMetrics$rmse,4),collapse=', ')),fill=TRUE)
	cat('---------------------------------------',fill=TRUE)
}

							##--------------------------------##
							##--- Plot Concordance Indices ---##
							##--------------------------------##
modelNames 		<- c('GPNonSurvNoCens','GPSurvSqExp','GPSurvSqExpRSFS.ensemble') # Reordered version of modelsList
modelNamesPlot 	<- c('GP','GPS3SqExp','GPS3SqExpRSFS')
toInvert 		<- c(TRUE,TRUE,TRUE)
c.index.mean 	<- numeric()
c.index.mat 	<- matrix(0,nrow=nReps,ncol=length(modelNames))
c.index.mat[,1] <- 1-c.index.GPNonSurvNoCens
c.index.mat[,2] <- 1-c.index.GPSurvSqExp
c.index.mat[,3] <- sapply(1:nReps,function(x) 1-ensembleStructure[[x]]$ensembleMetrics$c.index)

if(savePlots){
	pdf(file=paste0(getwd(),"/",outerFolder,"/",unid,'/',unid,'PlotCIndexAllModelsBoxplot.pdf'),width=8, height=6)
	pdf.output <- dev.cur()
		layout(rbind(1,2), heights=c(10,1))
		boxplot(c.index.mat,col=add.alpha('chartreuse3',0.8),xaxt='n',xlab='',ylab='Concordance Index',ylim=c(0.4,1))
		axis(1,at=1:length(modelNames),labels=modelNamesPlot,las=2)
		layout(1)
	dev.off(pdf.output)
}

##-------------------------------------------------------------------------------------##
##---------------------------------- PostProcessing -----------------------------------##
##---------- To be applied after running all ks up to length(dimensionToRun) ----------##
##-------------------------------------------------------------------------------------##
# unids <- c('','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','')
# forBreakdown 			<- c(1,2,3,5,6,9,10,13,14,17,18,21,22,26,27,32)

# c.index.GPNonSurvNoCens.all 	<- matrix(0,nrow=length(dimensionsToRun),ncol=nReps)
# c.index.GPSurvSqExp.all 		<- matrix(0,nrow=length(dimensionsToRun),ncol=nReps)
# c.index.GPSurvRSFS.all 			<- matrix(0,nrow=length(dimensionsToRun),ncol=nReps)
# c.index.GPNonSurvNoCens.mean 	<- numeric()
# c.index.GPSurvSqExp.mean 		<- numeric()
# c.index.GPSurvRSFS.mean 		<- numeric()
# c.index.GPNonSurvNoCens.median 	<- numeric()
# c.index.GPSurvSqExp.median 		<- numeric()
# c.index.GPSurvRSFS.median 		<- numeric()
# c.index.GPSurvRSFS 				<- list()
# bic.GPSurvRSFS 					<- list()
# ensembleStructure 				<- rep(list(list()),length(dimensionsToRun))
# trainingTestStructureForNReps 	<- list()

# modelNamesPlot 					<- c('GP','GPS3SqExp','GPS3RSFS')
# toInvert 						<- c(TRUE,TRUE,TRUE)
# toCount 						<- 'features'

# for(i in 1:length(dimensionsToRun)){
# 	k <- ks[i]
# 	load(file=paste0(getwd(),'/Runs/',unids[i],'/',unids[i],'_outputStructureGPNonSurvNoCens_Workspace.RData'))
# 	load(file=paste0(getwd(),'/Runs/',unids[i],'/',unids[i],'_outputStructureGPSurvSqExp_Workspace.RData'))
# 	load(file=paste0(getwd(),'/Runs/',unids[i],'/',unids[i],'_outputStructureGPSurvSqExpRSFS_Workspace.RData'))
# 	assign(paste0(unids[i],'outputStructureGPNonSurvNoCens'),outputStructureGPNonSurvNoCens)
# 	assign(paste0(unids[i],'outputStructureGPSurvSqExp'),outputStructureGPSurvSqExp)
# 	assign(paste0(unids[i],'outputStructureGPSurvSqExpRSFS'),outputStructureGPSurvSqExpRSFS)
# 	rm('outputStructureGPNonSurvNoCens','outputStructureGPSurvSqExp','outputStructureGPSurvSqExpRSFS')

# 	c.index.GPSurvRSFS[[i]] 				<-  matrix(rep(NA,nBootstrapsToRun[i]*nReps),ncol=nBootstrapsToRun[i])
# 	bic.GPSurvRSFS[[i]] 					<-  matrix(rep(NA,nBootstrapsToRun[i]*nReps),ncol=nBootstrapsToRun[i])
# 	for(j in 1:nReps){
# 		c.index.GPNonSurvNoCens.all[i,j] 	<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPNonSurvNoCens[[',j,']]$c.index')))
# 		c.index.GPSurvSqExp.all[i,j] 		<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPSurvSqExp[[',j,']]$c.index')))
# 		for(k in 1:nBootstrapsToRun[i]){
# 			c.index.GPSurvRSFS[[i]][j,k] 	<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPSurvSqExpRSFS[[',j,']][[',k,']]$c.index')))
# 			bic.GPSurvRSFS[[i]][j,k] 		<- CalculateBIC(eval(parse(text=paste0(unids[i],'outputStructureGPSurvSqExpRSFS[[',j,']][[',k,']]'))),toCount)
# 		}
# 		trainingTestStructureForNReps[[j]]  <- eval(parse(text=paste0(unids[i],'outputStructureGPSurvSqExp[[',j,']]$trainingTestStructure')))
# 		ensembleStructure[[i]][[j]] 		<- CalculateEnsembleResults(eval(parse(text=paste0(unids[i],'outputStructureGPSurvSqExpRSFS[[',j,']]'))),trainingTestStructureForNReps[[j]],bic.GPSurvRSFS[[i]][j,],c.index.GPSurvRSFS[[i]][j,])
# 		c.index.GPSurvRSFS.all[i,j] 		<- 1-ensembleStructure[[i]][[j]]$ensembleMetrics$c.index
# 	}
# 	c.index.GPNonSurvNoCens.mean[i] 		<- mean(c.index.GPNonSurvNoCens.all[i,],na.rm=TRUE)
# 	c.index.GPSurvSqExp.mean[i] 			<- mean(c.index.GPSurvSqExp.all[i,],na.rm=TRUE)
# 	c.index.GPSurvRSFS.mean[i] 				<- mean(c.index.GPSurvRSFS.all[i,],na.rm=TRUE)
# 	c.index.GPNonSurvNoCens.median[i] 		<- median(c.index.GPNonSurvNoCens.all[i,],na.rm=TRUE)
# 	c.index.GPSurvSqExp.median[i] 			<- median(c.index.GPSurvSqExp.all[i,],na.rm=TRUE)
# 	c.index.GPSurvRSFS.median[i] 			<- median(c.index.GPSurvRSFS.all[i,],na.rm=TRUE)
# 	rm(list=paste0(unids[i],'outputStructureGPNonSurvNoCens'))
# 	rm(list=paste0(unids[i],'outputStructureGPSurvSqExp'))
# 	rm(list=paste0(unids[i],'outputStructureGPSurvSqExpRSFS'))
# 	save(c.index.GPSurvRSFS.all,file='c.index.GPSurvRSFS.all.RData')
# 	save(c.index.GPNonSurvNoCens.all,file='c.index.GPNonSurvNoCens.all.RData')
# 	save(c.index.GPSurvSqExp.all,file='c.index.GPSurvSqExp.all.RData')
# }

# interp.GPNonSurvNoCens 	<- interp(x=extraDimensionsToRun,y=dimensionsToRun,z=c.index.GPNonSurvNoCens.mean)
# interp.GPSurvSqExp 		<- interp(x=extraDimensionsToRun,y=dimensionsToRun,z=c.index.GPSurvSqExp.mean)
# interp.GPSurvRSFS 		<- interp(x=extraDimensionsToRun,y=dimensionsToRun,z=c.index.GPSurvRSFS.mean)
# interpPlotRange 	<- range(c(interp.GPNonSurvNoCens$z,interp.GPSurvSqExp$z,interp.GPSurvRSFS$z),finite=TRUE)

# if(plotSaveOptions$savePlots){
# 	pdf(paste0('Runs/','PlotAllModelsInterpolatedContour.pdf'),width=12,height=7,onefile=TRUE)
# 	pdf.output <- dev.cur()
# 		filled.contour(interp.GPNonSurvNoCens$x,interp.GPNonSurvNoCens$y,interp.GPNonSurvNoCens$z,zlim=interpPlotRange,col=topo.colors(length(pretty(interpPlotRange,20))-1),
# 			main='',ylab='Total no. of dimensions',xlab='No. of non-informative dimensions')
# 		filled.contour(interp.GPSurvSqExp$x,interp.GPSurvSqExp$y,interp.GPSurvSqExp$z,zlim=interpPlotRange,col=topo.colors(length(pretty(interpPlotRange,20))-1),
# 			main='',ylab='Total no. of dimensions',xlab='No. of non-informative dimensions')
# 		filled.contour(interp.GPSurvRSFS$x,interp.GPSurvRSFS$y,interp.GPSurvRSFS$z,zlim=interpPlotRange,col=topo.colors(length(pretty(interpPlotRange,20))-1),
# 			main='',ylab='Total no. of dimensions',xlab='No. of non-informative dimensions')
# 	dev.off(pdf.output)

# 	pdf(paste0('Runs/','PlotAllModelsScatter.pdf'),width=10,height=8,onefile=TRUE)
# 	pdf.output <- dev.cur()
# 		scatter_fill(dimensionsToRun,extraDimensionsToRun,c.index.GPNonSurvNoCens.all,nlevels=15,xlim=c(min(dimensionsToRun),max(dimensionsToRun)),ylim=c(min(extraDimensionsToRun),max(extraDimensionsToRun)),zlim=interpPlotRange,
# 			xlab='Total no. of dimensions',ylab='No. of uninformative dimensions',main=paste0(modelNamesPlot[1],', c index'),pch=20,cex=1)
# 		scatter_fill(dimensionsToRun,extraDimensionsToRun,c.index.GPSurvSqExp.all,nlevels=15,xlim=c(min(dimensionsToRun),max(dimensionsToRun)),ylim=c(min(extraDimensionsToRun),max(extraDimensionsToRun)),zlim=interpPlotRange,
# 			xlab='Total no. of dimensions',ylab='No. of uninformative dimensions',main=paste0(modelNamesPlot[2],', c index'),pch=20,cex=1)
# 		scatter_fill(dimensionsToRun,extraDimensionsToRun,c.index.GPSurvRSFS.all,nlevels=15,xlim=c(min(dimensionsToRun),max(dimensionsToRun)),ylim=c(min(extraDimensionsToRun),max(extraDimensionsToRun)),zlim=interpPlotRange,
# 			xlab='Total no. of dimensions',ylab='No. of uninformative dimensions',main=paste0(modelNamesPlot[3],', c index'),pch=20,cex=1)
# 	dev.off(pdf.output)
# }

# breakDown 	<- matrix(forBreakdown,nrow=length(unique(dimensionsToRun)),byrow=TRUE) 
# plotMar 	<- c(2,4,1,2)+0.1

# of(plotSaveOptions$savePlots){
# 	pdf(paste0('Runs/','PlotAllModelsBoxplot.pdf'),width=12,height=12,onefile=TRUE)
# 	pdf.output <- dev.cur()
# 		layout(matrix(c(1,1,2:(2*length(unique(dimensionsToRun))+1),2*length(unique(dimensionsToRun))+2,2*length(unique(dimensionsToRun))+2),nrow=length(unique(dimensionsToRun))+2,byrow=TRUE),widths=c(1,9),heights=c(1,rep(3,length(unique(dimensionsToRun))),1))
# 		par(mar=c(0,0,0,0))
# 		plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 		text(x=0.5,y=0.5,labels=modelNamesPlot[1],font=2)
# 		par(mar=c(5,4,4,2)+0.1)
# 		for(i in 1:length(unique(dimensionsToRun))){
# 			par(mar=c(0,0,0,0))
# 			plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 			text(x=0.5,y=0.5,labels=paste('d =',unique(dimensionsToRun)[i]))
# 			par(mar=plotMar)
# 			boxplot(t(c.index.GPNonSurvNoCens.all[breakDown[i,1]:breakDown[i,2],,drop=FALSE]),at=extraDimensionsToRun[breakDown[i,1]:breakDown[i,2]],xlim=c(min(extraDimensionsToRun),max(extraDimensionsToRun)),ylab='Concordance index',ylim=c(0.4,1.0),xaxt='n')
# 			axis(1,at=unique(extraDimensionsToRun),labels=unique(extraDimensionsToRun))
# 			abline(h=0.5,col='lightgrey',lty=4)
# 		}
# 		par(mar=c(0,0,0,0))
# 		plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 		text(x=0.5,y=0.5,labels='No. uninformative dimensions')
# 		par(mar=c(5,4,4,2)+0.1)
	
# 		layout(matrix(c(1,1,2:(2*length(unique(dimensionsToRun))+1),2*length(unique(dimensionsToRun))+2,2*length(unique(dimensionsToRun))+2),nrow=length(unique(dimensionsToRun))+2,byrow=TRUE),widths=c(1,9),heights=c(1,rep(3,length(unique(dimensionsToRun))),1))
# 		par(mar=c(0,0,0,0))
# 		plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 		text(x=0.5,y=0.5,labels=modelNamesPlot[2],font=2)
# 		par(mar=c(5,4,4,2)+0.1)
# 		for(i in 1:length(unique(dimensionsToRun))){
# 			par(mar=c(0,0,0,0))
# 			plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 			text(x=0.5,y=0.5,labels=paste('d =',unique(dimensionsToRun)[i]))
# 			par(mar=plotMar)
# 			boxplot(t(c.index.GPSurvSqExp.all[breakDown[i,1]:breakDown[i,2],,drop=FALSE]),at=extraDimensionsToRun[breakDown[i,1]:breakDown[i,2]],xlim=c(min(extraDimensionsToRun),max(extraDimensionsToRun)),ylab='Concordance index',ylim=c(0.4,1.0),xaxt='n')
# 			axis(1,at=unique(extraDimensionsToRun),labels=unique(extraDimensionsToRun))
# 			abline(h=0.5,col='lightgrey',lty=4)
# 		}
# 		par(mar=c(0,0,0,0))
# 		plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 		text(x=0.5,y=0.5,labels='No. uninformative dimensions')
# 		par(mar=c(5,4,4,2)+0.1)
	
# 		layout(matrix(c(1,1,2:(2*length(unique(dimensionsToRun))+1),2*length(unique(dimensionsToRun))+2,2*length(unique(dimensionsToRun))+2),nrow=length(unique(dimensionsToRun))+2,byrow=TRUE),widths=c(1,9),heights=c(1,rep(3,length(unique(dimensionsToRun))),1))
# 		par(mar=c(0,0,0,0))
# 		plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 		text(x=0.5,y=0.5,labels=modelNamesPlot[3],font=2)
# 		par(mar=c(5,4,4,2)+0.1)
# 		for(i in 1:length(unique(dimensionsToRun))){
# 			par(mar=c(0,0,0,0))
# 			plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 			text(x=0.5,y=0.5,labels=paste('d =',unique(dimensionsToRun)[i]))
# 			par(mar=plotMar)
# 			boxplot(t(c.index.GPSurvRSFS.all[breakDown[i,1]:breakDown[i,2],,drop=FALSE]),at=extraDimensionsToRun[breakDown[i,1]:breakDown[i,2]],xlim=c(min(extraDimensionsToRun),max(extraDimensionsToRun)),ylab='Concordance index',ylim=c(0.4,1.0),xaxt='n')
# 			axis(1,at=unique(extraDimensionsToRun),labels=unique(extraDimensionsToRun))
# 			abline(h=0.5,col='lightgrey',lty=4)
# 		}
# 		par(mar=c(0,0,0,0))
# 		plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
# 		text(x=0.5,y=0.5,labels='No. uninformative dimensions')
# 		par(mar=c(5,4,4,2)+0.1)
# 	dev.off(pdf.output)
# }
