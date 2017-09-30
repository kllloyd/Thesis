#----------------------------------------------------------------------------------------#
# K Lloyd 2016_08_09
#----------------------------------------------------------------------------------------#
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# Models investigating feature selection
# All models applied to the same data. Models are:	GP, GPS3SqExp and GPS3SqExpRSFS
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Compare c index results of different models
# Number of bootstraps/feature subsets applied for GPS3SqExpRSFS varied
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
source('../toSource/MakeEnsemblePlots.R')


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

nReps 					<- 100
nBootstrapsToRun 		<- c(5,25,50,100,200,400,600,800,1000)
nBootstraps 			<- 1000

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
	logHypGenerate 		<- list('noise'=log(0.05),'func'=log(0.6),'length'=log(c(1.5,2,3,2.5,4,5,4.5,5.5,6,3,6.5,12,11,13.5,10)),'mean'=c(rep(0,dimension),0))
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
outputStructureGPSurvSqExpRSFS 		<- rep(list(list()),nReps)
c.index.GPNonSurvNoCens 			<- rep(NA,nReps)
c.index.GPSurvSqExp 				<- rep(NA,nReps)
c.index.GPSurvSqExpRSFS 			<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
rmse.GPNonSurvNoCens 				<- rep(NA,nReps)
rmse.GPSurvSqExp 					<- rep(NA,nReps)
rmse.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
brier.GPNonSurvNoCens 				<- rep(NA,nReps)
brier.GPSurvSqExp 					<- rep(NA,nReps)
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
parameterStructure$logHypStart 			<- list('noise'=log(0.2),'func'=log(0.8),'length'=log(2),'mean'=c(rep(0,dimension),0))
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
##------------------------------- Run GPS3 SqExp Model  -------------------------------##
##-------------------------------------------------------------------------------------##
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


##-------------------------------------------------------------------------------------##
##----------- Run GPS3 SqExp Model Applied to Bootstrapped Feature Subsets ------------##
##-------------------------------------------------------------------------------------##
subsetDimension 					<- 3 		# e.g. c(3,4)
toCount 							<- 'features'
bic.GPSurvSqExpRSFS 				<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
ensembleStructure 					<- list()
subsetIndices 						<- matrix(rep(NA,nBootstraps*nReps),ncol=nBootstraps)
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
	ensembleStructure.1000[[i]] <- CalculateEnsembleResults(outputStructureGPSurvSqExpRSFS[[i]],trainingTestStructureForNReps[[i]],bic.GPSurvSqExpRSFS[i,],c.index.GPSurvSqExpRSFS[i,])
}
file.rename(from=paste0('Runs/',unid,'/GPSurvCorrV'),to=paste0('Runs/',unid,'/GPSurvSqExpRSFS'))
save(list='outputStructureGPSurvSqExpRSFS',file=paste0('Runs','/',unid,'/',unid,'_','outputStructureGPSurvSqExpRSFS','_','Workspace.RData'))


##-------------------------------------------------------------------------------------##
##----------------- Sample from 1000 bootstraps to use other numbers ------------------##
##-------------------------------------------------------------------------------------##

# Results with GPNonSurvNoCens and GPSurvSqExp #
c.index.GPSurvSqExp 					<- sapply(1:nReps,function(x) outputStructureFullGPSurvSqExp[[x]]$c.index)
c.index.GPNonSurvNoCens 				<- sapply(1:nReps,function(x) outputStructureFullGPNonSurvNoCens[[x]]$c.index)

# Results with GPSurvSqExpRSFS, 1000 bootstraps #
c.index.GPSurvSqExpRSFS.1000.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.1000[[x]]$ensembleMetrics$c.index)
c.index.GPSurvSqExpRSFS.1000.uniform 	<- sapply(1:nReps,function(x) ensembleStructure.1000[[x]]$ensembleMetricsUnifWeight$c.index)

# Results with GPSurvSqExpRSFS, 800 bootstraps #
structureRSFS.800 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=800,trainingTestStructureForNReps)
ensembleStructure.800 					<- structureRSFS.800$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.800.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.800[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 600 bootstraps #
structureRSFS.600 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=750,trainingTestStructureForNReps)
ensembleStructure.600 					<- structureRSFS.600$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.600.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.600[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 400 bootstraps #
structureRSFS.400 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=400,trainingTestStructureForNReps)
ensembleStructure.400 					<- structureRSFS.400$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.400.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.400[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 200 bootstraps #
structureRSFS.200 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=200,trainingTestStructureForNReps)
ensembleStructure.200 					<- structureRSFS.200$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.200.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.200[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 100 bootstraps #
structureRSFS.100 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=100,trainingTestStructureForNReps)
ensembleStructure.100 					<- structureRSFS.100$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.100.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.100[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 50 bootstraps #
structureRSFS.50 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=50,trainingTestStructureForNReps)
ensembleStructure.50 					<- structureRSFS.50$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.50.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.50[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 25 bootstraps #
structureRSFS.25 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=25,trainingTestStructureForNReps)
ensembleStructure.25 					<- structureRSFS.25$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.25.ensemble 	<- sapply(1:nReps,function(x) ensembleStructure.25[[x]]$ensembleMetrics$c.index)

# Results with GPSurvSqExpRSFS, 5 bootstraps #
structureRSFS.5 						<- SampleRSFSBootstraps(unids,nBootstrapsFull=1000,nBootstraps=5,trainingTestStructureForNReps)
ensembleStructure.5 					<- structureRSFS.5$ensembleStructure.sampled
c.index.GPSurvSqExpRSFS.5.ensemble 		<- sapply(1:nReps,function(x) ensembleStructure.5[[x]]$ensembleMetrics$c.index)


##-------------------------------------------------------------------------------------##
##------------------------------ Plot concordance index -------------------------------##
##-------------------------------------------------------------------------------------##
c.index.mat <- cbind(c.index.GPNonSurvNoCens,c.index.GPSurvSqExp,
					c.index.GPSurvSqExpRSFS.5.ensemble,c.index.GPSurvSqExpRSFS.25.ensemble,c.index.GPSurvSqExpRSFS.50.ensemble,
					c.index.GPSurvSqExpRSFS.100.ensemble,c.index.GPSurvSqExpRSFS.200.ensemble,c.index.GPSurvSqExpRSFS.400.ensemble,c.index.GPSurvSqExpRSFS.600.ensemble,
					c.index.GPSurvSqExpRSFS.800.ensemble,c.index.GPSurvSqExpRSFS.1000.ensemble)
nBootstrapsToRun <- c(5,25,50,100,200,400,600,800,1000)

							##--------------------------------##
							##-------- Plot with RSFS --------##
							##----------- Ch5 Fig 4 ----------##
							##--------------------------------##
oldMar <- par('mar')
pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotCIndexRSFSNBoostraps.pdf'),width=10,height=6)
	pdf.output <- dev.cur()
	layout(matrix(c(1,2,3,4),ncol=2,byrow=TRUE),widths=c(1,7),heights=c(9,1))
	par(mar=c(5.1,4.1,4.1,0.1))
	boxplot(1-c.index.mat[,c(1,2)],at=c(2,4),xaxt='n',col=c('green4','chartreuse4'),ylim=c(0.4,1),xlim=c(1,5),axes=FALSE)
	axis(1,at=c(2,4),labels=c('GP','GPS3SqExp'),las=2)
	axis(2,at=seq(from=0.4,to=1.0,length.out=7),labels=seq(from=0.4,to=1.0,length.out=7),las=1)
	title(ylab='Concordance Index')
	vps1 <- do.call(vpStack, baseViewports())
	par(mar=c(5.1,0.1,4.1,2.1))
	boxplot(1-c.index.mat[,3:11],at=c(nBootstrapsToRun/10),xaxt='n',ylim=c(0.4,1),yaxt='n',axes=FALSE,col='white')
	axis(1,at=c(nBootstrapsToRun/10),labels=c(nBootstrapsToRun),las=1)
	title(xlab='GPS3SqExpRSFS, Number of bootstraps')
	lines(nBootstrapsToRun/10,c(median(1-c.index.GPSurvSqExpRSFS.5.ensemble),median(1-c.index.GPSurvSqExpRSFS.25.ensemble),median(1-c.index.GPSurvSqExpRSFS.50.ensemble),
	                            median(1-c.index.GPSurvSqExpRSFS.100.ensemble),median(1-c.index.GPSurvSqExpRSFS.200.ensemble),median(1-c.index.GPSurvSqExpRSFS.400.ensemble),median(1-c.index.GPSurvSqExpRSFS.600.ensemble),
	                            median(1-c.index.GPSurvSqExpRSFS.800.ensemble),median(1-c.index.GPSurvSqExpRSFS.1000.ensemble)),col='chartreuse3',lwd=1.5)
	boxplot(1-c.index.mat[,3:11],at=c(nBootstrapsToRun/10),xaxt='n',col=c(rep('chartreuse3',10)),ylim=c(0.4,1),yaxt='n',axes=FALSE,add=TRUE)
	vps3 <- do.call(vpStack, baseViewports())
	pushViewport(vps1)
	Y1 <- convertY(unit(1.024,"native"), "npc")
	Y2 <- convertY(unit(0.4-0.024,"native"), "npc")
	popViewport(3)
	grid.move.to(x = unit(0, "npc"), y = Y1, vp = vps1)
	grid.line.to(x = unit(1, "npc"), y = Y1, vp = vps3, gp = gpar(col = "black"))
	grid.move.to(x = unit(1.4, "npc"), y = Y1, vp = vps1)
	grid.line.to(x = unit(0, "npc"), y = Y1, vp = vps3, gp = gpar(col = "white",lwd=2))
	grid.move.to(x = unit(0, "npc"), y = Y2, vp = vps1)
	grid.line.to(x = unit(1, "npc"), y = Y2, vp = vps3, gp = gpar(col = "black"))
	grid.move.to(x = unit(1.4, "npc"), y = Y2, vp = vps1)
	grid.line.to(x = unit(0, "npc"), y = Y2, vp = vps3, gp = gpar(col = "white",lwd=2))
	grid.move.to(x = unit(0.85,"native"), y = unit(0, "npc"), vp = vps1)
	grid.line.to(x = unit(0.85,"native"), y = unit(1, "npc"), vp = vps1, gp = gpar(col = "black"))
	grid.move.to(x = unit(104.5,"native"), y = unit(0, "npc"), vp = vps3)
	grid.line.to(x = unit(104.5,"native"), y = unit(1, "npc"), vp = vps3, gp = gpar(col = "black"))
	layout(1)
	par(mar=oldMar)
dev.off(pdf.output)


							##--------------------------------##
							##-------- Plot weightings -------##
							##--------- Ch5 AppFig 4 ---------##
							##--------------------------------##
source('../toSource/add.alpha.R')
pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotWeightings.pdf'),width=6,height=6)
	pdf.output <- dev.cur()
	boxplot(sapply(1:10,function(x) ensembleStructure.5[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.25[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.50[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.100[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.200[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.400[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.600[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.800[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
	boxplot(sapply(1:10,function(x) ensembleStructure.1000[[x]]$weightings),ylab='Model weightings',xlab='Repeat',col=add.alpha('skyblue',0.5))
dev.off(pdf.output)

							##--------------------------------##
							##-- Plot with uniform weighting -##
							##--------- Ch5 AppFig 5 ---------##
							##--------------------------------##
c.index.GPSurvSqExpRSFS.800.uniform <- sapply(1:nReps,function(x) ensembleStructure.800[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.600.uniform <- sapply(1:nReps,function(x) ensembleStructure.600[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.400.uniform <- sapply(1:nReps,function(x) ensembleStructure.400[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.200.uniform <- sapply(1:nReps,function(x) ensembleStructure.200[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.100.uniform <- sapply(1:nReps,function(x) ensembleStructure.100[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.50.uniform 	<- sapply(1:nReps,function(x) ensembleStructure.50[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.25.uniform 	<- sapply(1:nReps,function(x) ensembleStructure.25[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.10.uniform 	<- sapply(1:nReps,function(x) ensembleStructure.10[[x]]$ensembleMetricsUnifWeight$c.index)
c.index.GPSurvSqExpRSFS.5.uniform 	<- sapply(1:nReps,function(x) ensembleStructure.5[[x]]$ensembleMetricsUnifWeight$c.index)

pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotCIndexRSFSUniformNBoostraps.pdf'),width=10,height=6)
pdf.output <- dev.cur()
	layout(rbind(1,2), heights=c(8,1))
	boxplot(cbind(1-c.index.GPSurvSqExpRSFS.5.ensemble,1-c.index.GPSurvSqExpRSFS.5.uniform,
				1-c.index.GPSurvSqExpRSFS.25.ensemble,1-c.index.GPSurvSqExpRSFS.25.uniform,
				1-c.index.GPSurvSqExpRSFS.50.ensemble,1-c.index.GPSurvSqExpRSFS.50.uniform,
				1-c.index.GPSurvSqExpRSFS.100.ensemble,1-c.index.GPSurvSqExpRSFS.100.uniform,
				1-c.index.GPSurvSqExpRSFS.200.ensemble,1-c.index.GPSurvSqExpRSFS.200.uniform,
				1-c.index.GPSurvSqExpRSFS.400.ensemble,1-c.index.GPSurvSqExpRSFS.400.uniform,
				1-c.index.GPSurvSqExpRSFS.600.ensemble,1-c.index.GPSurvSqExpRSFS.600.uniform,
				1-c.index.GPSurvSqExpRSFS.800.ensemble,1-c.index.GPSurvSqExpRSFS.800.uniform,
				1-c.index.GPSurvSqExpRSFS.1000.ensemble,1-c.index.GPSurvSqExpRSFS.1000.uniform),
				at=c(5-5,5+5,25-5,25+5,50-5,50+5,100-5,100+5,200-5,200+5,400-5,400+5,600-5,600+5,800-5,800+5,1000-5,1000+5)/10,
				xaxt='n',ylim=c(0.4,1),col='white',axes=FALSE)
	lines((nBootstrapsToRun-5)/10,c(median(1-c.index.GPSurvSqExpRSFS.5.ensemble),median(1-c.index.GPSurvSqExpRSFS.25.ensemble),median(1-c.index.GPSurvSqExpRSFS.50.ensemble),
		                            median(1-c.index.GPSurvSqExpRSFS.100.ensemble),median(1-c.index.GPSurvSqExpRSFS.200.ensemble),median(1-c.index.GPSurvSqExpRSFS.400.ensemble),median(1-c.index.GPSurvSqExpRSFS.600.ensemble),
		                            median(1-c.index.GPSurvSqExpRSFS.800.ensemble),median(1-c.index.GPSurvSqExpRSFS.1000.ensemble)),col='chartreuse3',lwd=1.5)
	lines((nBootstrapsToRun+5)/10,c(median(1-c.index.GPSurvSqExpRSFS.5.uniform),median(1-c.index.GPSurvSqExpRSFS.25.uniform),median(1-c.index.GPSurvSqExpRSFS.50.uniform),
		                            median(1-c.index.GPSurvSqExpRSFS.100.uniform),median(1-c.index.GPSurvSqExpRSFS.200.uniform),median(1-c.index.GPSurvSqExpRSFS.400.uniform),median(1-c.index.GPSurvSqExpRSFS.600.uniform),
		                            median(1-c.index.GPSurvSqExpRSFS.800.uniform),median(1-c.index.GPSurvSqExpRSFS.1000.uniform)),col='chartreuse4',lwd=1.5)
	boxplot(cbind(1-c.index.GPSurvSqExpRSFS.5.ensemble,1-c.index.GPSurvSqExpRSFS.5.uniform,
				1-c.index.GPSurvSqExpRSFS.25.ensemble,1-c.index.GPSurvSqExpRSFS.25.uniform,
				1-c.index.GPSurvSqExpRSFS.50.ensemble,1-c.index.GPSurvSqExpRSFS.50.uniform,
				1-c.index.GPSurvSqExpRSFS.100.ensemble,1-c.index.GPSurvSqExpRSFS.100.uniform,
				1-c.index.GPSurvSqExpRSFS.200.ensemble,1-c.index.GPSurvSqExpRSFS.200.uniform,
				1-c.index.GPSurvSqExpRSFS.400.ensemble,1-c.index.GPSurvSqExpRSFS.400.uniform,
				1-c.index.GPSurvSqExpRSFS.600.ensemble,1-c.index.GPSurvSqExpRSFS.600.uniform,
				1-c.index.GPSurvSqExpRSFS.800.ensemble,1-c.index.GPSurvSqExpRSFS.800.uniform,
				1-c.index.GPSurvSqExpRSFS.1000.ensemble,1-c.index.GPSurvSqExpRSFS.1000.uniform),
				at=c(5-5,5+5,25-5,25+5,50-5,50+5,100-5,100+5,200-5,200+5,400-5,400+5,600-5,600+5,800-5,800+5,1000-5,1000+5)/10,
				xaxt='n',col=c(rep(c('chartreuse3','chartreuse4'),9)),ylim=c(0.4,1),add=TRUE)
	axis(1,at=c(nBootstrapsToRun/10),labels=c(nBootstrapsToRun),las=1)
	title(xlab='Number of bootstraps')
	par(mar=c(0,0,0,0))
	plot.new()
	legend('center',c('GPS3SqExp exp(-BIC) weighting','GPS3SqExp uniform weighting'),col=c('chartreuse3','chartreuse4'),pch=c(15,15),lty=c(1,1),lwd=c(1.5,1.5),ncol=2,bty ="n")
	par(mar=oldMar)
	layout(1)
dev.off(pdf.output)


##-------------------------------------------------------------------------------------##
##-------------------------------- Compute statistics ---------------------------------##
##-------------------------------------------------------------------------------------##
							##--------------------------------##
							##-------- Plot with RSFS --------##
							##----------- Ch5 Fig 4 ----------##
							##--------------------------------##
wilcoxTestOut_GP_GPS3SqExp 		<- list()
wilcoxTestOut_GP_RSFS 			<- list()
wilcoxTestOut_GPS3SqExp_RSFS 	<- list()

WStatistic_GP_GPS3SqExp 	<- numeric(1)
CI.95.1_GP_GPS3SqExp 		<- numeric(1)
CI.95.2_GP_GPS3SqExp 		<- numeric(1)
p.value_GP_GPS3SqExp 		<- numeric(1)
cohens.D_GP_GPS3SqExp 		<- numeric(1)

WStatistic_GP_RSFS 			<- numeric(9)
CI.95.1_GP_RSFS 			<- numeric(9)
CI.95.2_GP_RSFS 			<- numeric(9)
p.value_GP_RSFS 			<- numeric(9)
cohens.D_GP_RSFS 			<- numeric(9)

WStatistic_GPS3SqExp_RSFS 	<- numeric(9)
CI.95.1_GPS3SqExp_RSFS 		<- numeric(9)
CI.95.2_GPS3SqExp_RSFS 		<- numeric(9)
p.value_GPS3SqExp_RSFS 		<- numeric(9)
cohens.D_GPS3SqExp_RSFS 	<- numeric(9)

wilcoxTestOut_GP_GPS3SqExp 	<- wilcox.test(c.index.GPNonSurvNoCens,c.index.GPSurvSqExp,alternative='greater',conf.int=TRUE)
WStatistic_GP_GPS3SqExp 	<- wilcoxTestOut_GP_GPS3SqExp$statistic
CI.95.1_GP_GPS3SqExp 		<- wilcoxTestOut_GP_GPS3SqExp$conf.int[1]
CI.95.2_GP_GPS3SqExp 		<- wilcoxTestOut_GP_GPS3SqExp$conf.int[2]
p.value_GP_GPS3SqExp 		<- wilcoxTestOut_GP_GPS3SqExp$p.value
cohens.D_GP_GPS3SqExp 		<- cohensD(c.index.GPNonSurvNoCens,c.index.GPSurvSqExp)

modelNames 					<- c('GPSurvSqExpRSFS.5.ensemble','GPSurvSqExpRSFS.25.ensemble','GPSurvSqExpRSFS.50.ensemble','GPSurvSqExpRSFS.100.ensemble',
								'GPSurvSqExpRSFS.200.ensemble','GPSurvSqExpRSFS.400.ensemble','GPSurvSqExpRSFS.600.ensemble','GPSurvSqExpRSFS.800.ensemble',
								'GPSurvSqExpRSFS.1000.ensemble')
for(i in 1:length(modelNames)){
	wilcoxTestOut_GP_RSFS[[i]] 	<- wilcox.test(c.index.GPNonSurvNoCens,eval(parse(text=paste0('c.index.',modelNames[i]))),alternative='greater',conf.int=TRUE)
	WStatistic_GP_RSFS[i] 		<- wilcoxTestOut_GP_RSFS[[i]]$statistic
	CI.95.1_GP_RSFS[i] 			<- wilcoxTestOut_GP_RSFS[[i]]$conf.int[1]
	CI.95.2_GP_RSFS[i] 			<- wilcoxTestOut_GP_RSFS[[i]]$conf.int[2]
	p.value_GP_RSFS[i] 			<- wilcoxTestOut_GP_RSFS[[i]]$p.value
	cohens.D_GP_RSFS[i] 		<- cohensD(c.index.GPNonSurvNoCens,eval(parse(text=paste0('c.index.',modelNames[i]))))
}
for(i in 1:length(modelNames)){
	wilcoxTestOut_GPS3SqExp_RSFS[[i]] 	<- wilcox.test(c.index.GPSurvSqExp,eval(parse(text=paste0('c.index.',modelNames[i]))),alternative='greater',conf.int=TRUE)
	WStatistic_GPS3SqExp_RSFS[i] 		<- wilcoxTestOut_GPS3SqExp_RSFS[[i]]$statistic
	CI.95.1_GPS3SqExp_RSFS[i] 			<- wilcoxTestOut_GPS3SqExp_RSFS[[i]]$conf.int[1]
	CI.95.2_GPS3SqExp_RSFS[i] 			<- wilcoxTestOut_GPS3SqExp_RSFS[[i]]$conf.int[2]
	p.value_GPS3SqExp_RSFS[i] 			<- wilcoxTestOut_GPS3SqExp_RSFS[[i]]$p.value
	cohens.D_GPS3SqExp_RSFS[i] 			<- cohensD(c.index.GPSurvSqExp,eval(parse(text=paste0('c.index.',modelNames[i]))))
}

p.value 		<- c(p.value_GP_GPS3SqExp,p.value_GP_RSFS,p.value_GPS3SqExp_RSFS)
p.value.Holm 	<- p.adjust(p.value,'holm')