#----------------------------------------------------------------------------------------#
# K Lloyd 2017_09_01
#----------------------------------------------------------------------------------------#
# Applying Gaussian process to synthetic data generated with censored survival times in the training set.
# All models applied to the same data. Models are:	GP, GPS1, GPS2, GPS3, RF, RFSurvival and Coxph
# Results saved to folder 'Runs' within working directory.
#---------------------------------------------------------------------------------------#
# Investigating changing noise hyperparameter and training set size
#---------------------------------------------------------------------------------------#


##-------------------------------------------------------------------------------------##
##---------------------------------- Load Libraries -----------------------------------##
##-------------------------------------------------------------------------------------##
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
source('../toSource/GetDataSetFeatures.R')
source('../toSource/ImputeMissingData.R')
source('../toSource/LogPriorX.R')
source('../toSource/MakeSyntheticData.R')
source('../toSource/MeanFunc.R')
source('../toSource/NormaliseExpressionData.R')
source('../toSource/PlotGPLaplaceHypObjFhat.R')
source('../toSource/PlotKaplanMeier.R')
source('../toSource/PreLearnHyperparam.R')
source('../toSource/PrintOptions.R')
source('../toSource/PrintOptions2.R')
source('../toSource/RemoveCensored.R')
source('../toSource/RenameGeneFeaturesAccordingToDataSetCuratedOvarian.R')
source('../toSource/SetParametersCh4Exp3.R')


##-------------------------------------------------------------------------------------##
##------------------------------ Folder & Run Parameters ------------------------------##
##-------------------------------------------------------------------------------------##
set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid())%%2^31))
nReps <- 30


##-------------------------------------------------------------------------------------##
##------------------------------------ Initialise -------------------------------------##
##-------------------------------------------------------------------------------------##
nTrainingToRun 					<- c(100,150,200,350,500,750,1000)
hypGenerateNoiseToRun 			<- c(0.01,0.05,0.10,0.30,0.50,0.70,0.90)
unids 							<- rep(NA,length(nTrainingToRun)*length(hypGenerateNoiseToRun))
count 							<- 1

for(k in 1:length(nTrainingToRun)){
	for(j in 1:length(hypGenerateNoiseToRun)){
		## Make empty vectors and lists ##
		outputStructureGPSurvNoCorr 	<- list()
		c.index.GPSurvNoCorr 			<- rep(NA,nReps)
		rmse.GPSurvNoCorr 				<- rep(NA,nReps)
		outputStructureGPSurvCorrL 		<- list()
		c.index.GPSurvCorrL 			<- rep(NA,nReps)
		rmse.GPSurvCorrL 				<- rep(NA,nReps)
		outputStructureGPSurvCorrV 		<- list()
		c.index.GPSurvCorrV 			<- rep(NA,nReps)
		rmse.GPSurvCorrV 				<- rep(NA,nReps)
		outputStructureRFSurvival 		<- list()
		c.index.RFSurvival 				<- rep(NA,nReps)
		rmse.RFSurvival 				<- rep(NA,nReps)
		outputStructureCoxph 			<- list()
		c.index.Coxph 					<- rep(NA,nReps)
		rmse.Coxph 						<- rep(NA,nReps)
		trainingTestStructureForNReps 	<- list()

		## Generate parameters ##
		allParameterStructures 					<- SetParametersCh4Exp3(nTraining=nTrainingToRun[k],hypGenerateNoise=hypGenerateNoiseToRun[j],censoredProportion=0.75,nReps=nReps)
	
		## Make data ##
		trainingTestStructureForNReps 	 		<- GenerateData(allParameterStructures$dataOptionsStructure,allParameterStructures$plotSaveOptions$outerFolder,nReps)
		for(i in 1:nReps){
			trainingTestStructureForNReps[[i]] 	<- NormaliseExpressionData(trainingTestStructureForNReps[[i]],normaliseFlag=TRUE,winsoriseFlag=FALSE) # normaliseFlag TRUE -> genes normalised to (mean=0,sd=1), winsoriseFlag TRUE -> outside (0.05,0.95) quantiles clipped to here
			cat('Data set for run',i,'normalised',fill=TRUE)
			colnames(trainingTestStructureForNReps[[i]]$trainingData) 	<- paste0('x',1:allParameterStructures$dataOptionsStructure$dimension)
			colnames(trainingTestStructureForNReps[[i]]$testData) 		<- paste0('x',1:allParameterStructures$dataOptionsStructure$dimension)
		}
	
		## Run GPSurvNoCorr ##
		for(i in 1:nReps){
			trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
			outputStructureGPSurvNoCorr[[i]] 	<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
			c.index.GPSurvNoCorr[i] 			<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$c.index)!=0,outputStructureGPSurvNoCorr[[i]]$c.index,NA)
			rmse.GPSurvNoCorr[i] 				<- ifelse(length(outputStructureGPSurvNoCorr[[i]]$rmse)!=0,outputStructureGPSurvNoCorr[[i]]$rmse,NA)
		}
		assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvNoCorr'),outputStructureGPSurvNoCorr)
		save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvNoCorr'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'GPSurvNoCorrWorkspace.RData'))

		## Run GPSurvCorrL ##
		allParameterStructures$parameterStructure$noiseCorr <- 'noiseCorrLearned'
		for(i in 1:nReps){
			trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
			outputStructureGPSurvCorrL[[i]] 	<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
			c.index.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$c.index)!=0,outputStructureGPSurvCorrL[[i]]$c.index,NA)
			rmse.GPSurvCorrL[i] 				<- ifelse(length(outputStructureGPSurvCorrL[[i]]$rmse)!=0,outputStructureGPSurvCorrL[[i]]$rmse,NA)
		}
		assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvCorrL'),outputStructureGPSurvCorrL)
		save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvCorrL'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'GPSurvCorrLWorkspace.RData'))

		## Run GPSurvCorrV ##
		allParameterStructures$parameterStructure$noiseCorr <- 'noiseCorrVec'
		for(i in 1:nReps){
			trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
			outputStructureGPSurvCorrV[[i]] 	<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
			c.index.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$c.index)!=0,outputStructureGPSurvCorrV[[i]]$c.index,NA)
			rmse.GPSurvCorrV[i] 				<- ifelse(length(outputStructureGPSurvCorrV[[i]]$rmse)!=0,outputStructureGPSurvCorrV[[i]]$rmse,NA)
		}
		assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvCorrV'),outputStructureGPSurvCorrV)
		save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureGPSurvCorrV'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'GPSurvCorrVWorkspace.RData'))

		## Run RFSurvival ##
		allParameterStructures$parameterStructure$noiseCorr <- FALSE
		for(i in 1:nReps){
			trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
			outputStructureRFSurvival[[i]] 		<- ApplyRFSurvival(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
			c.index.RFSurvival[i] 				<- ifelse(length(outputStructureRFSurvival[[i]]$c.index)!=0,outputStructureRFSurvival[[i]]$c.index,NA)
			rmse.RFSurvival[i] 					<- ifelse(length(outputStructureRFSurvival[[i]]$rmse)!=0,outputStructureRFSurvival[[i]]$rmse,NA)
		}

		assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureRFSurvival'),outputStructureRFSurvival)
		save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureRFSurvival'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'RFSurvivalWorkspace.RData'))
	
		## Run Coxph ##
		allParameterStructures$parameterStructure$noiseCorr <- FALSE
		for(i in 1:nReps){
			trainingTestStructure 				<- trainingTestStructureForNReps[[i]]
			outputStructureCoxph[[i]] 			<- ApplyCoxph(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
			c.index.Coxph[i] 					<- ifelse(length(outputStructureCoxph[[i]]$c.index)!=0,outputStructureCoxph[[i]]$c.index,NA)
			rmse.Coxph[i] 						<- ifelse(length(outputStructureCoxph[[i]]$rmse)!=0,outputStructureCoxph[[i]]$rmse,NA)
		}

		assign(paste0(allParameterStructures$parameterStructure$unid,'outputStructureCoxph'),outputStructureCoxph)
		save(list=paste0(allParameterStructures$parameterStructure$unid,'outputStructureCoxph'),file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'CoxphWorkspace.RData'))
	
		## Keep unids ##
		unids[count] 	<- allParameterStructures$parameterStructure$unid
		count 			<- count + 1

		save(list='trainingTestStructureForNReps',file=paste0('Runs/',allParameterStructures$parameterStructure$unid,'/',allParameterStructures$parameterStructure$unid,'trainingTestStructureForNReps.RData'))
		
		rm('outputStructureGPSurvNoCorr','c.index.GPSurvNoCorr','rmse.GPSurvNoCorr',
			'outputStructureGPSurvCorrL','c.index.GPSurvCorrL','rmse.GPSurvCorrL',
			'outputStructureGPSurvCorrV','c.index.GPSurvCorrV','rmse.GPSurvCorrV',
			'outputStructureRFSurvival','c.index.RFSurvival','rmse.RFSurvival',
			'outputStructureCoxph','c.index.Coxph','rmse.Coxph',
			'trainingTestStructureForNReps')
	}

}


##-------------------------------------------------------------------------------------##
##----------------------------------- Make Heatmaps -----------------------------------##
##-------------------------------------------------------------------------------------##
c.index.GPS1.all 		<- matrix(0,nrow=length(nTrainingToRun)*length(hypGenerateNoiseToRun),ncol=nReps)
c.index.GPS2.all 		<- matrix(0,nrow=length(nTrainingToRun)*length(hypGenerateNoiseToRun),ncol=nReps)
c.index.GPS3.all 		<- matrix(0,nrow=length(nTrainingToRun)*length(hypGenerateNoiseToRun),ncol=nReps)
c.index.RFSurvival.all 	<- matrix(0,nrow=length(nTrainingToRun)*length(hypGenerateNoiseToRun),ncol=nReps)
c.index.Coxph.all 		<- matrix(0,nrow=length(nTrainingToRun)*length(hypGenerateNoiseToRun),ncol=nReps)
c.index.GPS1.mean 		<- numeric()
c.index.GPS2.mean 		<- numeric()
c.index.GPS3.mean 		<- numeric()
c.index.RFSurvival.mean <- numeric()
c.index.Coxph.mean 		<- numeric()
for(i in 1:(length(nTrainingToRun)*length(hypGenerateNoiseToRun))){
	for(j in 1:nReps){
		c.index.GPS1.all[i,j] 		<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPSurvNoCorr[[',j,']]$c.index')))
		c.index.GPS2.all[i,j] 		<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPSurvCorrL[[',j,']]$c.index')))
		c.index.GPS3.all[i,j] 		<- 1-eval(parse(text=paste0(unids[i],'outputStructureGPSurvCorrV[[',j,']]$c.index')))
	
		c.index.RFSurvival.all[i,j] <- eval(parse(text=paste0(unids[i],'outputStructureRFSurvival[[',j,']]$c.index')))
		c.index.Coxph.all[i,j] 		<- eval(parse(text=paste0(unids[i],'outputStructureCoxph[[',j,']]$c.index')))
	}
	c.index.GPS1.mean[i] 		<- mean(c.index.GPS1.all[i,],na.rm=TRUE)
	c.index.GPS2.mean[i] 		<- mean(c.index.GPS2.all[i,],na.rm=TRUE)
	c.index.GPS3.mean[i] 		<- mean(c.index.GPS3.all[i,],na.rm=TRUE)

	c.index.RFSurvival.mean[i] 	<- mean(c.index.RFSurvival.all[i,],na.rm=TRUE)
	c.index.Coxph.mean[i] 		<- mean(c.index.Coxph.all[i,],na.rm=TRUE)
}

c.index.GPS1.mean 					<- matrix(c.index.GPS1.mean,nrow=length(nTrainingToRun))
c.index.GPS2.mean 					<- matrix(c.index.GPS2.mean,nrow=length(nTrainingToRun))
c.index.GPS3.mean 					<- matrix(c.index.GPS3.mean,nrow=length(nTrainingToRun))
c.index.RFSurvival.mean 			<- matrix(c.index.RFSurvival.mean,nrow=length(nTrainingToRun))
c.index.Coxph.mean 					<- matrix(c.index.Coxph.mean,nrow=length(nTrainingToRun))
colnames(c.index.GPS1.mean) 		<- nTrainingToRun
colnames(c.index.GPS2.mean) 		<- nTrainingToRun
colnames(c.index.GPS3.mean) 		<- nTrainingToRun
colnames(c.index.RFSurvival.mean) 	<- nTrainingToRun
colnames(c.index.Coxph.mean) 		<- nTrainingToRun
rownames(c.index.GPS1.mean) 		<- hypGenerateNoiseToRun
rownames(c.index.GPS2.mean) 		<- hypGenerateNoiseToRun
rownames(c.index.GPS3.mean) 		<- hypGenerateNoiseToRun
rownames(c.index.RFSurvival.mean) 	<- hypGenerateNoiseToRun
rownames(c.index.Coxph.mean) 		<- hypGenerateNoiseToRun
grd 								<- expand.grid(x=hypGenerateNoiseToRun, y=nTrainingToRun)

pdf(paste0('Runs/','GPSurvNoCorr','PlotCIndexMatrix.pdf'),width=10,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	oldpar <- par(no.readonly=TRUE)
	par(mar=c(5.5,4.5,2.5,6.0))
	image(x=nTrainingToRun,y=hypGenerateNoiseToRun,z=t(c.index.GPS1.mean),xlab='# Training Samples',ylab=expression(sigma[n]^{2}),col=colorRampPalette(c('white','chartreuse4'))(16),
		xaxt='n',yaxt='n',zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)),cex.lab=1.2)
	axis(1,at=nTrainingToRun,labels=nTrainingToRun,cex.axis=1.2)
	axis(2,at=hypGenerateNoiseToRun,labels=hypGenerateNoiseToRun,las=2,cex.axis=1.2)
	image.plot(y=hypGenerateNoiseToRun,x=nTrainingToRun,z=t(c.index.GPS1.mean),add=T,legend.only=TRUE,legend.lab='Mean Concordance Index',col=colorRampPalette(c('white','chartreuse4'))(16),
		legend.line=3.1,legend.shrink=0.6,legend.mar=5.8,zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)))
	text(rep(nTrainingToRun,length(nTrainingToRun)),rep(hypGenerateNoiseToRun,each=length(hypGenerateNoiseToRun)),sprintf('%.2f',round(t(c.index.GPS1.mean),2)),offset=0,cex=0.8)
	par(oldpar)
dev.off(pdf.output)

pdf(paste0('Runs/','GPSurvCorrL','PlotCIndexMatrix.pdf'),width=10,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	oldpar <- par(no.readonly=TRUE)
	par(mar=c(5.5,4.5,2.5,6.0))
	image(x=nTrainingToRun,y=hypGenerateNoiseToRun,z=t(c.index.GPS2.mean),xlab='# Training Samples',ylab=expression(sigma[n]^{2}),col=colorRampPalette(c('white','chartreuse4'))(16),
				xaxt='n',yaxt='n',zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)),cex.lab=1.2)
	axis(1,at=nTrainingToRun,labels=nTrainingToRun,cex.axis=1.2)
	axis(2,at=hypGenerateNoiseToRun,labels=hypGenerateNoiseToRun,las=2,cex.axis=1.2)
	image.plot(y=hypGenerateNoiseToRun,x=nTrainingToRun,z=t(c.index.GPS2.mean),add=T,legend.only=TRUE,legend.lab='Mean Concordance Index',col=colorRampPalette(c('white','chartreuse4'))(16),
		legend.line=3.1,legend.shrink=0.6,legend.mar=5.8,zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)))
	text(rep(nTrainingToRun,length(nTrainingToRun)),rep(hypGenerateNoiseToRun,each=length(hypGenerateNoiseToRun)),sprintf('%.2f',round(t(c.index.GPS2.mean),2)),offset=0,cex=0.8)
	par(oldpar)
dev.off(pdf.output)

pdf(paste0('Runs/','GPSurvCorrV','PlotCIndexMatrix.pdf'),width=10,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	oldpar <- par(no.readonly=TRUE)
	par(mar=c(5.5,4.5,2.5,6.0))
	image(x=nTrainingToRun,y=hypGenerateNoiseToRun,z=t(c.index.GPS3.mean),xlab='# Training Samples',ylab=expression(sigma[n]^{2}),col=colorRampPalette(c('white','chartreuse4'))(16),
				xaxt='n',yaxt='n',zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)),cex.lab=1.2)
	axis(1,at=nTrainingToRun,labels=nTrainingToRun,cex.axis=1.2)
	axis(2,at=hypGenerateNoiseToRun,labels=hypGenerateNoiseToRun,las=2,cex.axis=1.2)
	image.plot(y=hypGenerateNoiseToRun,x=nTrainingToRun,z=t(c.index.GPS3.mean),add=T,legend.only=TRUE,legend.lab='Mean Concordance Index',col=colorRampPalette(c('white','chartreuse4'))(16),
		legend.line=3.1,legend.shrink=0.6,legend.mar=5.8,zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)))
	text(rep(nTrainingToRun,length(nTrainingToRun)),rep(hypGenerateNoiseToRun,each=length(hypGenerateNoiseToRun)),sprintf('%.2f',round(t(c.index.GPS3.mean),2)),offset=0,cex=0.8)
	par(oldpar)
dev.off(pdf.output)

pdf(paste0('Runs/','RFSurvival','PlotCIndexMatrix.pdf'),width=10,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	oldpar <- par(no.readonly=TRUE)
	par(mar=c(5.5,4.5,2.5,6.0))
	image(x=nTrainingToRun,y=hypGenerateNoiseToRun,z=t(c.index.RFSurvival.mean),xlab='# Training Samples',ylab=expression(sigma[n]^{2}),col=colorRampPalette(c('white','chartreuse4'))(16),
				xaxt='n',yaxt='n',zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)),cex.lab=1.2)
	axis(1,at=nTrainingToRun,labels=nTrainingToRun,cex.axis=1.2)
	axis(2,at=hypGenerateNoiseToRun,labels=hypGenerateNoiseToRun,las=2,cex.axis=1.2)
	image.plot(y=hypGenerateNoiseToRun,x=nTrainingToRun,z=t(c.index.RFSurvival.mean),add=T,legend.only=TRUE,legend.lab='Mean Concordance Index',col=colorRampPalette(c('white','chartreuse4'))(16),
		legend.line=3.1,legend.shrink=0.6,legend.mar=5.8,zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)))
	text(rep(nTrainingToRun,length(nTrainingToRun)),rep(hypGenerateNoiseToRun,each=length(hypGenerateNoiseToRun)),sprintf('%.2f',round(t(c.index.RFSurvival.mean),2)),offset=0,cex=0.8)
	par(oldpar)
dev.off(pdf.output)

pdf(paste0('Runs/','Coxph','PlotCIndexMatrix.pdf'),width=10,height=7,onefile=TRUE)
pdf.output <- dev.cur()
	oldpar <- par(no.readonly=TRUE)
	par(mar=c(5.5,4.5,2.5,6.0))
	image(x=nTrainingToRun,y=hypGenerateNoiseToRun,z=t(c.index.Coxph.mean),xlab='# Training Samples',ylab=expression(sigma[n]^{2}),col=colorRampPalette(c('white','chartreuse4'))(16)
		,		xaxt='n',yaxt='n',zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)),cex.lab=1.2)
	axis(1,at=nTrainingToRun,labels=nTrainingToRun,cex.axis=1.2)
	axis(2,at=hypGenerateNoiseToRun,labels=hypGenerateNoiseToRun,las=2,cex.axis=1.2)
	image.plot(y=hypGenerateNoiseToRun,x=nTrainingToRun,z=t(c.index.Coxph.mean),add=T,legend.only=TRUE,legend.lab='Mean Concordance Index',col=colorRampPalette(c('white','chartreuse4'))(16),
		legend.line=3.1,legend.shrink=0.6,legend.mar=5.8,zlim=range(c(c.index.GPS1.mean,c.index.GPS2.mean,c.index.GPS3.mean,c.index.Coxph.mean,c.index.RFSurvival.mean)))
	text(rep(nTrainingToRun,length(nTrainingToRun)),rep(hypGenerateNoiseToRun,each=length(hypGenerateNoiseToRun)),sprintf('%.2f',round(t(c.index.Coxph.mean),2)),offset=0,cex=0.8)
	par(oldpar)
dev.off(pdf.output)