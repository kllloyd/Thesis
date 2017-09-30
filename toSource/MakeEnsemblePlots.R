MakeEnsemblePlots <- function(outputStructure,trainingTestStructure,bic,c.index,ensembleStructure){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Generating plots relevant to ensemble predictions
	# Requires UpSetR package
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	source('../toSource/PlotKaplanMeier.R')

	# Extract results for plotting #
	unid 						<- outputStructure[[1]][[1]]$parameterStructure$unid
	ensemblePredictions 		<- lapply(ensembleStructure,'[[',4)
	ensembleMetrics 			<- lapply(ensembleStructure,'[[',5)
	nReps 						<- length(trainingTestStructure)
	nBootstraps 				<- length(outputStructure[[1]])
	features 					<- colnames(trainingTestStructure[[1]]$trainingData)
	dimension 					<- length(features)
	featuresPerBootstrap 		<- rep(list(list()),nReps)
	featuresBinary 				<- rep(list(matrix(rep(0,nBootstraps*dimension),nrow=nBootstraps)),nReps)
	featuresBinary 				<- lapply(featuresBinary, function(x) {colnames(x) <- features;x})
	resultsAll 					<- list()
	resultsUnique 				<- list()
	marginalisedProb 			<- matrix(rep(NA,nReps*dimension),ncol=dimension)
	for(i in 1:nReps){
		for(j in 1:nBootstraps){
			featuresPerBootstrap[[i]][[j]] <- colnames(outputStructure[[i]][[j]]$trainingTestStructure$trainingData)
			for(k in 1:dimension){
				if(features[k]%in%featuresPerBootstrap[[i]][[j]]) featuresBinary[[i]][j,k] <- 1
			}
		}
		resultsAll[[i]] 			<- as.data.frame(cbind(featuresBinary[[i]],bic[i,],c.index[i,]))
		colnames(resultsAll[[i]]) 	<- c(paste0('x',1:dimension),'bic','c.index')
		resultsAll[[i]]$c.index 	<- 1-resultsAll[[i]]$c.index
		resultsUnique[[i]] 			<- resultsAll[[i]]
		# resultsUnique[[i]] 			<- resultsAll[[i]][!duplicated(resultsAll[[i]][,1:10]),]
	}
	for(i in 1:nReps){	
		for(k in 1:dimension){
			forMarginalising 		<- exp(-resultsUnique[[i]]$bic)
			forMarginalising[is.infinite(forMarginalising)&forMarginalising>0] <- 10^300
			forMarginalising[is.infinite(forMarginalising)&forMarginalising<0] <- -10^300
			marginalisedProb[i,k] 	<- forMarginalising%*%resultsUnique[[i]][,k]
		}
	}

	# Plot summary of bootstraps #
	library(UpSetR)
	pdf(paste0('Runs/',unid,'/',unid,'PlotAllFeatureSubsets.pdf'),width=20,height=10,onefile=TRUE)
	pdf.output <- dev.cur()
	for(i in 1:nReps){
		upset(resultsAll[[i]],nsets=dimension,nintersects=NA,sets=rev(paste0('x',1:dimension)),keep.order=TRUE,boxplot.summary=c('bic','c.index'))
	}
	dev.off(pdf.output)

	# Plot BIC and C Index #
	pdf(paste0('Runs/',unid,'/',unid,'PlotHistograms.pdf'),width=12,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
	for(i in 1:nReps){
		layout(matrix(c(1,2),ncol=2),widths=c(1,1))
		hist(resultsUnique[[i]]$c.index,main='Histogram of concordance index',xlab='Concordance index')
		# abline(v=c(0.55,0.77),col='chartreuse3',lty=4)
		hist(resultsUnique[[i]]$bic,main='Histogram of BIC',xlab='BIC')
		# abline(v=c(-250,200),col='chartreuse3',lty=4)
	}
	dev.off(pdf.output)

	# Plot BIC vs C Index #
	pdf(paste0('Runs/',unid,'/',unid,'PlotBICcIndex.pdf'),width=12,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
	for(i in 1:nReps){
		plot(resultsUnique[[i]]$bic,resultsUnique[[i]]$c.index,xlab='BIC',ylab='C Index')
	}
	dev.off(pdf.output)

	# Plot BIC marginalised over features #
	pdf(paste0('Runs/',unid,'/',unid,'PlotMarginalisedBIC.pdf'),width=12,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
	for(i in 1:nReps){
		plot(marginalisedProb[i,],ylab='p(D|M,F)',xaxt='n')
		axis(1,at=c(1:dimension),labels=paste0('x',1:dimension))
	}
	dev.off(pdf.output)

	# Plot ensemble predictions vs test targets #
	pdf(paste0('Runs/',unid,'/',unid,'PlotEnsemblePredictions.pdf'),width=12,height=8,onefile=TRUE)
	pdf.output <- dev.cur()
	for(i in 1:nReps){
		plot(trainingTestStructure[[i]]$testTargetsPreCensoring,ensemblePredictions[[i]],main=paste0('c index = ',round(ensembleMetrics[[i]]$c.index,4),', rmse = ',round(ensembleMetrics[[i]]$rmse,4)))
	}
	dev.off(pdf.output)

	# Plot Kaplan-Meier plots of ensemble predictions #
	# pdf(paste0('Runs/',unid,'/',unid,'PlotEnsembleKaplanMeier.pdf'),width=12,height=8,onefile=TRUE)
	# pdf.output <- dev.cur()
	# for(i in 1:nReps){
	# 	PlotKaplanMeier(ensemblePredictions[[i]],trainingTestStructure[[i]]$testTargets,trainingTestStructure[[i]]$testEvents,'GPSurvCorrV')
	# }
	# dev.off(pdf.output)
}