CalculateEnsembleResults <- function(outputStructure,trainingTestStructure,bic,c.index){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Calculates ensemble predictions using results of subset bootstraps by GPSurvSqExpRSFS
	# Predictions weighted by rescaled exp(-BIC)
	# Also produces predictions using uniform weighting, ensemblePredictionsUnifWeight
	# marginalisedProb margainalises over models to give probability of feature being important
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	source('../toSource/CalculateMetrics.R')

	# Extract full feature list #
	features 				<- colnames(trainingTestStructure$trainingData)
	dimension 				<- length(features)

	# Extract feature list per bootstrap #
	nBootstraps 			<- length(outputStructure)
	featuresPerBootstrap 	<- list()
	for(j in 1:nBootstraps){
		featuresPerBootstrap[[j]] <- colnames(outputStructure[[j]]$trainingTestStructure$trainingData)
	}

	# Convert feature lists to 0/1 features for not-included/included #
	featuresBinary 				<- matrix(rep(0,nBootstraps*dimension),nrow=nBootstraps)
	colnames(featuresBinary) 	<- features
	for(j in 1:nBootstraps){
		for(k in 1:dimension){
			if(features[k]%in%featuresPerBootstrap[[j]]) featuresBinary[j,k] <- 1
		}
	}

	# Collate binary features, BIC, C Index #
	resultsAll 			<- as.data.frame(cbind(featuresBinary,bic,c.index))
	resultsAll$c.index 	<- 1-resultsAll$c.index 

	# Remove repeated subsets #
	# resultsUnique 		<- resultsAll[!duplicated(resultsAll[,1:10]),]
	resultsUnique 		<- resultsAll

	# Marginalise over models #
	marginalisedProb 	<- numeric(dimension)
	for(k in 1:dimension){
		forMarginalising 	<- exp(-resultsUnique$bic)
		forMarginalising[is.infinite(forMarginalising)&forMarginalising>0] <- 10^300
		forMarginalising[is.infinite(forMarginalising)&forMarginalising<0] <- -10^300
		marginalisedProb[k] <- forMarginalising%*%resultsUnique[,k]
	}

	# For reference #
	averageBIC 		<- numeric(dimension)
	averageCIndex 	<- numeric(dimension)
	for(k in 1:dimension){
		averageBIC[k] 		<- (resultsUnique$bic%*%resultsUnique[,k])/length(resultsUnique$bic)
		averageCIndex[k] 	<- (resultsUnique$c.index%*%resultsUnique[,k])/length(resultsUnique$c.index)
	}

	# Create ensemble predictions #
	marLik 		<- exp(-resultsUnique$bic)
	marLik[is.infinite(marLik)&marLik>0] <- 10^300
	marLik[is.infinite(marLik)&marLik<0] <- -10^300
	weightings 	<- marLik+min(marLik)
	if(sum(weightings)==0&all(weightings==0)) weightings <- rep(1/nBootstraps,length(weightings)) else weightings <- weightings/sum(weightings)

	predictions <- matrix(rep(NA,dim(trainingTestStructure$testData)[1]*nBootstraps),ncol=nBootstraps)
	for(i in 1:nBootstraps){
		predictions[,i] <- outputStructure[[i]]$funcMeanPred
	}

	ensemblePredictions <- numeric(dim(trainingTestStructure$testData)[1])
	for(m in 1:dim(trainingTestStructure$testData)[1]){
		ensemblePredictions[m] <- predictions[m,]%*%weightings
	}

	# Create ensemble predictions, uniformly weighted #
	ensemblePredictionsUnifWeight <- numeric(dim(trainingTestStructure$testData)[1])
	for(m in 1:dim(trainingTestStructure$testData)[1]){
		ensemblePredictionsUnifWeight[m] <- sum(predictions[m,]*(1/nBootstraps))
	}
	
	# Calculate ensemble metrics #
	ensembleMetrics 			<- CalculateMetrics(ensemblePredictions,trainingTestStructure$testTargets,trainingTestStructure$testEvents,NA,NA,NA)
	ensembleMetricsUnifWeight 	<- CalculateMetrics(ensemblePredictionsUnifWeight,trainingTestStructure$testTargets,trainingTestStructure$testEvents,NA,NA,NA)
	
	# Return #
	toReturn 		<- list('resultsAll'=resultsAll,'resultsUnique'=resultsUnique,'weightings'=weightings,
							'ensemblePredictions'=ensemblePredictions,'ensembleMetrics'=ensembleMetrics,
							'ensemblePredictionsUnifWeight'=ensemblePredictionsUnifWeight,'ensembleMetricsUnifWeight'=ensembleMetricsUnifWeight,
							'marginalisedProb'=marginalisedProb)
	return(toReturn)
}