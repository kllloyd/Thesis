ExtractSubsetOfFeatures <- function(trainingTestStructure,subsetDimension){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Extracts subset of full set of features, for use with StepBIC
	#--------------------------------------------------------------------------------------------------------------------------------------------#	
	featureIndices 										<- sample(1:(dim(trainingTestStructure$trainingData)[2]),subsetDimension,replace=FALSE)
	trainingTestStructure$trainingData 					<- trainingTestStructure$trainingData[,featureIndices,drop=FALSE]
	trainingTestStructure$testData 						<- trainingTestStructure$testData[,featureIndices,drop=FALSE]
	trainingTestStructure$dimension 					<- subsetDimension

	return(trainingTestStructure)
}