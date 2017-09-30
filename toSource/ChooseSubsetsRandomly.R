ChooseSubsetsRandomly <- function(trainingTestStructure,subsetDimension,subsetIndex){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Generating subsets containing a given number of features from the full training and test sets
	# Generates nBootstraps from total number choose(dimension,subsetDimension), as specified by subsetIndex, using combinadic
	# If subsetDimension is not a single value, need to re-index
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	source('../toSource/combinadic.R')

	allFeatures 			<- colnames(trainingTestStructure$trainingData)
	dimension 				<- length(allFeatures)
	numCombinations 		<- choose(dimension,subsetDimension)
	if(subsetIndex<=numCombinations[1]){
		featureIndices 		<- combinadic(dimension,subsetDimension[1],subsetIndex)
		dimension 			<- subsetDimension[1]
	} else {
		for(i in 2:length(subsetDimension)){
			if(sum(numCombinations[1:(i-1)])<subsetIndex&subsetIndex<(sum(numCombinations[1:i])+1)){
				subsetIndexShifted 		<- subsetIndex-sum(numCombinations[1:(i-1)])
				featureIndices 			<- combinadic(dimension,subsetDimension[i],subsetIndexShifted)
				dimension 				<- subsetDimension[i]
			}
		}
	}

	trainingTestStructure$trainingData 	<- trainingTestStructure$trainingData[,featureIndices,drop=FALSE]
	trainingTestStructure$testData 		<- trainingTestStructure$testData[,featureIndices,drop=FALSE]
	trainingTestStructure$dimension 	<- dimension

	return(trainingTestStructure)
}
