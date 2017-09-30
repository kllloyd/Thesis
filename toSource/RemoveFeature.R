RemoveFeature <- function(trainingTestStructure,feature,allParameterStructures){
	#---------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2017_08_09
	#---------------------------------------------------------------------------------------------------------------------------------#
	# Function removes features from the trainingTestStructure to make a smaller subset
	# Also changes trainingTestStructure$dimension, dataOptionsStructure$dimension and parameterStructure$logHypStart
	#---------------------------------------------------------------------------------------------------------------------------------#

	featureIndex 						<- which(colnames(trainingTestStructure$trainingData)%in%feature)
	trainingTestStructure$trainingData 	<- trainingTestStructure$trainingData[,-featureIndex,drop=FALSE]
	trainingTestStructure$testData 		<- trainingTestStructure$testData[,-featureIndex,drop=FALSE]

	allParameterStructures$dataOptionsStructure$dimension 	<- allParameterStructures$dataOptionsStructure$dimension-1
	trainingTestStructure$dimension 						<- trainingTestStructure$dimension-1
	if(!allParameterStructures$parameterStructure$covFuncForm=='ARD'){
		allParameterStructures$parameterStructure$logHypStart$length 	<- log(0.9)
		allParameterStructures$parameterStructure$logHypStart$mean 		<- c(rep(0,allParameterStructures$dataOptionsStructure$dimension),0)
	} else {
		allParameterStructures$parameterStructure$logHypStart$length 	<- rep(log(0.9),allParameterStructures$dataOptionsStructure$dimension)
		allParameterStructures$parameterStructure$logHypStart$mean 		<- c(rep(0,allParameterStructures$dataOptionsStructure$dimension),0)
	}

	toReturn <- list('trainingTestStructure'=trainingTestStructure,'allParameterStructures'=allParameterStructures)
	return(toReturn)
}