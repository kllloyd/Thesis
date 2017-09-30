AddFeature <- function(trainingTestStructure,feature,allParameterStructures,trainingTestStructureComplete,trainingIndices,testIndices){
	#---------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2017_08_09
	#---------------------------------------------------------------------------------------------------------------------------------#
	# Function adds features from the full trainingTestStructure to a subset trainingTestStructure
	# Also changes trainingTestStructure$dimension, dataOptionsStructure$dimension and parameterStructure$logHypStart
	#---------------------------------------------------------------------------------------------------------------------------------#

	featureIndex 						<- which(colnames(trainingTestStructureComplete$trainingData)%in%feature)
	trainingDataTemp 					<- trainingTestStructureComplete$trainingData[trainingIndices,featureIndex,drop=FALSE]
	colnames(trainingDataTemp) 			<- feature
	testDataTemp 						<- trainingTestStructureComplete$trainingData[testIndices,featureIndex,drop=FALSE]
	colnames(testDataTemp) 				<- feature

	if(colnames(trainingTestStructure$trainingData)[1]!='constant'){
		trainingTestStructure$trainingData 						<- cbind(trainingTestStructure$trainingData,trainingDataTemp)
		trainingTestStructure$testData 							<- cbind(trainingTestStructure$testData,testDataTemp)
		allParameterStructures$dataOptionsStructure$dimension 	<- allParameterStructures$dataOptionsStructure$dimension+1
		trainingTestStructure$dimension 						<- trainingTestStructure$dimension+1
	} else {
		trainingTestStructure$trainingData 						<- trainingDataTemp
		trainingTestStructure$testData 							<- testDataTemp
		allParameterStructures$dataOptionsStructure$dimension 	<- 1
		trainingTestStructure$dimension 						<- 1
	}
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