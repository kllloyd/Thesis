RunModelCombinationsBackward <- function(trainingTestStructureAll,allParameterStructuresFull){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Applies GPSurvBICBackward, removing and replacing each feature in turn and calculating BICs
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	modelFeatures 				<- colnames(trainingTestStructureAll$trainingData)
	nFeatures					<- length(modelFeatures)

	outputStructureGPSurv 	<- list()
	bic 					<- rep(NA,nFeatures)
	names(bic) 				<- modelFeatures
	for(i in 1:nFeatures){
		reducedFeatureStructure 		<- RemoveFeature(trainingTestStructureAll,modelFeatures[i],allParameterStructuresFull)
		allParameterStructures 			<- reducedFeatureStructure$allParameterStructures
		trainingTestStructure 			<- reducedFeatureStructure$trainingTestStructure
		outputStructureGPSurv[[i]] 		<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		bic[i] 							<- CalculateBIC(outputStructureGPSurv[[i]],'features')
		cat((dim(trainingTestStructure$trainingData)[2]),'feature model BIC =',bic[i],', features =',colnames(trainingTestStructure$trainingData),fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	}

	toReturn <- list('outputStructure'=outputStructureGPSurv,'bic'=bic)
	return(toReturn)
}