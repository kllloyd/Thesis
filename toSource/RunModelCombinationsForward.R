RunModelCombinationsForward <- function(trainingTestStructureSmall,allParameterStructuresOld,trainingTestStructureComplete,trainingIndices,testIndices){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Applies GPSurvBICForward, adding and removing each feature in turn and calculating BICs
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	modelFeatures 			<- colnames(trainingTestStructureComplete$trainingData)
	existingFeatures 		<- colnames(trainingTestStructureSmall$trainingData)
	featuresToRun 			<- modelFeatures[!(modelFeatures%in%existingFeatures)]
	nFeatures 				<- length(featuresToRun)

	outputStructureGPSurv 	<- list()
	bic 					<- rep(NA,nFeatures)
	names(bic) 				<- featuresToRun
	for(i in 1:nFeatures){
		increasedFeatureStructure 		<- AddFeature(trainingTestStructureSmall,featuresToRun[i],allParameterStructuresOld,trainingTestStructureComplete,trainingIndices,testIndices)
		allParameterStructures 			<- increasedFeatureStructure$allParameterStructures
		trainingTestStructure 			<- increasedFeatureStructure$trainingTestStructure
		outputStructureGPSurv[[i]] 		<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		bic[i] 							<- CalculateBIC(outputStructureGPSurv[[i]],'features')
		cat((dim(trainingTestStructure$trainingData)[2]),'feature model BIC =',bic[i],', features =',colnames(trainingTestStructure$trainingData),fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	}
	
	toReturn <- list('outputStructure'=outputStructureGPSurv,'bic'=bic)
	return(toReturn)
}