ApplyGPSurvBIC <- function(trainingTestStructureForNReps,allParameterStructures){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Applying BIC stepwise feature selection with GP for survival model
	# Splits the training set into new training set and test set for feature selection.
	# When feature subset chosen, model fitted on full training set and used to predict on unseen test set
	# Two stepwise modes available:
	# 	* forward - Forward selection
	# 	* backward - Backwards elimination
	# Mode indicated as variable allParameterStructures$parameterStructure$direction
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	#-------------------------------------------------------------------------------------------------------#
	#-------------------------------- Load required libraries and functions --------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	source('../toSource/CalculateBIC.R')
	source('../toSource/RunModelCombinationsBackward.R')
	source('../toSource/RunModelCombinationsForward.R')
	source('../toSource/RemoveFeature.R')
	source('../toSource/AddFeature.R')

	#-------------------------------------------------------------------------------------------------------#
	#----------------------------------- Extract variables and initalise -----------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	timeStart 											<- Sys.time()
	allParameterStructures$plotSaveOptions$savePlots 	<- FALSE
	allParameterStructures$plotSaveOptions$saveFiles 	<- FALSE
	count 												<- 1
	direction 											<- allParameterStructures$parameterStructure$direction
	# Set up new training and test sets as subsets of training set #
	trainingTestStructure 								<- trainingTestStructureForNReps
	trainingIndices 									<- sample(1:(dim(trainingTestStructureForNReps$trainingData)[1]),round(0.9*dim(trainingTestStructureForNReps$trainingData)[1]),replace=FALSE)
	testIndices 										<- sample((1:(dim(trainingTestStructureForNReps$trainingData)[1]))[-trainingIndices])
	trainingTestStructure$trainingData 					<- trainingTestStructureForNReps$trainingData[trainingIndices,,drop=FALSE]
	trainingTestStructure$trainingTargets 				<- trainingTestStructureForNReps$trainingTargets[trainingIndices,,drop=FALSE]
	trainingTestStructure$events 						<- trainingTestStructureForNReps$events[trainingIndices,,drop=FALSE]
	trainingTestStructure$testData 						<- trainingTestStructureForNReps$trainingData[testIndices,,drop=FALSE]
	trainingTestStructure$testTargets 					<- trainingTestStructureForNReps$trainingTargets[testIndices,,drop=FALSE]
	trainingTestStructure$testEvents 					<- trainingTestStructureForNReps$events[testIndices,,drop=FALSE]
	trainingTestStructure$trainingTargetsPreCensoring 	<- trainingTestStructureForNReps$trainingTargetsPreCensoring[trainingIndices,,drop=FALSE]
	trainingTestStructure$testTargetsPreCensoring 		<- trainingTestStructureForNReps$trainingTargetsPreCensoring[testIndices,,drop=FALSE]
	trainingTestStructure$nTraining 					<- length(trainingIndices)
	trainingTestStructure$nTest 						<- length(testIndices)
	logHypStart 										<- allParameterStructures$parameterStructure$logHypStart

	#-------------------------------------------------------------------------------------------------------#
	#------------------------------------- Implement feature selection -------------------------------------#
	#-------------------------------------------------------------------------------------------------------#
	# Implement stepwise BIC, starting with full or empty model #
	if(direction=='forward'){
		trainingTestStructure$trainingData 				<- matrix(rep(1,dim(trainingTestStructure$trainingData)[1])+rnorm(dim(trainingTestStructure$trainingData)[1],0,0.001),ncol=1)
		colnames(trainingTestStructure$trainingData) 	<- 'constant'
		trainingTestStructure$testData 					<- matrix(rep(1,dim(trainingTestStructure$testData)[1])+rnorm(dim(trainingTestStructure$testData)[1],0,0.001),ncol=1)
		colnames(trainingTestStructure$testData) 		<- 'constant'
		trainingTestStructure$dimension 				<- 1
		dataOptionsStructureForward 					<- allParameterStructures$dataOptionsStructure
		dataOptionsStructureForward$dimension 			<- 1
		# allParameterStructures$parameterStructure$logHypStart <- list('noise'=logHypStart$noise,'func'=logHypStart$func,'length'=logHypStart$length,'mean'=logHypStart$mean[c(1,length(logHypStart$mean))])
		outputStructureForward 							<- ApplyGP(trainingTestStructure,dataOptionsStructureForward,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		bic 											<- CalculateBIC(outputStructureForward,'features')
		cat('Empty model BIC =',bic,', features =',colnames(trainingTestStructure$trainingData),fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	} else if(direction=='backward'){
		outputStructureBackward 						<- ApplyGP(trainingTestStructure,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
		bic 											<- CalculateBIC(outputStructureBackward,'features')
		cat('Full model BIC =',bic,', features =',colnames(trainingTestStructure$trainingData),fill=TRUE)
		cat('-------------------------------------------------------',fill=TRUE)
	}
	continueFlag 										<- TRUE
	while(continueFlag&count<=trainingTestStructureForNReps$dimension){
		count 												<- count+1
		allParameterStructures$plotSaveOptions$savePlots 	<- FALSE
		allParameterStructures$plotSaveOptions$saveFiles 	<- FALSE
		if(direction=='forward'){
			# Calculate BIC when adding each feature #
			outputStructure 				<- RunModelCombinationsForward(trainingTestStructure,allParameterStructures,trainingTestStructureForNReps,trainingIndices,testIndices)
			bicsThisCycle 					<- outputStructure$bic
			bicMin 							<- min(bicsThisCycle)
			bic 							<- rbind(bic,bicMin)
			continueFlag 					<- (bic[count]<=bic[count-1])|(count==1)
			# Add most useful feature #
			if(continueFlag){
				featureToAdd 				<- names(bicsThisCycle)[which.min(bicsThisCycle)]
				# browser()
				increasedFeatureStructure 	<- AddFeature(trainingTestStructure,featureToAdd,allParameterStructures,trainingTestStructureForNReps,trainingIndices,testIndices)
				allParameterStructures 		<- increasedFeatureStructure$allParameterStructures
				trainingTestStructure 		<- increasedFeatureStructure$trainingTestStructure
			}
		} else if(direction=='backward'){
			# Calculate BIC when removing each feature #
			outputStructure 				<- RunModelCombinationsBackward(trainingTestStructure,allParameterStructures)
			bicsThisCycle 					<- outputStructure$bic
			bicMax 							<- max(bicsThisCycle)
			bic 							<- rbind(bic,bicMax)
			continueFlag 					<- (bic[count]>=bic[count-1])
			# Remove least useful feature #
			if(continueFlag){
				featureToRemove 			<- names(bicsThisCycle)[which.min(bicsThisCycle)]
				reducedFeatureStructure 	<- RemoveFeature(trainingTestStructure,featureToRemove,allParameterStructures)
				allParameterStructures 		<- reducedFeatureStructure$allParameterStructures
				trainingTestStructure 		<- reducedFeatureStructure$trainingTestStructure
			}
		}
	}

	#-------------------------------------------------------------------------------------------------------#
	#--------------- Fit final model, using full training set, and predict on unseen test set --------------#
	#-------------------------------------------------------------------------------------------------------#
	allParameterStructures$plotSaveOptions$savePlots 	<- TRUE
	allParameterStructures$plotSaveOptions$saveFiles 	<- FALSE
	trainingTestStructureFinal 							<- trainingTestStructureForNReps
	trainingTestStructureFinal$trainingData 			<- trainingTestStructureFinal$trainingData[,colnames(trainingTestStructureFinal$trainingData)%in%colnames(trainingTestStructure$trainingData),drop=FALSE]
	trainingTestStructureFinal$testData 				<- trainingTestStructureFinal$testData[,colnames(trainingTestStructureFinal$testData)%in%colnames(trainingTestStructure$testData),drop=FALSE]
	trainingTestStructureFinal$dimension 				<- trainingTestStructure$dimension
	outputStructureGPSurvBIC 							<- ApplyGP(trainingTestStructureFinal,allParameterStructures$dataOptionsStructure,allParameterStructures$parameterStructure,allParameterStructures$plotSaveOptions)
	cat('Chosen features =',colnames(trainingTestStructureFinal$trainingData),fill=TRUE)
	cat('-------------------------------------------------------',fill=TRUE)
	
	timeEnd 											<- Sys.time()
	outputStructureGPSurvBIC$timeTaken 					<- difftime(timeEnd,timeStart, units='min')
	outputStructureGPSurvBIC$bic 						<- bic
	outputStructureGPSurvBIC$finalCycleBIC 				<- bicsThisCycle

	return(outputStructureGPSurvBIC)
}