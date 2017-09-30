ModifyFeatureListsInformedARD2 <- function(dataOptionsStructure,parameterStructure,trainingTestStructure){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# If covFuncForm is InformedARD and IARD2 is wanted, need to modify gene lists and make new length hyperparameters for the new lists
	# New lists will be contructed from intesections of old lists -> make new hyperparameter mean of old hyperparameters
	# For each element, look at which lists it is in
	# Each combination of lists becomes a new possible list: sum_{n=1}^{N} N!/(n!(N-n)!)
	# If no elements in a list, throw it away
	# New hyperparameter is the mean of the hyperparameters of the contributing lists
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	lengthHyp 		<- parameterStructure$logHypStart$length
	funcHyp 		<- parameterStructure$logHypStart$func
	extraParam 		<- parameterStructure$extraParam
	covFuncForm		<- parameterStructure$covFuncForm
	dimension 		<- dataOptionsStructure$dimension
	featureNames 	<- colnames(trainingTestStructure$trainingData)

	# Make all possible new lists from old lists
	numNewLists 		<- sum(sapply(1:length(extraParam),function(x) factorial(length(extraParam))/(factorial(x)*factorial(length(extraParam)-x))))
	listNames 			<- t(expand.grid(lapply(numeric(length(extraParam)), function(x) c(TRUE,FALSE)))[-(numNewLists+1),])
	colnames(listNames) <- paste0('newlist',1:numNewLists)
	rownames(listNames) <- paste0('oldlist',1:length(extraParam))

	# Calculate new length hyperparameter for all possible new lists
	newLengthHyp 		<- sapply(1:numNewLists, function(x) log(mean(exp(lengthHyp[listNames[,x]]))))
	if(covFuncForm=='InformedARDV3') newFuncHyp <- sapply(1:numNewLists, function(x) log(mean(exp(funcHyp[listNames[,x]]))))

	allFeatures 				<- matrix(NA,dimension,length(extraParam))
	allFeaturesNewListNo 		<- numeric(dimension)
	for(i in 1:dimension){
	  	for(j in 1:length(extraParam)){
	  		# for each feature make a vector of list names for each feature
	  		allFeatures[i,j] 	<- ifelse(featureNames[i]%in%extraParam[[j]],paste0('oldlist',j),NA)
	  	}
	  	# for each feature assign a new list number based on combination of list names
		allFeaturesNewListNo[i] <- colnames(listNames)[sapply(1:dim(listNames)[2],function(x) setequal(allFeatures[i,][!is.na(allFeatures[i,])],rownames(listNames)[listNames[,x]]))]
	}

	newExtraParam <- list()
	for(i in 1:numNewLists){
		if(any(allFeaturesNewListNo%in%colnames(listNames)[i])){
			# For each new list, find all features on it
			newExtraParam[[i]] 	<- featureNames[allFeaturesNewListNo%in%colnames(listNames)[i]]
		}
	}
	newLengthHyp 				<- newLengthHyp[lapply(newExtraParam,length)>0]
	if(covFuncForm=='InformedARDV3') newFuncHyp <- newFuncHyp[lapply(newExtraParam,length)>0]
	newExtraParam 				<- newExtraParam[lapply(newExtraParam,length)>0]
	names(newExtraParam) 		<- paste0('list',1:length(newExtraParam))

	# to return
	parameterStructure$logHypStart$length 	<- newLengthHyp
	if(covFuncForm=='InformedARDV3') parameterStructure$logHypStart$func 	<- newFuncHyp
	parameterStructure$extraParam 			<- newExtraParam

	return(parameterStructure)
}