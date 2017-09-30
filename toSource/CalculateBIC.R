CalculateBIC <- function(outputStructure,toCount){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Calcuates BIC, either counting features or hyperparamters
	# Uses negative log marginal likelihood output by GP or GPSurv models, passed via outputStructure
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	objectiveTable 	<- outputStructure$objectiveTable
	negLogMarLik 	<- objectiveTable[length(objectiveTable)]
	n 				<- dim(outputStructure$trainingTestStructure$trainingData)[1]
	if(toCount=='features'){
		p 			<- dim(outputStructure$trainingTestStructure$trainingData)[2]
	} else if(toCount=='hyperparameters'){
		p 			<- length(unlist(outputStructure$logHypChosen))
	}

	BIC 			<- 2*negLogMarLik + p*log(n)

	return(BIC)
}