SubsetIndices <- function(dimension,subsetDimension,nBootstraps){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Generating subsets containing a given number of features from the full training and test sets
	# Generates nBootstraps from total number choose(dimension,subsetDimension)
	# Two methods:
	# 	* Small numbers of possible subsets - generate all and select randomly
	# 	* Large number of possible subsets - select randomly from uniform distribution
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	numCombinations 		<- choose(dimension,subsetDimension)
	totalNumCombinations 	<- sum(numCombinations)
	if(totalNumCombinations<10^10){
		subsetIndices 		<- sample(totalNumCombinations,nBootstraps,replace=FALSE)	
	} else {
		subsetIndices 		<- ceiling(runif(nBootstraps,min=1,max=totalNumCombinations))
	}

	return(subsetIndices)
}
