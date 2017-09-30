RandomGeneSelectionCuratedOvarian <- function(dataSet,nRand,notGenes){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Selects nRand random genes from all available in selected CuratedOvarian dataset
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	library(curatedOvarianData)
	data(list=paste0(dataSet,'_eset'))
	expressionData 				<- exprs(get(paste0(dataSet,'_eset')))
	remainingExpressionFeatures <- rownames(expressionData)[!(rownames(expressionData)%in%notGenes)]
	chosenExpressionFeatures 	<- remainingExpressionFeatures[sample(1:length(remainingExpressionFeatures),round(nRand+0.3*nRand,0))]

	renamedGeneFeatures 		<- RenameGeneFeaturesAccordingToDataSetCuratedOvarian(expressionData,chosenExpressionFeatures)
	allGenes 					<- renamedGeneFeatures$geneListPresentInDataSet

	randomGenes 				<- allGenes[sample(1:length(allGenes),nRand)]

	return(randomGenes)
}