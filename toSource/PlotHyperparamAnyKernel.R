PlotHyperparamAnyKernel <- function(outputStructureGP){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Plots hyperparamter values for each cycle of a GPSurv model with any covariance kernel 
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	modelType 		<- outputStructureGP$parameterStructure$modelType
	burnIn 			<- outputStructureGP$parameterStructure$burnIn

	noiseHypChosen 	<- outputStructureGP$logHypChosen$noise
	funcHypChosen 	<- outputStructureGP$logHypChosen$func
	lengthHypChosen <- outputStructureGP$logHypChosen$length

	colFunc 		<- colorRampPalette(c('lightcoral', 'darkred'))
	colLength 		<- colorRampPalette(c('greenyellow', 'darkgreen'))
	cols 			<- c('black',colFunc(length(lengthHypChosen)),colLength(length(lengthHypChosen)))
	ltys 			<- c(1,rep(2,length(funcHypChosen)),rep(5,length(lengthHypChosen)))
	main 			<- paste0('noise hyp = ',paste(round(exp(c(outputStructureGP$dataOptionsStructure$logHypGenerate$noise)),4),collapse=', '),'/n',
							', func hyp = ',paste(round(exp(c(outputStructureGP$dataOptionsStructure$logHypGenerate$func)),4),collapse=', '),'/n',
							', length hyp = ',paste(round(exp(c(outputStructureGP$dataOptionsStructure$logHypGenerate$length)),4),collapse=', '))

	if(modelType=='survival'){
		matplot(exp(outputStructureGP$logHypTable[,1:(1+length(funcHypChosen)+length(lengthHypChosen))]),col=cols,lty=ltys,type='l',ylab ='Hyperparameters',xlab='Cycle number',
				main=main)
		legend('topleft',legend=c(expression(sigma[n]^2),sapply(1:length(funcHypChosen), function(i) {as.expression(substitute(sigma[f[B]]^2,list(B = i)))}),paste0('l',1:length(lengthHypChosen))),
				col=cols,pch=c(NA,rep(NA,length(funcHypChosen)),rep(NA,length(lengthHypChosen))),lty=ltys)
	} else {
		if(burnIn){
			matplot(exp(rbind(unlist(t(outputStructureGP$parameterStructure$logHypStart[1:3])),
								unlist(t(outputStructureGP$parameterStructure$logHypBurnIn[1:3])),
								unlist(t(outputStructureGP$logHypChosen[1:3]))))
					/exp(unlist(outputStructureGP$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
		} else {
			matplot(exp(rbind(unlist(t(outputStructureGP$logHypChosen[1:3])),unlist(t(outputStructureGP$parameterStructure$logHypStart[1:3]))))
					/exp(unlist(outputStructureGP$dataOptionsStructure$logHypGenerate)[1:3]),
					type='l',ylab ='Hyperparameters (norm wrt generating values)',xlab='Cycle number')
		}
		legend('topleft',legend=c(expression(sigma[n]^2),expression(sigma[f]^2),'l'),col=cols,pch=c(NA,NA,NA),lty=ltys)
	}
}
