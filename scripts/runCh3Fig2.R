#--------------------------------------------------------------------------------------------------------------------------------------------#
# K Lloyd 2017_08_09
#--------------------------------------------------------------------------------------------------------------------------------------------#
# Diagrams of hyperpriors for noise, function and length hyperparameters
#--------------------------------------------------------------------------------------------------------------------------------------------#

unid 				<- outputStructureGPSurvNoCorr[[1]]$parameterStructure$unid
oldMar 				<- par('mar')

# Noise, func and length #
x 					<- seq(0.01,20,length.out=200)
trainingData 		<- outputStructureGPSurvNoCorr[[1]]$trainingTestStructure$trainingData
ks 					<- c(2,2)
ts 					<- c(1,diff(quantile(unlist(trainingData),probs=c(0.05,0.95))))
xlab 				<- c(expression(sigma[n]^2 ~ or ~ sigma[f]^2 ~ hyperparameter ~ value),expression(Length ~ hyperparameter ~ value))

pdf(paste0('Runs/',unid,'/','PlotHyperpriors1.pdf'),width=8,height=6,onefile=TRUE)
pdf.output <- dev.cur()
	for(i in c(1:2)){
		par(mar=c(5.1,4.4,4.1,2.1))
		plot(x,exp((function(k,t,x) -log(gamma(k)) - k*log(t) + (k-1)*log(x) -x/t)(k=ks[i],t=ts[i],x=x)),type='l',xlab=xlab[i],ylab='Hyperprior',las=1,lwd=1.3)
		par(mar=oldMar)
	}
dev.off(pdf.output)
