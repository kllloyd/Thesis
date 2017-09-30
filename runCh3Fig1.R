#----------------------------------------------------------------------------------------#
# K Lloyd 2017_08_09
#----------------------------------------------------------------------------------------#
# Diagram of adjusted mean and variance by truncation of prediction normal distribution
# Equal area under the curve
#----------------------------------------------------------------------------------------#
source('../toSource/AdjustTrainingSurvivalMeanVariance.R')
library(zoo)

# Generate untruncated curve #
data 					<- seq(-10,10,length.out=5000)
idD 					<- order(data)
origMean 				<- 0
origSD					<- 2
densityCurve 			<- dnorm(data,mean=origMean,sd=origSD)

# Truncated curve #
truncPoint 				<- 1.5
truncatedData 			<- data[which(data>=truncPoint)]
truncatedDensityCurve 	<- densityCurve[which(data>=truncPoint)]
idTD 					<- order(truncatedData)
truncatedAUC 			<- sum(diff(truncatedData[idTD])*rollmean(truncatedDensityCurve[idTD],2))

# Normal approximating truncated curve #
adjStructure 			<- AdjustTrainingSurvivalMeanVariance(truncPoint,origMean,origSD^2)
adjMean  				<- adjStructure$mu
adjSD 					<- sqrt(adjStructure$sigmaSquared)
approxDensityCurve 		<- dnorm(data,mean=adjMean,sd=adjSD)
approxAUC 				<- sum(diff(data[idD])*rollmean(approxDensityCurve[idD],2))

normalisedDensityCurve 	<- (approxDensityCurve/approxAUC)*truncatedAUC

# Plot #
pdf(paste0('Runs/','AdjustingMeanVariance.pdf'),width=10,height=8,onefile=TRUE)
pdf.output <- dev.cur()
	cols <- c('black','skyblue','black','purple')
	ltys <- c(4,NA,1,1)
	plot(data,densityCurve,type='l',col=cols[1],lty=ltys[1],xlab='Log time',ylab='Probability density')
	segments(x0=origMean,x1=origMean,y0=0,y1=max(densityCurve),col=cols[1],lty=ltys[1])
	polygonX <- c(truncPoint,truncatedData,truncatedData[length(truncatedData)])
	polygonY <- c(0,truncatedDensityCurve,truncatedDensityCurve[length(truncatedDensityCurve)])
	polygon(polygonX,polygonY,col=cols[2],border=NA)
	segments(x0=truncPoint,x1=truncPoint,y0=0,y1=densityCurve[data>0][which.min(abs(data[data>0]-truncPoint))],col=cols[3],lty=ltys[3])
	lines(truncatedData,truncatedDensityCurve,col=cols[3],lty=ltys[3])
	lines(data,normalisedDensityCurve,col=cols[4],lty=ltys[4])
	segments(x0=adjMean,x1=adjMean,y0=0,y1=max(normalisedDensityCurve),col=cols[4],lty=ltys[4])
	legend('topright',c('pre-truncation distribution','truncated distribution','approximate distribution'),col=cols[c(1,2,4)],lty=ltys[c(1,2,4)],pch=c(NA,15,NA))
dev.off(pdf.output)
