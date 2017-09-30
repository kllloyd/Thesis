#--------------------------------------------------------------------------------------------------------------------------------------------#
# K Lloyd 2017_08_09
#--------------------------------------------------------------------------------------------------------------------------------------------#
# Diagram of random functions drawn from prior and posterior of GP
# Recreated from Rassmussen and Williams (2006) Gaussian processes for machine learning, Figure 2.2
# Modified from http://katbailey.github.io/post/gaussian-processes-for-dummies/
#--------------------------------------------------------------------------------------------------------------------------------------------#

source('../toSource/CovFunc.R')
source('../toSource/add.alpha.R')

# Hyperparameters and parameters #
hypNoise 		<- 0.01
hypFunc 		<- 0.99
hypLength 		<- 0.7
trainingData 	<- NA
trainingTargets <- NA
extraParam 		<- NA
nDraws 			<- 5 	# Number of functions drawn

# Generate training and test data #
nTrain			<- 7
trainingData 	<- matrix(runif(n=nTrain,min=-5,max=5),ncol=1)
trainingTargets <- sin(trainingData)
nTest 			<- 100
testData 		<- seq(from=-5,to=5,length.out=nTest) 	# Used to make curves

# Covariance matrices #
K 				<- CovFunc(trainingData,trainingData,trainingData,trainingTargets,extraParam,hypFunc,hypLength,'SqExp') + diag(rep(1e-10,nTrain))
L 				<- t(chol(K))
KStar 			<- CovFunc(trainingData,testData,trainingData,trainingTargets,extraParam,hypFunc,hypLength,'SqExp')
KStarStar 		<- CovFunc(testData,testData,trainingData,trainingTargets,extraParam,hypFunc,hypLength,'SqExp') + diag(rep(1e-10,nTest))
LStarStar 		<- t(chol(KStarStar))

# Calculate mean and sd of predictions #
v 				<- solve(L,KStar)
alpha 			<- solve(t(L),(solve(L,(trainingTargets))))
mu 				<- t(KStar)%*%alpha
s2 				<- diag(KStarStar-t(v)%*%v)
stdv 			<- sqrt(s2)

# Draw from prior: ~ K** N(0,1)
fPrior 			<- list()
for(i in 1:nDraws){
	fPrior[[i]] <- LStarStar%*%matrix(rnorm(nTest),ncol=1)
}

# Plot prior draws #
cols 			<- colorRampPalette(c('green', 'purple'))(nDraws)
pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotDrawsFromPrior.pdf'),width=6,height=6)
pdf.output <- dev.cur()
	plot(1,1,xlim=c(min(testData),max(testData)),ylim=range(c(unlist(fPosterior),mu-2*stdv,mu+2*stdv),unlist(fPrior)),col='white',xlab='Data, x',ylab='Outcome, f(x)',xaxs="i")
	polygon(c(testData,rev(testData)),c(-rep(2*sqrt(hypFunc),nTest),rep(2*sqrt(hypFunc),nTest)),col=add.alpha('skyblue',0.5),border=FALSE)
	for(i in 1:nDraws){
		lines(testData,fPrior[[i]],col=cols[i])
	}
dev.off(pdf.output)

# Draw from posterior: ~ mean + (K**-vv) N(0,1)
LPost  			<- t(chol(KStarStar-t(v)%*%v+1e-10*diag(nTest)))
fPosterior 		<- list()
for(i in 1:nDraws){
	fPosterior[[i]] <- mu+LPost%*%matrix(rnorm(nTest),ncol=1)
}

# Plot posterior draws #
pdf(file=paste0(getwd(),"/",'Runs',"/",'PlotDrawsFromPosterior.pdf'),width=6,height=6)
pdf.output <- dev.cur()
	plot(1,1,xlim=c(min(testData),max(testData)),ylim=range(c(unlist(fPosterior),mu-2*stdv,mu+2*stdv),unlist(fPrior)),col='white',xlab='Data, x',ylab='Outcome, f(x)',xaxs="i")
	polygon(c(testData,rev(testData)),c(mu-2*stdv,rev(mu+2*stdv)),col=add.alpha('skyblue',0.5),border=FALSE)
	for(i in 1:nDraws){
		lines(testData,fPosterior[[i]],col=cols[i])
	}
	points(trainingData,trainingTargets,pch=20)
	lines(testData,mu,lty=2,col='black')
dev.off(pdf.output)

