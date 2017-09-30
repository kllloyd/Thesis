LogPriorX <- function(hyp,x,data,targetsUncensored,targetsCensored,target,calcDeriv){

	switch(hyp,
		'noise'	={k <- 2
				  t <- 1},
		'func'	={k <- 2
				  t <- 1},
		'length'={k <- 2
				  t <- diff(quantile(unlist(data),probs=c(0.05,0.95)))},
		'targets'	={if(mean(exp(targetsUncensored))>exp(target)){
						x <- exp(x)-exp(target)
						t <- (-(mean(exp(targetsUncensored))-exp(target))+sqrt((mean(exp(targetsUncensored))-exp(target))^2+4*4*sd(exp(c(targetsUncensored)))^2))/2
						k <- ((mean(exp(targetsUncensored))-exp(target))/t)+1
					  } else {
						x <- exp(x)-exp(target)
						k <- 0.99
						t <- sqrt(4*sd(exp(c(targetsUncensored)))^2/k)
					  }}
	)

	if(any(x<0)){
		logPriorX <- -10^300
	} else {
		logPriorX <- -log(gamma(k)) - k*log(t) + (k-1)*log(x) -x/t
	}

	logPriorX <- ifelse(is.infinite(logPriorX)&logPriorX<0,-10^300,logPriorX) 		# Prevent returning Inf or -Inf
	logPriorX <- ifelse(is.infinite(logPriorX)&logPriorX>0,10^300,logPriorX)

	if(calcDeriv){
		dLogPriorXdX <- (k-1)*(x^-1)-t
	}

	if(!calcDeriv){
		return(logPriorX)
	} else {
		return(dLogPriorXdX)
	}
}
