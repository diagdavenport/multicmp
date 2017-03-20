#################################################################################
## Program name: suma.R                                                        ##
## Authors: Kimberly Sellers and Darcy Steeg Morris                            ##
## Date: 8/20/2015                                                             ##
## Decription: Calculates terms of the inner sum (over a) of the BCMP density. ##
## Inputs: BCMP parameters (bivprob) and outer sum summation index (n).        ##
#################################################################################

suma <- function(bivprob, x, y, n){
	p <- bivprob/sum(bivprob)
	result <- rep(0,length(x))
	for (i in 1:length(x)) {
		a <- (n-x[i]-y[i]):n
		resulti <- rep(0, length(a))
		for (place in 1:length(a)) {
			terms <- c(a[place], n-a[place]-y[i], n-a[place]-x[i], x[i]+y[i]+a[place]-n)
			if(any(terms < 0)) resulti[place] <- 0
			else resulti[place] <- dmultinom(terms, prob=c(p[1],p[3],p[2],p[4]))
		}
		result[i] = sum(resulti) 
	}
	return(result)
}

