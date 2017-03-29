####################################################################
## Program name: ComputeConstantBCMPests.R                        ##
## Authors: Kimberly Sellers and Darcy Steeg Morris               ##
## Date: 8/20/2015                                                ##
## Decription: Estimates BCMP parameters via ML.                  ##
## Inputs: Data (x & y).                                          ##
## Output: MLEs (lambda, nu, p00, p10, p01 and p11), negative LL. ##
####################################################################

## Added 3/22/17: LRT (nu=1) and bivariate Poisson model ##

ComputeConstantBCMPests <- function(data,max,startvalues=NULL) {
	if (dim(data)[2] != 2) {
		stop('data must have 2 columns')
	}
	x <-  data[,1]
	y <- data[,2]
	maxit <- max
	# Likelihood Function #
	minusloglike <- function(par){
		lambda <- par[1]
		nu <- par[2]
		p <- par[3:6]
		nll <- -1*sum(log(dbivCMP(lambda=lambda,nu=nu,bivprob=p,x=x,y=y,maxit=maxit)))
		return(nll)
	}
	
	invisible(capture.output(bp <- simple.bp(x,y)))
	
	# Starting Values #
	if (is.null(startvalues)) {
		p10 = .25
		p01 = .25
		p00 = .25
		p11 = .25
		lambda_start = 1
		nu_start = 1
	}
	if (is.null(startvalues)==0) {
		p10 = startvalues[5]
		p01 = startvalues[4]
		p00 = startvalues[3]
		p11 = startvalues[6]
		lambda_start = startvalues[1]
		nu_start = startvalues[2]
	}
	
	# NLMINB to find MLEs #
	cat("Iterating...", "\n")
	BCMPests <- nlminb(start=c(lambda_start,nu_start,p00,p01,p10,p11), minusloglike, lower = c(rep(0,6)), upper = c(Inf,Inf,Inf,Inf,Inf,Inf),control=list(trace=10,iter.max=1000))
	LRT_bpd <- -2*(bp$loglikelihood[length(bp$loglikelihood)] - (-1*BCMPests$obj))
	p_bpd <- 1 - pchisq(abs(LRT_bpd),df=1)
	
	# prepare tabular data for printing
	cat("\n", "The parameter estimates are as follows:","\n")
	par.est <- data.frame("Parameter" = c("lambda", "nu", "p00", "p10", "p01", "p11") , "MLE" = c(BCMPests$par[1:2],BCMPests$par[3:6]/sum(BCMPests$par[3:6])) , "SE" = rep("?", 6))
	print(par.est)

	cat("\n" , "Hypothesis test results:","\n")
	hyp.tst <- data.frame("Likelihood ratio test" = LRT_bpd, "p-value" = p_bpd)
	print(hyp.tst)
	
	cat("\n")
	
	# finally return results
	return(list(par=c(BCMPests$par[1:2],BCMPests$par[3:6]/sum(BCMPests$par[3:6])),negll=BCMPests$obj,LRT_bpd=LRT_bpd,p_bpd=p_bpd))
	
}

