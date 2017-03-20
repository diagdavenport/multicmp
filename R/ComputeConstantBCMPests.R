####################################################################
## Program name: ComputeConstantBCMPests.R                        ##
## Authors: Kimberly Sellers and Darcy Steeg Morris               ##
## Date: 8/20/2015                                                ##
## Decription: Estimates BCMP parameters via ML.                  ##
## Inputs: Data (x & y).                                          ##
## Output: MLEs (lambda, nu, p00, p10, p01 and p11), negative LL. ##
####################################################################

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
	BCMPests <- nlminb(start=c(lambda_start,nu_start,p00,p01,p10,p11), minusloglike, lower = c(rep(0,6)), upper = c(Inf,Inf,Inf,Inf,Inf,Inf),control=list(trace=1,iter.max=1000))
	return(list(par=c(BCMPests$par[1:2],BCMPests$par[3:6]/sum(BCMPests$par[3:6])),negll=BCMPests$obj))
}
