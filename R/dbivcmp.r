######################################################################
## Program name: dbivCMP.R                                          ##
## Authors: Kimberly Sellers and Darcy Steeg Morris                 ##
## Date: 8/20/2015                                                  ##
## Decription: Calculates BCMP density.                             ##
## Inputs: BCMP parameters (lambda, nu & bivprob) and data (x & y). ##
## Output: P(X=x,Y=y)                                               ##
######################################################################

dbivCMP <- function(lambda, nu, bivprob, x, y, maxit) {

	# lambda (constant) -- mean under Poisson model
	# nu (constant) -- dispersion parameter
	# bivprob -- c(p00, p01, p10, p11)
	# maxit -- stopping point making infinite sum finite

	# Compute the Z function #
	forans <- rep(0,maxit+1)
	for (j in 1:maxit) {
		temp <- rep(0,j)
		for (i in 1:j) {
			temp[i] <- lambda/(i^nu)
		}
		forans[j+1] <- prod(temp)
	}
	forans[1] <- 1
	z = sum(forans)

	# Calculate the (in)finite sum #
	forans <- matrix(0,nrow=length(x),ncol=maxit+1)
	for (n in 1:maxit){
		temp <- rep(0,n)
		for (i in 1:n){
			temp[i] <- lambda/(i^nu)
		}
		forans[,n+1] <- prod(temp)*suma(bivprob, x, y, n)  
	}
	forans[,1] <- suma(bivprob, x, y, 0)
	ans <- (1/z) * apply(forans,1,sum)
	return(ans)
}
