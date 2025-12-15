#Simulations for entropy-complexity test of serial independence

#All possible SOPs:
getPerms <- function(x) {
    if (length(x) == 1) {
        return(x)
    }
    else {
        res <- matrix(nrow = 0, ncol = length(x))
        for (i in seq_along(x)) {
            res <- rbind(res, cbind(x[i], Recall(x[-i])))
        }
        return(res)
    }
}


#Entropy:
H <- function(p){
	value <- 0
	p0 <- p[p>0] #to avoid log(0)
	value <- value - sum(p0*log(p0))
	value
}


#Normalized entropy:
H.norm <- function(p){
	H(p) / log(length(p))
}


#JS divergence:
D <- function(p){
	mfak <- length(p)
	u <- rep(1/mfak, mfak)
	
	H((p+u)/2) - H(p)/2 - log(mfak)/2
}


#Maximal JS divergence:
D.max <- function(mfak){
	(-(mfak+1)/2/mfak * log(mfak+1) + log(2*mfak) - log(mfak)/2)
}


#Normalized JS divergence:
D.norm <- function(p){
	mfak <- length(p)
	
	D(p) / D.max(mfak)
}


#(Normalized) complexity:
C <- function(p) H.norm(p) * D.norm(p)








tabd <- 1 #only delay 1

m <- 3
mfak <- factorial(m)

# getPerms(1:m)
     # [,1] [,2] [,3]
# [1,]    1    2    3
# [2,]    1    3    2
# [3,]    2    1    3
# [4,]    2    3    1
# [5,]    3    1    2
# [6,]    3    2    1



#Lexicographic arrangement as in manuscript:
Sm <- getPerms(1:m)
np <- dim(Sm)[1]

#Same as vector of strings:
sSm <- apply(Sm, 1, function(x) paste(x, collapse=""))
#sSm
#"123" "132" "213" "231" "312" "321"


#Covariance matrix from Weiß (2022):
Sigma3 <- rbind(c(46,-23,-23,7,7,-14),c(-23,28,10,-20,-2,7),c(-23,10,28,-2,-20,7),c(7,-20,-2,28,10,-23),c(7,-2,-20,10,28,-23),c(-14,7,7,-23,-23,46)) / 360

#Corresponding eigenvalues:
eigen3 <- c( (2+sqrt(2))/12, 2/15, 1/10, (2-sqrt(2))/12 )



#Compute asymptotics first:
#install.packages("CompQuadForm")
library("CompQuadForm")

level <- 0.05

qup <- uniroot(function(x) level-davies(x, eigen3, lim=1e7, acc=1e-8)$Qq, c(1,2))$root
qup #1.484225

#Function for computing P-values:
Pvalue <- function(stat) davies(stat, eigen3, lim=1e7, acc=1e-8)$Qq




#############################
#Tests for serial dependence:
#############################

reps <- 1e4

tabT <- c(100, 250, 500, 1000) #
nT <- length(tabT)

Hmax <- log(mfak)
D0 <- 1/D.max(mfak)

#AR parameter:
a1 <- 0.5
sig <- sqrt(1/(1-a1^2))

#Ma parameter:
b1 <- 0.8

#TEAR(1) parameter:
prerun <- 100
ea1 <- 0.15

set.seed(123)

for(iT in 1:nT){
	T <- tabT[iT]

	tabDecisions <- rep(0, 6)

for(r in 1:reps){
	#H0: iid
	data <- rnorm(T)
	
	# #H1: AR(1)
	# X <- rnorm(1,0,sig)
	# eps <- rnorm(T)
	# data <- rep(NA, T)
	# for(t in 1:T){
		# X <- a1*X+eps[t]
		# data[t] <- X
	# }

	# #H1: QMA(1)
	# eps <- rnorm(T+1)
	# data <- eps[-1]+b1*eps[1:T]^2
	
	# #H1: TEAR(1)
	# eps <- rexp(prerun, 1)
	# tabR <- rbinom(prerun, 1, ea1)
	# X <- 1
	# for(t in 1:prerun) X <- (1-ea1)*eps[t] + tabR[t]*X
	
	# eps <- rexp(T, 1)
	# tabR <- rbinom(T, 1, ea1)
	# data <- rep(NA, T)
	# for(t in 1:T){
		# X <- (1-ea1)*eps[t] + tabR[t]*X
		# data[t] <- X
	# }
	
	
	for(d in tabd){
		n <- T-(m-1)*d

		data0 <- data[1:n]
		for(k in 1:(m-1)){
			data0 <- cbind(data0, data[(k*d+1):(k*d+n)])
		}

		#Generate ordinal patterns as rank vectors:
		opts <- t(matrix(c(apply(data0, 1, rank, ties.method="first")), nrow=m))

		#rm(data0)

		#transform into vector of strings:
		sopts <- apply(opts, 1, function(x) paste(x, collapse=""))

		#rm(opts)

		#So sopts is categorical time series with range sSm.

		#Numeric coding by 1-m!:
		datanum <- match(sopts, sSm)

		# #Binarization:
		# bincodes <- diag(1,np)
		# databin <- bincodes[datanum,] #assign row vectors
		
		# hatp <- colMeans(databin)
		
		hatp <- (as.vector(table(c(datanum, 1:mfak))) - 1) / n
		
		#Statistics:
		Hstat <- H(hatp)
		Dstat <- D(hatp)
		Cstat <- C(hatp)
		
		#Test decisions:
		if(n/mfak*(Hmax-Hstat+4*Dstat) > qup){
			tabDecisions[1] <- tabDecisions[1] + 1
		}
		if(n*8/5/mfak*(Hmax-Hstat+Dstat) > qup){
			tabDecisions[2] <- tabDecisions[2] + 1
		}
		if(n*8/sqrt(17)/mfak*sqrt((Hmax-Hstat)^2+Dstat^2) > qup){
			tabDecisions[3] <- tabDecisions[3] + 1
		}
		
		if(n/mfak*(Hmax-Hstat+4*Cstat/D0) > qup){
			tabDecisions[4] <- tabDecisions[4] + 1
		}
		if(n*8/mfak/(4/Hmax+D0)*(1-Hstat/Hmax+Cstat) > qup){
			tabDecisions[5] <- tabDecisions[5] + 1
		}
		if(n*8/mfak/sqrt((4/Hmax)^2+D0^2)*sqrt((1-Hstat/Hmax)^2+Cstat^2) > qup){
			tabDecisions[6] <- tabDecisions[6] + 1
		}
	}#for d
		
}#for reps

print(c(T, tabDecisions))
	
}#for T

#H0: iid
# 100 549 549 553 534 534 534
# 250 485 481 479 478 478 478
# 500 509 506 507 505 504 504
# 1000  502  502  505  501  501  501

#H1: AR(1)
#  100 2563 2572 2575 2524 2524 2530
#  250 6146 6174 6184 6110 6110 6114
#  500 9331 9335 9335 9328 9328 9328
# 1000 9997 9997 9997 9997 9997 9997

#H1: QMA(1)
#  100 1601 1599 1594 1558 1558 1559
#  250 3433 3437 3443 3411 3411 3414
#  500 6036 6049 6053 6014 6015 6015
# 1000 8861 8867 8869 8859 8859 8859

#H1: TEAR(1)
#  100 2478 2482 2485 2435 2435 2436
#  250 5265 5267 5267 5241 5241 5241
#  500 8347 8354 8355 8333 8333 8333
# 1000 9882 9883 9883 9882 9882 9882


