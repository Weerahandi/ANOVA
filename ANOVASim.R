# Simulation Study for ANOVA under unequal variances
# Compile as source("AnovaSim.R")
# Run as # out=simulate(c(2,3,2),c(1,2,3),5000)
# Or as out=simulate(c(10,10),c(1,2),1000)
# Or as out=simulate(c(4,4,4,3),c(7,41,38,430),1000)

# Test performance of Competing Tests
# n is the vector of sample sizes, VarE is the vector of error variances
simulate=function(n=c(2,3,2),VarE=c(1,1,1),nr=5000)
{
# set.seed(4230)
	k=ng=length(n)
	k1 = k - 1 ; k21 = k^2-1
	r=n-1 # degrees of freedom
	sqn = sqrt(n)
	sig = sqrt(VarE)
	
	mu=1 # common mean is also not important 
	pvalsW=pvalsF=pvalsP=rep(0,nr)
# Pre-generate data as Kris did
	XB = matrix(0, nr,k); SQ = matrix(0, nr, k)
	z = matrix(0, nr, k); chi = matrix(0, nr, k); 
# Generate observed sampe means and smaple sds as Kris did
	T=matrix(0,nr,k)
for(i in 1:k){ 
	SQ[,i] = VarE[i]*rchisq(nr, r[i])/r[i]
	XB[,i] = mu + rnorm(nr)*sig[i]/sqn[i]
# Also generate random numbers at the same time
	z[,i] = rnorm(nr)
	chi[, i] = rchisq(nr, r[i])
	T[,i] = rt(nr,r[i])
}
	sd = sqrt(SQ)

	nroot=sqrt(n)
# Generate sample means and sample varinces
	Tobs=rep(0,nr)
     	for (m in 1:nr) { # We are doing one observed data at a time and using independent T raadom values
		xb=XB[m,] #	rnorm(ng,mu,sqrt(VarE/n))# Sample means for ng groups
		s2=SQ[m,] #	VarE*rchisq(ng,r)/r	# Sample variances for ng groups
		sd=sqrt(s2)
		wts=nos2=n/s2
		w=sumnos2=sum(nos2)
		sqnos2=sqrt(nos2)
		xbs=sum(nos2*xb)/sum(nos2) # common sample mean
		Tobs[m]=sum(nos2*(xb-xbs)^2) # Observed statistic
# Compute p-value using a sample of 1,000 random numbers from each distribution	
		M=1000 # Montecarlo samples to evaluate p-value; this need not be large
		ind=sample(c(1:nr),M)
# Using Kris-like method to minimize computing time compute p-values
		hitsF=hitsF=hitsP=rep(0,M)
# First Compute Welch test, requiring no Monte Carlo samples
		wpart = (1-wts/w)^2
		Fnume =  sum(wts*(xb-xbs)^2)/k1
		Fpart = sum((1/r)*wpart)
		Fdeno = 1 + 2*(k-2)*Fpart/k21
		Fstat = Fnume/Fdeno
		df = k21/3 /sum(wpart/r)
		wpval = 1-pf(Fstat,k1,df)
		for (i in 1:M) {	
# Next compute Fiducial test
			Tvec=T[ind[i],] 
			Ttest=sum(Tvec^2)-sum(Tvec*sqnos2)^2/sumnos2
			if (Ttest>Tobs[m]) hitsF[i]=1	
# Now compute p-value of PB test
			Zvec=	z[ind[i],]	
			Chivec=chi[ind[i],]
			Ptest=sum(Zvec^2*r/Chivec)-(sum(sqn*Zvec*r/sd/Chivec))^2/sum(n*r/s2/Chivec)
			if (Ptest>Tobs[m]) hitsP[i]=1	
		}
		pvalsW[m]=wpval
		pvalsF[m]=mean(hitsF) 
		pvalsP[m]=mean(hitsP)
	}
print("Type I error of Welch, Fiducial and PB Tests when intended size os .05")
print(mean(pvalsW<.05))
print(mean(pvalsF<.05))
print(mean(pvalsP<.05))

return(cbind.data.frame(pvalsW,pvalsF,pvalsP))	
}

