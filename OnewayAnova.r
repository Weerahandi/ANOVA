
# Oneway ANOVA under unqual Variances
# source as source("onewayanova.r")
# run as in # out=oneway.anova(extra~group,sleep)

oneway.anova<-function(formula, data, nM=1000000) 
{
    	if (missing(formula) || (length(formula) != 3L)) 
        	stop("'formula' missing or incorrect")
    	dp = as.character(formula)
    	if (length(dp) != 3L) 
        	stop("a two-sided formula is required")
    	regname = paste(dp[[2L]], "on", dp[[3L]])
    	m = match.call(expand.dots = FALSE)
    	if (is.matrix(eval(m$data, parent.frame()))) 
        	m$data = as.data.frame(data)
    	m$var.equal = NULL
    	m[[1L]] = quote(stats::model.frame)
    	mf = eval(m, parent.frame())
    	response = attr(attr(mf, "terms"), "response")
    	y = mf[[response]]
    	if (length(mf[-response]) > 1L) 
        g = factor(do.call("interaction", mf[-response]))
    	else g = factor(mf[[-response]])
    	k = nlevels(g)
	k1 = k - 1 ; k21 = k^2-1

    	if (k < 2L) 
        stop("not enough groups")
   	ns = tapply(y, g, length)
    	if (any(ns < 2)) 
        stop("not enough observations")
    	rs = ns-1
    	xb = tapply(y, g, mean) # Sample means
    	sq = tapply(y, g, var) # Sample variances
  	wts = ns/sq
    	xbs = sum(ns*xb/sq)/sum(ns/sq) # TO CHECK
	rhs = sum(wts*(xb-xbs)^2)
# Save summary Stats
	Stats = cbind.data.frame(xb,sq,ns)
	names(Stats) = c("Sample Means", "Sample Variances", "Sample Sizes")
# Compute Welch test
	w = sum(wts)
	wpart = (1-wts/w)^2
	Fnume =  sum(wts*(xb-xbs)^2)/k1
	Fpart = sum((1/rs)*wpart)
	Fdeno = 1 + 2*(k-2)*Fpart/k21
	Fstat = Fnume/Fdeno
	df = k21/3 /sum(wpart/rs)
	wpval = 1-pf(Fstat,k1,df)	
	
# Compue Genralized Test using random T variables: Generalized Test Equivalent to Fiducial Test
	Deno = w # sum(ns/sq)
	sqw = sqrt(wts)
# Compue Genralized Test using random Chi and Z variables: Generalized Test Equivalent to PB Test
	sqn = sqrt(ns)
	sqr = sqrt(rs)
	T1 = T2 = P1 = P2 = P3 = 0

	for (i in 1:k) {
		ztmp = rnorm(nM)
		Zvec = (ztmp - mean(ztmp))/sd(ztmp) # Bias corrected random numbers
		Chivec = rchisq(nM, rs[i])
		Chivec = rs[i] * Chivec / mean(Chivec) # Chi-aqured has only one parameter to be corrected	
		Tvec = sqr[i]*Zvec /sqrt(Chivec)
		Tvec2 = Tvec^2 # Square of T
# Comput Gen. Test based on Tvec
		T1 = T1 + Tvec2
		T2 = T2 + Tvec*sqw[i]	
		P2 = P2 + sqw[i]*Zvec*rs[i] / Chivec
		P3 = P3 + wts[i]*rs[i]/ Chivec
	}

	Del = T1 - T2^2/Deno
	g1pval = mean(Del > rhs)

	P1 = T1
	pb = P1 - P2^2 / P3

	g2pval = mean(pb > rhs)

	pvals = round(c(wpval, g1pval,g2pval),3)
	methods = c("Welch Test (approximate)", "Generalized Test Equivalent to Fiducial Test (size-assured, but conservative)",  
               	"Generalized Test Equivalent to Parametric Bootstrap Test (size close to intended)")
	pVals = cbind.data.frame(methods, pvals)
	names(pVals) = c("Method", "p-value")
	print(pVals)
	outlist=list(Stats=Stats, pValues=pVals)	
	return(outlist)
}
