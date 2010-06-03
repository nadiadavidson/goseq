makespline=function (x, y, newX=NULL, nKnots = 6) 
{
    #Should not be used outside of goseq package.  Refer to the help pages for pcls in the mgcv package for more general 
    #contstrained spline fitting.
	ww=order(x)
	size=ceiling(length(y)/10)
	low=sum(y[ww][1:size])
	hi=sum(y[ww][(length(y)-size):length(y)])
	#Is the trend decreasing
	if(hi<=low){
		#print("Decreasing")
		#Reform as a montonicly increasing fit
		x=10^10-x
	}
    nKnots <- round(nKnots)
    if (is.null(newX)) 
        newX <- x
    #We need to force the fit to go through (0,0), mono.con constructs a requirment that the LOWEST VALUE go through 0
    #In order to make the lowest value 0, we add in a dummy point.  Hopefully this does not skew things too greatly
    #We use an inequality constraint rather than equality to only skew the fit enough that the spline is always >0
    #This won't work if we have a monotonically decreasing function, so check that
    if(hi>low){
	    x=c(0,x)
	    y=c(0,y)
    }
    f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"),family=binomial())
    dat <- data.frame(x = x, y = y)
    sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
    if (length(sm$xp) < 6) 
        warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
    #The lower and upper parameters here don't seem to build the right constraints 100% of the time, so add them manually
    F <- mono.con(sm$xp,TRUE)
    F$A=rbind(F$A,c(1,rep(0,ncol(F$A)-1)),c(rep(0,ncol(F$A)-1),-1))
    F$b=c(F$b,0,1)
    G <- list(X = sm$X, C = matrix(0, 0, 0), sp = f.ug$sp, p = sm$xp, 
        y = y, w = y * 0 + 1, Ain =F$A, bin = F$b, S = sm$S, 
        off = 0)
    p <- pcls(G)
    fv <- Predict.matrix(sm, data.frame(x = newX)) %*% p
    fv <- as.vector(fv)
   return(fv)
}
