makespline=function (x, y, newX=NULL, nKnots = 6) 
{
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
   f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"),family=binomial())
    dat <- data.frame(x = x, y = y)
    sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
    if (length(sm$xp) < 6) 
        warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
    F <- mono.con(sm$xp,TRUE,lower=0,upper=1)
    G <- list(X = sm$X, C = matrix(0, 0, 0), sp = f.ug$sp, p = sm$xp, 
        y = y, w = y * 0 + 1, Ain =F$A, bin = F$b, S = sm$S, 
        off = 0)
    p <- pcls(G)
    fv <- Predict.matrix(sm, data.frame(x = newX)) %*% p
    fv <- as.vector(fv)
   return(fv)
}
