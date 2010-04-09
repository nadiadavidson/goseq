nullp = function(DEgenes, genome, id, bias.data=NULL, plot.fit=TRUE){
	#Generates the Null distribution

	#Input Checking
	if(!is.null(bias.data) & length(bias.data)!=length(DEgenes)){
		stop("bias.data vector must have the same length as DEgenes vector!")
	}
	bias.data=unfactor(bias.data)


	#Look it up if necessary
	if(is.null(bias.data)){
		bias.data=getlength(names(DEgenes),genome,id)
	}

	#Do the fit.
	#May not have bias data for some of the entries, return NA at the relevant position
	pwf=rep(NA,length(DEgenes))
	w=!is.na(bias.data)
	pwf[w]=makespline(bias.data[w],DEgenes[w])

	#Plot
	if(plot.fit){
		o=order(bias.data[w])
		j=min(100,floor(sum(w)*.08))
		resid=1
		rang=max(pwf,na.rm=TRUE)-min(pwf,na.rm=TRUE)
		oldwarn=options()$warn
		options(warn=-1)
		while(j<=floor(sum(w)*.1) & resid/rang>.001){
			binsize=j
			max=floor(sum(w)/binsize)
			de=vector()
			binlen=vector()
			for(i in 1:max){
				a=sum(DEgenes[w][o][(1+(i-1)*binsize):(i*binsize)])
				binlen=c(binlen,median(as.numeric(bias.data[w][o][(1+(i-1)*binsize):(i*binsize)])))
				de=c(de,a/binsize)
			}
			resid=sum((de-approx(bias.data[w][o],pwf[w][o],binlen)$y)^2)/length(binlen)
			j=j+100
		}
		options(warn=oldwarn)
		plot(binlen,de,xlab=paste("Bias Data in ",binsize," gene bins.",sep=""),ylab="Proportion DE")
		lines(bias.data[w][o],pwf[w][o],col=3,lwd=2)
	}

	#Warn if we couldn't find length data for at least 80% of genes.  This is also in goseq function, so leave it out
	nafrac=(sum(is.na(pwf))/length(pwf))*100
	#if(nafrac>.25){
	#	warning(paste("Missing length data for ",round(nafrac*100),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
	#}

	#Return
	pwf
}


