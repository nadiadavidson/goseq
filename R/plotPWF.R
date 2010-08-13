#Function for plotting the PWF to exhibit length bias
plotPWF = function(pwf,binsize="auto"){
		w=!is.na(pwf$bias.data)
		o=order(pwf$bias.data[w])
		if(binsize=="auto"){
			binsize=max(1,min(100,floor(sum(w)*.08)))
			rang=max(pwf$pwf,na.rm=TRUE)-min(pwf$pwf,na.rm=TRUE)
			resid=rang
			oldwarn=options()$warn
			options(warn=-1)
			while(binsize<=floor(sum(w)*.1) & resid/rang>.001){
				binsize=binsize+100
				splitter=ceiling(1:length(pwf$DEgenes[w][o])/binsize)
				de=sapply(split(pwf$DEgenes[w][o],splitter),mean)
				binlen=sapply(split(as.numeric(pwf$bias.data[w][o]),splitter),mean)
				resid=sum((de-approx(pwf$bias.data[w][o],pwf$pwf[w][o],binlen)$y)^2)/length(binlen)
			}
			options(warn=oldwarn)
		}else{
			splitter=ceiling(1:length(pwf$DEgenes[w][o])/binsize)
			de=sapply(split(pwf$DEgenes[w][o],splitter),mean)
			binlen=sapply(split(as.numeric(pwf$bias.data[w][o]),splitter),mean)
		}
		plot(binlen,de,xlab=paste("Bias Data in ",binsize," gene bins.",sep=""),ylab="Proportion DE")
		lines(pwf$bias.data[w][o],pwf$pwf[w][o],col=3,lwd=2)
}
