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

	#Make a suitable object
	out=data.frame(DEgenes=DEgenes,bias.data=bias.data,pwf=pwf,stringsAsFactors=FALSE)
	rownames(out)=names(DEgenes)

	#Plot
	if(plot.fit){
		plotPWF(out)
	}

	#Warn if we couldn't find length data for at least 80% of genes.  This is also in goseq function, so leave it out
	#nafrac=(sum(is.na(pwf))/length(pwf))*100
	#if(nafrac>.25){
	#	warning(paste("Missing length data for ",round(nafrac*100),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
	#}

	out
}


