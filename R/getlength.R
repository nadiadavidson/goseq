getlength=function(genes,genome,id){
	#Test if the package exists
	if(!any(library()$results[,"Package"]=="geneLenDataBase")){
		stop("Package geneLenDataBase is required for automatic gene length lookup.")
	}
	#The local repository of Length stuff
	repo=grep(".*\\..*\\.LENGTH",as.data.frame(data(package="geneLenDataBase")$results,stringsAsFactors=FALSE)$Item,ignore.case=TRUE,value=TRUE)
	#Is this genome/id in the local repository?
	if(paste(genome,id,"LENGTH",sep='.')%in%repo){
		message("Loading ",genome," length data...")
		data(list=paste(genome,id,"LENGTH",sep='.'),package='geneLenDataBase')
		data=split(get(paste(genome,id,"LENGTH",sep='.'))$Length,get(paste(genome,id,"LENGTH",sep='.'))$Gene)
	}else{
		#If not, then the user will have to supply it themselves
		stop("Length information for genome ",genome," and gene ID ",id," is not in the geneLenDataBase database.  You will have to specify bias.data manually.")
	}
	len=data.frame(Gene=names(data),Median=sapply(data,median),Min=sapply(data,min),Max=sapply(data,max),Count=sapply(data,length))
	return(len$Median[match(genes,len$Gene)])
}
