#############################################################################
#Description: Fetches gene length data for the genome and id specified.
#Notes:
#Author: Matthew Young
#Date Modified: 13/12/2010

getlength=function(genes,genome,id){
	#geneLenDataBase contains all the gene length data, check it is installed
	if(!any(library()$results[,"Package"]=="geneLenDataBase")){
		stop("Package geneLenDataBase is required for automatic gene length lookup.")
	}
	#Get a list of the length data files that are in the package
	repo=grep(".*\\..*\\.LENGTH",as.data.frame(data(package="geneLenDataBase")$results,stringsAsFactors=FALSE)$Item,ignore.case=TRUE,value=TRUE)
	#Is this genome/id in the local repository?
	if(paste(genome,id,"LENGTH",sep='.')%in%repo){
		#Yes it is, so load it and group transcript lengths into genes
		message("Loading ",genome," length data...")
		data(list=paste(genome,id,"LENGTH",sep='.'),package='geneLenDataBase')
		data=split(get(paste(genome,id,"LENGTH",sep='.'))$Length,get(paste(genome,id,"LENGTH",sep='.'))$Gene)
	}else{
		#If not, then the user will have to supply it themselves
		stop("Length information for genome ",genome," and gene ID ",id," is not in the geneLenDataBase database.  You will have to specify bias.data manually.")
	}
	#In order to assign just one length to each gene, we store the Median, Min Max and count of transcripts for each gene for later use
	len=data.frame(Gene=names(data),Median=sapply(data,median),Min=sapply(data,min),Max=sapply(data,max),Count=sapply(data,length))
	#Get the names of the genes used in the genome package
	gene_names=get(paste(genome,id,"LENGTH",sep="."))$Gene
	#How many of them match the names the user specified with "genes"
	matched_frac=sum(genes%in%gene_names)/length(genes)
	#None? Throw an error.
	if(matched_frac==0){
		stop("The gene names specified do not match the gene names for genome ",genome," and ID ",id,".\n\tGene names given were: ",paste(genes[1:10],collapse=", "),"\n\tRequired gene names are: ",paste(gene_names[1:10],collapse=", "))
	}
	#Less than 60%? Give a warning and continue
	if(matched_frac<.6){
		warning("More than 40% of gene names specified did not match the gene names for genome ",genome," and ID ",id,".  No length data will be available for these genes.\n\tGene names which failed to match were: ",paste(genes[!genes%in%gene_names][1:10],collapse=", "),"\n\tRequired gene names are: ",paste(gene_names[1:10],collapse=", "))
	}
	#Finally, return the median transcript length for each matching gene
	return(len$Median[match(genes,len$Gene)])
}
