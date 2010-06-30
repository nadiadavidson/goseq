goseq=function(DEgenes,pwf,genome,id,gene2cat=NULL,test.cats=c("GO:CC","GO:BP","GO:MF"),method="Wallenius",repcnt=2000){
	#Gene Ontology testing
	#Do some validation of input variables
	if(any(!test.cats%in%c("GO:CC","GO:BP","GO:MF","KEGG"))){
		stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
	}
	if((missing(genome) | missing(id))){
		if(is.null(gene2cat)){
			stop("You must specify the genome and gene ID format when automatically fetching gene to GO category mappings.")
		}
		genome='dummy'
		id='dummy'
	}
	if(!any(method%in%c("Wallenius","Sampling","Hypergeometric"))){
		stop("Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric.")
	}
	if(length(pwf)!=length(DEgenes)){
		stop("The PWF vector must have the same length as DEgenes vector!")
	}
	if(!is.null(gene2cat) && (!is.data.frame(gene2cat) & !is.list(gene2cat))){
		stop("Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again.")
	}
	if(missing(pwf)){
		if(method=="Hypergeometric"){
			pwf=rep(1,length(DEgenes))
		}else{
			stop("The PWF vector must be supplied!")
		}
	}
		
	#Factors are evil
	pwf=unfactor(pwf)
	DEgenes=unfactor(DEgenes)
	gene2cat=unfactor(gene2cat)

	#Fetch the GO categories if necessary, otherwise do the necessary formatting
	if(is.null(gene2cat)){
		#When we fetch the data using getgo it will be in the list format
		#Fetch the data
		message("Fetching GO annotations...")
		gene2cat=getgo(names(DEgenes),genome,id,fetch.cats=test.cats)
		names(gene2cat)=names(DEgenes)
		#Do the two rebuilds to remove any nulls
		cat2gene=reversemapping(gene2cat)
		gene2cat=reversemapping(cat2gene)
	}else{
		#The user is asked to input categories in the flat data.frame format, so need to convert to list
		message("Using manually entered categories.")
		#Is it a flat mapping?
		if(class(gene2cat)!="list"){
			#Work out which column contains the genes
			genecol=as.numeric(apply(gene2cat,2,function(u){sum(u%in%names(DEgenes))}))
			genecol=which(genecol!=0)
			if(length(genecol)>1){
				genecol=genecol[order(-genecol)[1]]
				warning(paste("More than one possible gene column found in gene2cat, using the one headed",colnames(gene2cat)[genecol]))
			}
			if(length(genecol)==0){
				genecol=1
				warning(paste("Gene column could not be identified in gene2cat conclusively, using the one headed",colnames(gene2cat)[genecol]))
			}
			othercol=1
			if(genecol==1){othercol=2}
			#Now put it into our delicious listy format
			gene2cat=split(gene2cat[,othercol],gene2cat[,genecol])
			#Do the appropriate builds
			cat2gene=reversemapping(gene2cat)
			gene2cat=reversemapping(cat2gene)
		}
		#Is the variable specified actually cat2gene?  If so convert it.
		if(sum(unique(unlist(gene2cat,use.names=FALSE))%in%names(DEgenes))>sum(unique(names(gene2cat))%in%names(DEgenes))){
			gene2cat=reversemapping(gene2cat)
		}
		#Alright, we're garunteed a list now
		#Now we have to throw out any genes which aren't in req
		gene2cat=gene2cat[names(gene2cat)%in%names(DEgenes)]
		#Rebuild because it's a fun thing to do
		cat2gene=reversemapping(gene2cat)
		gene2cat=reversemapping(cat2gene)
	}
	
	#Fill in the NA values with the median
	w=!is.na(pwf)
	common=pwf[w][order(pwf[w])][floor(length(pwf[w])/2)]
	nafrac=(sum(is.na(pwf))/length(pwf))*100
	#if(nafrac>25){
	#	warning(paste("Missing length data for ",round(nafrac),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
	#}
	pwf[is.na(pwf)]=common

	#Now do the calculation of p-values
	#Always need these, so build them
	cats=names(cat2gene)
	n=sum(DEgenes)
	N=length(DEgenes)
	de=names(DEgenes)[as.logical(DEgenes)]
	stop=length(cats)
	gpval=data.frame(category=cats,upval=0,dpval=0,stringsAsFactors=FALSE)
	if(method=="Sampling"){
		pvalmask=rep(0,length(cat2gene))
		a=table(unlist(gene2cat[de],FALSE,FALSE))
		pvalmask[match(names(a),cats)]=a
		pvalmask=as.integer(pvalmask)
		#Need to do this to ensure NULL category entries for "empty" genes and to ensure order matches pwf
		gtg=gene2cat
		gtg=gtg[names(DEgenes)]
		names(gtg)=names(DEgenes)
		gtg[names(gene2cat)]=gene2cat[names(gene2cat)]
		message("Running the simulation...")
		lookup=matrix(0,nrow=repcnt,ncol=length(cat2gene))
		for(i in 1:repcnt){
			#Will return character(0) if the list is empty, which should then go through without altering anything
			a=table(as.character(unlist(gtg[order(runif(N)^(1/pwf),decreasing=TRUE)[1:n]],FALSE,FALSE)))
			lookup[i,match(names(a),cats)]=a
			pp(repcnt)
		}
		message("Calculating the p-values...")
		for(i in 1:stop){
			gpval[i,2:3]=c((sum(lookup[,i]>=pvalmask[i])+1)/(repcnt+1),(sum(lookup[,i]<=pvalmask[i])+1)/(repcnt+1))
			pp()
		}
	}else if(method=="Wallenius"){
		message("Calculating the p-values...")
		degenesnum=which(DEgenes==1)
		cat2genenum=relist(match(unlist(cat2gene),names(DEgenes)),cat2gene)
		alpha=sum(pwf)
		gpval[,2:3]=unlist(lapply(cat2genenum,function(u){
			subcount=sum(degenesnum%in%u)
			totcount=length(u)
			a=mean(pwf[u])
			weight=(a*(N-totcount))/(alpha-totcount*a)
			c(dWNCHypergeo(subcount,totcount,N-totcount,n,weight)+pWNCHypergeo(subcount,totcount,N-totcount,n,weight,lower.tail=FALSE),pWNCHypergeo(subcount,totcount,N-totcount,n,weight))
			}))[c(which(1:(2*stop)%%2==1),which(1:(2*stop)%%2==0))]
	}else if(method=="Hypergeometric"){
		message("Calculating the p-values...")
		degenesnum=which(DEgenes==1)
		cat2genenum=relist(match(unlist(cat2gene),names(DEgenes)),cat2gene)
		gpval[,2:3]=unlist(lapply(cat2genenum,function(u){
			subcount=sum(degenesnum%in%u)
			totcount=length(u)
			c(dhyper(subcount,totcount,N-totcount,n)+phyper(subcount,totcount,N-totcount,n,lower.tail=FALSE),phyper(subcount,totcount,N-totcount,n))
			}))[c(which(1:(2*stop)%%2==1),which(1:(2*stop)%%2==0))]
	}
	gpval=gpval[order(gpval$upval),]
	return(gpval)
}
