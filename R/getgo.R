getgo=function(genes,genome,id,fetch.cats=c("GO:CC","GO:BP","GO:MF")){
	if(any(!fetch.cats%in%c("GO:CC","GO:BP","GO:MF","KEGG"))){
		stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
	}
	#Convert from genome ID to org.__.__.db format
	orgstring=as.character(.ORG_PACKAGES[grep(gsub("[0-9]+",'',genome),names(.ORG_PACKAGES),ignore.case=TRUE)])
	#Multimatch or no match
	if(length(orgstring)!=1){
		stop("Couldn't grab GO categories automatically.  Please manually specify.")
	}
	library(paste(orgstring,"db",sep='.'),character.only=TRUE)
	coreid=strsplit(orgstring,"\\.")[[1]][3]

	#Now we need to convert it into the currently used ID type
	idname=as.character(.ID_MAP[grep(id,names(.ID_MAP),ignore.case=TRUE)])
	#Multimatch or no match
	if(length(idname)!=1){
		stop("Couldn't grab GO categories automatically.  Please manually specify.")
	}

	#Most of the core IDs are eg, for now assume that is the case
	#This really needs to be re-written.  We want two IDs, coreid and sourceid and then have two routines.  One for when coreid==sourceid and another for coreid!=sourceid
	if(idname=="Entrez Gene ID"){
		if(coreid=="eg"){
			go=list()
			if(length(grep("^GO",fetch.cats))!=0){
				x=as.list(get(paste(orgstring,"GO",sep='')))
				tmp=unlist(x)
				goids=tmp[grep("GOID$",names(tmp))]
				type=tmp[grep("Ontology$",names(tmp))]
				keep=goids[as.character(type)%in%gsub("^GO:",'',fetch.cats)]
				go=split(as.character(keep),gsub("\\.GO.*",'',names(keep)))
			}
			if(length(grep("^KEGG",fetch.cats))!=0){
				x=as.list(get(paste(orgstring,"PATH",sep='')))
				if(length(go)==0){
					curr.go=NULL
				}else{
					curr.go=unlist(go,use.names=FALSE)
					names(curr.go)=rep(names(go),times=as.numeric(summary(go)[,1]))
				}
				new.go=unlist(x,use.names=FALSE)
				names(new.go)=rep(names(x),times=as.numeric(summary(x)[,1]))
				tot=c(new.go,curr.go)
				go=split(as.character(tot),names(tot))
			}
		}else{
			stop("Couldn't grab GO categories automatically.  Please manually specify.")
		}
	}else if(idname=="Ensembl gene ID"){
		if(coreid=="eg"){
			#Super awesome-o efficient, but not readable
			en2eg=as.list(get(paste(orgstring,"ENSEMBL2EG",sep='')))
			eg2go=list()
			if(length(grep("^GO",fetch.cats))!=0){
				x=as.list(get(paste(orgstring,"GO",sep='')))
				tmp=unlist(x)
				goids=tmp[grep("GOID$",names(tmp))]
				type=tmp[grep("Ontology$",names(tmp))]
				keep=goids[as.character(type)%in%gsub("^GO:",'',fetch.cats)]
				eg2go=split(as.character(keep),gsub("\\.GO.*",'',names(keep)))
			}
			if(length(grep("^KEGG",fetch.cats))!=0){
				x=as.list(get(paste(orgstring,"PATH",sep='')))
				if(length(go)==0){
					curr.go=NULL
				}else{
					curr.go=unlist(go,use.names=FALSE)
					names(curr.go)=rep(names(go),times=as.numeric(summary(go)[,1]))
				}
				new.go=unlist(x,use.names=FALSE)
				names(new.go)=rep(names(x),times=as.numeric(summary(x)[,1]))
				tot=c(new.go,curr.go)
				eg2go=split(as.character(tot),names(tot))
			}
			#eg2go=lapply(as.list(get(paste(orgstring,"GO",sep=''))),names)
			len_eg2go=as.numeric(summary(eg2go)[,1])
			len_en2eg=as.numeric(summary(en2eg)[,1])
			flat_en2eg=unlist(en2eg,use.names=FALSE)
			names(flat_en2eg)=rep(names(en2eg),times=len_en2eg)
			w=match(flat_en2eg,names(eg2go))
			len_eg=len_eg2go[w]
			len_eg[is.na(len_eg)]=0
			names(len_eg)=rep(names(en2eg),times=len_en2eg)
			len_en=split(as.numeric(len_eg),names(len_eg))
			len_en=sapply(len_en,sum,na.rm=TRUE)
			flat_eg2go=unlist(eg2go[w],use.names=FALSE)
			go=split(flat_eg2go,rep(names(len_en),times=as.numeric(len_en)))
			go=lapply(go,unique)
			#Obvious but slow way, probably better for the vignette
			#x=as.list(get(paste(orgstring,"ENSEMBL2EG",sep='')))[genes]
			#mapkeys=lapply(as.list(get(paste(orgstring,"GO",sep='')))[unlist(x,use.names=FALSE)],names)
			#grepGO=function(id, mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
			#go=lapply(x,grepGO,mapkeys)
		}else{
			stop("Couldn't grab GO categories automatically.  Please manually specify.")
		}
	}else if(idname=="Gene Symbol"){
		if(coreid=="eg"){
			#Super awesome-o efficient, but not readable
			en2eg=as.list(get(paste(orgstring,"SYMBOL2EG",sep='')))
			eg2go=list()
			if(length(grep("^GO",fetch.cats))!=0){
				x=as.list(get(paste(orgstring,"GO",sep='')))
				tmp=unlist(x)
				goids=tmp[grep("GOID$",names(tmp))]
				type=tmp[grep("Ontology$",names(tmp))]
				keep=goids[as.character(type)%in%gsub("^GO:",'',fetch.cats)]
				eg2go=split(as.character(keep),gsub("\\.GO.*",'',names(keep)))
			}
			if(length(grep("^KEGG",fetch.cats))!=0){
				x=as.list(get(paste(orgstring,"PATH",sep='')))
				if(length(go)==0){
					curr.go=NULL
				}else{
					curr.go=unlist(go,use.names=FALSE)
					names(curr.go)=rep(names(go),times=as.numeric(summary(go)[,1]))
				}
				new.go=unlist(x,use.names=FALSE)
				names(new.go)=rep(names(x),times=as.numeric(summary(x)[,1]))
				tot=c(new.go,curr.go)
				eg2go=split(as.character(tot),names(tot))
			}
			#eg2go=lapply(as.list(get(paste(orgstring,"GO",sep=''))),names)
			len_eg2go=as.numeric(summary(eg2go)[,1])
			len_en2eg=as.numeric(summary(en2eg)[,1])
			flat_en2eg=unlist(en2eg,use.names=FALSE)
			names(flat_en2eg)=rep(names(en2eg),times=len_en2eg)
			w=match(flat_en2eg,names(eg2go))
			len_eg=len_eg2go[w]
			len_eg[is.na(len_eg)]=0
			names(len_eg)=rep(names(en2eg),times=len_en2eg)
			len_en=split(as.numeric(len_eg),names(len_eg))
			len_en=sapply(len_en,sum,na.rm=TRUE)
			flat_eg2go=unlist(eg2go[w],use.names=FALSE)
			go=split(flat_eg2go,rep(names(len_en),times=as.numeric(len_en)))
			go=lapply(go,unique)
			#Obvious but slow way, probably better for the vignette
			#x=as.list(get(paste(orgstring,"ENSEMBL2EG",sep='')))[genes]
			#mapkeys=lapply(as.list(get(paste(orgstring,"GO",sep='')))[unlist(x,use.names=FALSE)],names)
			#grepGO=function(id, mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
			#go=lapply(x,grepGO,mapkeys)

		}else{
			stop("Couldn't grab GO categories automatically.  Please manually specify.")
		}
	}else if(idname=="HAVANA Pseudogene ID"){
		if(coreid=="eg"){
			stop("Couldn't grab GO categories automatically.  Please manually specify.")
		}else{
			stop("Couldn't grab GO categories automatically.  Please manually specify.")
		}
	}else{
		stop("Couldn't grab GO categories automatically.  Please manually specify.")
	}
	#Now look them up
	gene2cat=go[genes]
	return(gene2cat)
}
