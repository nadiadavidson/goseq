#  ZZZ.R

.onAttach <- function(libname, pkgname)
{
	if(.Platform$OS.type=="windows" && .Platform$GUI=="Rgui" ) {
		winMenuAddItem("Vignettes","goseq","shell.exec(system.file(\"doc\",\"goseq.pdf\",package=\"goseq\"))")
	}
}
#Mappings that are primarily required by getgo
.ORG_PACKAGES=paste("org.",c("Ag.eg","At.tair","Bt.eg","Ce.eg","Cf.eg","Dm.eg","Dr.eg","EcK12.eg","EcSakai.eg","Gg.eg","Hs.eg","Mm.eg","Mmu.eg","Pf.plasmo","Pt.eg","Rn.eg","Sc.sgd","Ss.eg","Xl.eg"),sep='')
names(.ORG_PACKAGES)=c("anoGam","Arabidopsis","bosTau","ce","canFam","dm","danRer","E. coli K12","E. coli Sakai","galGal","hg","mm","rheMac","Malaria","panTro","rn","sacCer","Pig","xenTro")
#These are the only formats supported by getgo at the moment...
.ID_MAP=c("Entrez Gene ID","Entrez Gene ID","Ensembl gene ID","Gene Symbol")
names(.ID_MAP)=c("knownGene","refGene","ensGene","geneSymbol")
