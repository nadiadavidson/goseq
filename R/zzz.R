#  ZZZ.R

.onAttach <- function(libname, pkgname)
{
	if(.Platform$OS.type=="windows" && .Platform$GUI=="Rgui" ) {
		winMenuAddItem("Vignettes","goseq","shell.exec(system.file(\"doc\",\"goseq.pdf\",package=\"goseq\"))")
	}
}

.ORG_PACKAGES=paste("org.",c("Ag.eg","At.tair","Bt.eg","Ce.eg","Cf.eg","Dm.eg","Dr.eg","EcK12.eg","EcSakai.eg","Gg.eg","Hs.eg","Mm.eg","Mmu.eg","Pf.plasmo","Pt.eg","Rn.eg","Sc.sgd","Ss.eg","Xl.eg"),sep='')
names(.ORG_PACKAGES)=c("A. gambiae","Arabidopsis","Cow","C. elegans","Dog","D. melanogaster","Zebrafish","E. coli K12","E. coli Sakai","Chicken","Human","Mouse","Rhesus","Malaria","Chimp","Rat","Yeast","Pig","X. tropicalis")
