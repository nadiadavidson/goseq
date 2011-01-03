\name{getgo}
\Rdversion{1.1}
\alias{getgo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fetch GO categories
}
\description{
Obtains all gene ontology (GO) categories associated with a set of genes using the relevant organism package.
}
\usage{
getgo(genes, genome, id,fetch.cats=c("GO:CC","GO:BP","GO:MF"))
}
\arguments{
  \item{genes}{
A vector or list of genes to get the associated GO categories.
}
  \item{genome}{
A string identifying the genome that \code{genes} refer to.  For a list of supported organisms run \code{\link{supportedGenomes}}.

}
  \item{id}{
A string identifying the gene identifier used by \code{genes}.  For a list of supported gene IDs run \code{\link{supportedGeneIDs}}.
}
  \item{fetch.cats}{
A vector specifying which categories to fetch the mapping between category names and genes for.  See details for vaild options.
}
}
\details{
This function attempts to make use of the organism packages (org.<Genome>.<GeneID>.db) to obtain the mapping between gene ID and GO categories.  As with \code{\link{getlength}} it is preferable that the same gene identifier system is used for both summarization and retrieving GO categories.

Valid options for the \code{fetch.cats} arguement are any combination of "GO:CC", "GO:BP", "GO:MF" & "KEGG".  The three GO terms refer to the Cellular Component, Biological Process and Molecular Function respectively.  "KEGG" refers to KEGG pathways.

Note that \code{getgo} is a convenience function, designed to make extracting mappings between GO categories and Gene ID easy.  For less common organisms and/or gene ID \code{getgo} may fail to return a mapping even when a legitimate mapping exists in the revelant organism package.  If \code{getgo} fails, you should always try to build the mapping yourself from the organism package (if one exists) before deciding that the information is unavailable.  Further information and examples of this can be found in the package Vignette.
}
\value{
A list where each entry is named by a gene and contains a vector of all the associated GO categories.  This can be used directly with the \code{gene2cat} option in \code{\link{goseq}}.
}
\author{
Matthew D. Young \email{myoung@wehi.edu.au}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{supportedGenomes}}, \code{\link{supportedGeneIDs}}, \code{\link{goseq}}
}
\examples{
genes <- c("ENSG00000124208", "ENSG00000182463", "ENSG00000124201", "ENSG00000124205", "ENSG00000124207")
getgo(genes,'hg19','ensGene')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
