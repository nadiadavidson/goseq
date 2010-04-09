\name{nullp}
\Rdversion{1.1}
\alias{nullp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability Weighting Function
}
\description{
Calculates a Probability Weighting Function for a set of genes based on a given set of biased data (usually gene length) and each genes status as differentially expressed or not.
}
\usage{
nullp(DEgenes, genome, id, bias.data=NULL,plot.fit=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DEgenes}{
A named binary vector where 1 represents DE, 0 not DE and the names are gene IDs.
}
  \item{genome}{
A string identifying the genome that \code{genes} refer to.  For a list of supported organisms run \code{\link{supportedGenomes}}.

}
  \item{id}{
A string identifying the gene identifier used by \code{genes}.  For a list of supported gene IDs run \code{\link{supportedGeneIDs}}.
}
  \item{bias.data}{
A numeric vector containing the data on which the DE may depend.  Usually this is the median transcript length of each gene in bp.  If set to \code{NULL} \code{nullp} will attempt to fetch length using \code{\link{getlength}}.
}
  \item{plot.fit}{
Plot the PWF or not?
}
}
\details{
It is essential that the entire analysis pipeline, from summarizing raw reads through to using \code{goseq} be done in just one gene identifier format.  If your data is in a different format you will need to obtain the gene lengths and supply them to the \code{nullp} function using the \code{bias.data} arguement.  Converting to a supported format from another format should be avoided whenever possible as this will almost always result in data loss.

\code{NA}s are allowed in the bias.data vector if you do not have information about a certain gene.  Setting a gene to \code{NA} is preferable to removing it from the analysis.

If \code{bias.data} is left as \code{NULL}, \code{nullp} attempts to use \code{\link{getlength}} to fetch GO catgeory to gene identifier mappings.

It is recommended you review the fit produced by the \code{nullp} function before proceeding by leaving \code{plot.fit} as \code{TRUE}.
}
\value{
A numeric vector containing the value on the probability weighting function for each gene.  This is usually passed to the function \code{goseq} via the \code{pwf} arguement.
}
\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010) \emph{Gene ontology analysis for RNA-seq: accounting for selection bias}
  Genome Biology
  Date: Feb 2010
  Vol: 11
  Issue: 2
  Pages: R14
}
\author{
Matthew D. Young \email{myoung@wehi.edu.au}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{supportedGenomes}}, \code{\link{supportedGeneIDs}}, \code{\link{goseq}}, \code{\link{getlength}}
}
\examples{
data(prostate)
pwf <- nullp(genes, 'hg19', 'ensGene')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
