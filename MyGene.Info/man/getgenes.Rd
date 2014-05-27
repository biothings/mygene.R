\name{getgenes}
\alias{getgenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the list of gene object for the given list of geneids.
}
\description{
This is a wrapper for POST query of "/gene" service.
}
\usage{
getgenes(list(geneids), fields = 'symbol,name,taxid,entrezgene', ...)
}

\arguments{
  \item{geneids}{
A list or comm-sep string entrez/ensembl gene ids
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass comma-separated species names ortaxonomy ids
}}

\references{
Ref: http://mygene.info/doc/annotation_service.html
}

\seealso{
\code{\link{getgene, query, querymany}}
}
\examples{
mg.getgenes(list(1017,1018))
mg.getgenes(list(1017,1018,'ENSG00000148795'), species='human')

}