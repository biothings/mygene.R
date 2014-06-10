\name{mg.getgenes}
\alias{mg.getgenes}

\title{
Return the list of gene objects for the given list of geneids.
}
\description{
This is a wrapper for POST query of "/gene" service.
}
\usage{
mg.getgenes(geneids, fields = 'symbol,name,taxid,entrezgene', ...)
}

\arguments{
  \item{geneids}{
A list, vector, or comm-sep string entrez/ensembl gene ids
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass comma-separated species names ortaxonomy ids
}
\item{return.as.list}{
If True, return gene objects as a list, otherwise return as DataFrame.
}}

\references{
Ref: http://mygene.info/doc/annotation_service.html
}

\seealso{
\code{\link{mg.getgene}}  
\code{\link{mg.query}} 
\code{\link{mg.querymany}}
}
\examples{
## Return the list of gene object for the given list of gene ids.
mg.getgenes(c(1017,1018))

## mix types of gene ids
mg.getgenes(c(1017,1018,'ENSG00000148795'), species='human')

}