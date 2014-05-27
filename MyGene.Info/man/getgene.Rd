\name{getgene}
\alias{getgene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the gene object for the given geneid.
}
\description{
This is a wrapper for GET query of "/gene/<geneid>" service.
}
\usage{
getgene(geneid, fields = 'symbol,name,taxid,entrezgene', ...)
}
\arguments{
  \item{geneid}{
Entrez/ensembl gene id
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass comma-separated species names ortaxonomy ids
}}
\references{
http://mygene.info/doc/annotation_service.html
}

\seealso{
\code{\link{getgenes, query, querymany}}
}
\examples{
mg.getgene(1017)
mg.getgene(1017, 'name,symbol,refseq')
}
