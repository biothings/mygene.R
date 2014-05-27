\name{mg.getgene}
\alias{mg.getgene}

\title{
Return the gene object for the given geneid.
}
\description{
This is a wrapper for GET query of "/gene/<geneid>" service.
}
\usage{
mg.getgene(geneid, fields = 'symbol,name,taxid,entrezgene', ...)
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
\code{\link{mg.getgenes}} 
\code{\link{mg.query}} 
\code{\link{mg.querymany}}
}
\examples{
## return the gene object for the given gene id
mg.getgene(1017)

## customize fields
mg.getgene(1017, fields='name,symbol,refseq')
}
