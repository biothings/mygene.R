\name{getGene}
\docType{methods}
\alias{getGene}


\title{
Return the gene object for the given geneid.
}
\description{
This is a wrapper for GET query of "/gene/<geneid>" service.
}
\usage{
getGene(geneid, fields = c("symbol", "name", "taxid", "entrezgene", ..., return.as)
}
\arguments{
  \item{geneid}{
Entrez/ensembl gene id
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass comma-separated species names or taxonomy ids
}
\item{return.as}{
"records"" (list), "text" (JSON).
}}
\references{
http://mygene.info/doc/annotation_service.html
}

\seealso{
\code{\link{getGenes}} 
\code{\link{query}} 
\code{\link{queryMany}}
}
\examples{
## return the gene object for the given gene id
getGene(1017)

## customize fields
getGene(1017, fields=c("name", "symbol", "refseq"), return.as="text")

## all fields
getGene(1017, fields="all")
}

