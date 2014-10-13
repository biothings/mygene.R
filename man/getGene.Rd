\name{getGene}
\docType{methods}
\alias{getGene}
\alias{getGene,MyGene-method}
\alias{getGene,missing-method}

\title{
Return the gene object for the given geneid.
}
\description{
This is a wrapper for GET query of "/gene/<geneid>" service.
}
\usage{
getGene(geneid, fields = c("symbol", "name", "taxid", "entrezgene"), 
        ..., return.as=c("records", "text"), mygene)
}
\arguments{
  \item{geneid}{
Entrez/ensembl gene id
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
\item{...}{
Includes \code{species} as well as several other fields. 
View available fields by calling \code{?metadata}.
Also, see http://docs.mygene.info/en/latest/doc/annotation_service.html for complete argument details and syntax.
}
\item{return.as}{
"records"" (list), "text" (JSON).
}
\item{mygene}{
A MyGene object that describes how to connect to data resources. See \code{\link{MyGene-class}}. If missing, default object will be used that accesses the main MyGene.info portal. Default is recommended.
}
}
\value{returns a gene object containing the queried annotations}

\references{
Ref: http://docs.mygene.info/en/latest/doc/annotation_service.html
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

