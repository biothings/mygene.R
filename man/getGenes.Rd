\name{getGenes}
\docType{methods}
\alias{getGenes}

\title{
Return the list of gene objects for the given list of geneids.
}
\description{
This is a wrapper for POST query of "/gene" service.
}
\usage{
getGenes(geneids, fields = c("symbol", "name", "taxid", "entrezgene", ..., 
          species=c("human", "mouse", "rat"), return.as="DataFrame")
}

\arguments{
  \item{geneids}{
A vector, list, or comm-sep string entrez/ensembl gene ids
}
  \item{fields}{
A vector of fields to return. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass species names or taxonomy ids
}
\item{return.as}{
"DataFrame" (Iranges, default), "records" (list), "text" (JSON).
}
\item{...}{
See http://docs.mygene.info/en/latest/doc/annotation_service.html for complete argument details and syntax.
}}

\references{
Ref: http://docs.mygene.info/en/latest/doc/annotation_service.html
}

\seealso{
\code{\link{getGene}}  
\code{\link{query}} 
\code{\link{queryMany}}
}
\examples{
## Return the list of gene object for the given list of gene ids.
getGenes(c(1017,1018))

## mix types of gene ids
getGenes(c(1017,1018,"ENSG00000148795"))

}