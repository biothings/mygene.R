\name{getGenes}
\docType{methods}
\alias{getGenes}
\alias{getGenes,MyGene-method}
\alias{getGenes,missing-method}

\title{
Return the list of gene objects for the given list of geneids.
}
\description{
This is a wrapper for POST query of "/gene" service.
}
\usage{
getGenes(geneids, fields = c("symbol", "name", "taxid", "entrezgene"), ..., 
          return.as=c("DataFrame", "records", "text"), mygene)
}

\arguments{
  \item{geneids}{
A vector, list, or comm-sep string entrez/ensembl gene ids
}
  \item{fields}{
A vector of fields to return. If fields=="all", all available fields are returned.
}
\item{...}{
Includes \code{species} as well as several other fields. 
View available fields by calling \code{?metadata}.
Also, see http://docs.mygene.info/en/latest/doc/annotation_service.html for complete argument details and syntax.
}
\item{return.as}{
"DataFrame" (default), "records" (list), "text" (JSON).
}
\item{mygene}{
A MyGene object that describes how to connect to data resources. See \code{\link{MyGene-class}}. If missing, default object will be used that accesses the main MyGene.info portal. Default is recommended.
}
}

\value{returns a gene object (DataFrame, list, or JSON text) containing the queried annotations}

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