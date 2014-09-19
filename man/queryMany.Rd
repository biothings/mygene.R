\name{queryMany}
\docType{methods}
\alias{queryMany}


\title{
Return the batch query result.
}
\description{
This is a wrapper for POST query of "/query" service.
}
\usage{
queryMany(qterms, scopes=NULL, fields=c("name", "symbol", "taxid", "entrezgene"), 
          species=c("human","mouse","rat"), ..., skip=0, entrezonly=FALSE, 
          return.as=c("DataFrame", "records", "text"), returnall=FALSE)
}

\arguments{
  \item{qterms}{
A vector or list, or string of comma-separated query terms
}
  \item{fields}{
a Vector of fields to return. If fields=="all", all available fields are returned.
}
  \item{scopes}{
Type of types of identifiers, either a list or a comma-separated fields to specify type of
input qterms, e.g. c("reporter", "ensembl.gene", "symbol")
refer to "http://mygene.info/doc/query_service.html#available_fields" for full list
of fields.
}
  \item{species}{
Optionally, you can pass species names or taxonomy ids
}
  \item{skip}{
The number of results to skip. Default: 0.
}
  \item{entrezonly}{
If True, return only matching entrez genes, otherwise, including matching Ensemble-only genes (those have no matching entrez genes).
}
\item{returnall}{
if True, return a list of all related data, including dup. and missing qterms. FALSE by default.
}
\item{return.as}{
"DataFrame" (IRanges, default), "records" (list), "text" (JSON).
}
\item{...}{
See http://docs.mygene.info/en/latest/doc/query_service.html for complete argument details and syntax.
}}
\references{
Ref: http://docs.mygene.info/en/latest/doc/query_service.html
}

\seealso{
\code{\link{query}}
\code{\link{getGene}}
\code{\link{getGenes}}
}
\examples{
## return the batch query result
queryMany(c('1053_at', '117_at', '121_at'), scopes="reporter", fields="ensembl.gene", 
          species="human", return.as="records")

queryMany(c('1053_at', '117_at', '121_at'), scopes="reporter", species=9606)

queryMany(c('DDX26B', 'CCDC83', 'MAST3', 'FLOT1'), scopes="symbol", fields="entrezgene", species="human")
}

