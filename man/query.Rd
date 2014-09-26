\name{query}
\docType{methods}
\alias{query}

\title{
Return  the query result.
}
\description{
This is a wrapper for GET query of "/query?q=<query>" service.
}
\usage{
query(q, fields=c("name", "symbol", "taxid", "entrezgene"), 
          species=c("human","mouse","rat"), size=10, skip=0, ..., entrezonly=FALSE,
          return.as=c("DataFrame", "records", "text"))
}
\arguments{
  \item{q}{
query term(s), see query syntax at http://mygene.info/doc/query_service.html#query-syntax
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all
available fields are returned.
}
  \item{species}{
Optionally, you can pass comma-separated species names or taxonomy ids
}
  \item{size}{
The maximum number of results to return (with a cap of 1000 at the
moment). Default: 10.
}
  \item{skip}{
The number of results to skip. Default: 0.
}
  \item{entrezonly}{
If True, return only matching entrez genes, otherwise, including matching 
Ensemble-only genes (those have no matching entrez genes).
}
\item{return.as}{
"DataFrame" (IRanges, default), "records" (list), or "text" (JSON).
}
\item{...}{
See http://docs.mygene.info/en/latest/doc/query_service.html for complete argument details and syntax.
}
}
\value{returns a gene object (list, data.frame, or JSON text) containing the queried annotations}
\references{
Ref: http://docs.mygene.info/en/latest/doc/query_service.html
}

\seealso{
\code{\link{queryMany}} 
\code{\link{getGene}} 
\code{\link{getGenes}}
}
\examples{
## return the query result
query("cdk2", size=5)

query("reporter:1000_at")

query("symbol:cdk2", return.as="text")

query(q="cyclin-dependent kinase", fields="uniprot")
}
