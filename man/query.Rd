\name{query}
\docType{methods}
\alias{query}
\alias{query,MyGene-method}
\alias{query,missing-method}

\title{
Return  the query result.
}
\description{
This is a wrapper for GET query of "/query?q=<query>" service.
}
\usage{
query(q, ..., return.as=c("DataFrame", "records", "text"), mygene)
}
\arguments{
  \item{q}{
query term(s), see query syntax at http://mygene.info/doc/query_service.html#query-syntax
}
  \item{...}{
Commonly queried fields include \code{species, fields, size} as well as several other fields. 
View available fields by calling \code{?metadata}.
Also, see http://docs.mygene.info/en/latest/doc/query_service.html for complete argument details and syntax.
}
  \item{return.as}{
"DataFrame" (default), "records" (list), or "text" (JSON).
}
  \item{mygene}{
A MyGene object that describes how to connect to data resources. See \code{\link{MyGene-class}}. If missing, default object will be used that accesses the main MyGene.info portal. Default is recommended.
}
}
\value{returns a gene object (DataFrame, list, or JSON text) containing the queried annotations}
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
