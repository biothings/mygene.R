\name{queryMany-methods}
\docType{methods}
\alias{queryMany-method}


\title{
Return the batch query result.
}
\description{
This is a wrapper for POST query of "/query" service.
}
\usage{
queryMany(qterms, scopes=NULL, ..., return.as)
}

\arguments{
  \item{qterms}{
A vector or list, or string of comma-separated query terms
}
  \item{fields}{
a Vector of fields to return. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass species names ortaxonomy ids
}
  \item{size}{
The maximum number of results to return (with a cap of 1000 at the moment). Default: 10.
}
  \item{skip}{
The number of results to skip. Default: 0.
}
  \item{sort}{
Prefix with "-" for descending order, otherwise in ascending order.Default: sort by matching scores in decending order.
}
  \item{entrezonly}{
If True, return only matching entrez genes, otherwise, including matching Ensemble-only genes (those have no matching entrez genes).
}
\item{returnall}{
if True, return a list of all related data, including dup. and missing qterms.
}
\item{verbose}{
If True (default), print out infomation about dup and missing qterms.
}
\item{return.as}{
"data.frame", "records" (list), "text" (JSON), otherwise return as "DataFrame" (IRanges).
}}
\references{
http://mygene.info/doc/query_service.html
}

\seealso{
\code{\link{query}}
\code{\link{getGene}}
\code{\link{getGenes}}
}
\examples{
## return the batch query result
queryMany(c(1018, 695), scopes="entrezgene", fields="ensembl.gene", species="human", return.as="records")

queryMany(c(1017, '695'), scopes='entrezgene', species=9606)

queryMany('1017,695', scopes='entrezgene', species='human')
}

