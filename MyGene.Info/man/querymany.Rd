\name{querymany}
\alias{querymany}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the batch query result.
}
\description{
This is a wrapper for POST query of "/query" service.
}
\usage{
querymany(qterms, scopes=NULL, ...)
}

\arguments{
  \item{qterms}{
A list of query terms, or a string of comma-separated query terms
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
  \item{species}{
Optionally, you can pass comma-separated species names ortaxonomy ids
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
}}
\references{
http://mygene.info/doc/query_service.html
}

\seealso{
\code{\link{query, getgene, getgenes}}
}
\examples{
mg.querymany(list(1017, '695'), scopes='entrezgene', species='human')
mg.querymany(list(1017, '695'), scopes='entrezgene', species=9606)
mg.querymany('1017,695', scopes='entrezgene', species='human')
}

