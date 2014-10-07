\name{queryMany}
\docType{methods}
\alias{queryMany}
\alias{queryMany,MyGene-method}
\alias{queryMany,missing-method}


\title{
Return the batch query result.
}
\description{
This is a wrapper for POST query of "/query" service.
}
\usage{
queryMany(qterms, scopes=NULL, ..., return.as=c("DataFrame", "records", "text"), mygene)
}

\arguments{
  \item{qterms}{
A vector or list, or string of comma-separated query terms
}
  \item{scopes}{
Type of types of identifiers, either a list or a comma-separated fields to specify type of
input qterms, e.g. c("reporter", "ensembl.gene", "symbol")
refer to "http://mygene.info/doc/query_service.html#available_fields" for full list
of fields.
}
  \item{...}{
Commonly queried fields include \code{species, fields, size} as well as several other fields. 
\code{returnall} returns a list of all related data including duplicated and missing qterms. False by default.
View available fields by calling \code{?metadata}.
Also, see http://docs.mygene.info/en/latest/doc/query_service.html for complete argument details and syntax.
}
  \item{return.as}{
"DataFrame" (IRanges, default), "records" (list), "text" (JSON).
}
  \item{mygene}{
If missing, default object will be used.
}
}
\value{returns a gene object (list, DataFrame, or JSON text) containing the queried annotations}
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

