\name{makeTxDbFromMyGene}
\docType{methods}
\alias{makeTxDbFromMyGene}

\title{
makeTxDbFromMyGene
}
\description{
makeTxDbFromMyGene allows the user to make a TxDb object in the Genomic Features package from a mygene "exons" query using a default MyGene object. 
}
\usage{
makeTxDbFromMyGene(gene.list, scopes, species, returnall=FALSE)
}
\arguments{
  \item{gene.list}{
A list, vector, or comma-separated string of query terms.
}
  \item{scopes}{
Type of types of identifiers, either a list, vector, or a comma-separated fields to specify type of
input qterms, e.g. c("reporter", "ensembl.gene", "symbol")
refer to "http://mygene.info/doc/query_service.html#available_fields" for full list
of fields.
}
  \item{species}{
Names or taxonomy ids
}
  \item{returnall}{
Logical, if TRUE, return list of genes without exons annotations. False by Default.
}}
\value{returns TxDb object}
\details{
makeTxDbFromMyGene invokes either the query or queryMany method and passes the response to construct a \code{TxDb} object. See \code{?TxDb} for utilizing transcript annotations.
}
\seealso{
\code{\link[GenomicFeatures]{makeTxDb}}
}
\examples{
xli <- c('DDX26B','CCDC83','MAST3', 'RPL11')
txdb <- makeTxDbFromMyGene(xli, scopes="symbol", species="human")

}
