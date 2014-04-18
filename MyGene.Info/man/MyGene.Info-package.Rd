\name{MyGene.Info-package}
\alias{MyGene.Info-package}
\alias{MyGene.Info}
\docType{package}
\title{
Access MyGene.Info services

}
\description{
MyGene.Info_ provides simple-to-use REST web services to query/retrieve gene annotation data. It's designed with simplicity and performance emphasized. *mygene*, is an easy-to-use R wrapper to access MyGene.Info_ services.
}
\details{
\tabular{ll}{
Package: \tab MyGene.Info\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-04-18\cr
License: \tab BSD\cr
Depends: \tab R(>=3.1) httr Rcurl jsonlite rjson\cr
}
An overview of how to use the package, including the most important functions ~~
}
\author{
Adam Mark, Chunlei Wu

Maintainer: <help@mygene.info>
Chunlei Wu 
}
\references{
Wu C, MacLeod I, Su AI (2013) BioGPS and MyGene.info: organizing online, gene-centric information. Nucl. Acids Res. 41(D1): D561-D565. 
http://mygene.info/doc/annotation_service.html
http://mygene.info/doc/query_service.html
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{

~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

mg.getgene(1017)
mg.getgene(1017, 'name,symbol,refseq')
mg.getgenes(list(1017,1018,'ENSG00000148795'))
mg.query('cdk2', size=5)
mg.query('symbol:cdk2', species='human')
mg.querymany(list(1017, '695'), scopes='entrezgene', species='human')
mg.querymany(list(1017, '695'), scopes='entrezgene', species=9606)
}
