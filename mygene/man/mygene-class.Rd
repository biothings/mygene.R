\name{mygene-class}
\Rdversion{1.1}
\docType{class}
\alias{MyGene-class}


\title{Class \code{"mygene"}}
\description{
R Client to access MyGene.Info services
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
Wu C, MacLeod I, Su AI (2013) BioGPS and MyGene.info: organizing online, gene-centric information. Nucl. Acids Res. 41(D1): D561-D565.
}
\author{
Adam Mark, Chunlei Wu
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{Url}:}{http://mygene.info/v2 }
    \item{\code{delay}:}{sleep time }
    \item{\code{step}:}{batch limit }

  }
}
\section{Methods}{
  \describe{
  \item{\code{mg.getgene(geneid, fields, ...)}:}{ Return the gene object for the given geneid }
  \item{\code{mg.getgenes(geneids, fields, ...)}:}{ Return the list of gene object for the given list of geneids. }
  \item{\code{mg.query(q, ...)}:}{ Return  the query result. }
  \item{\code{mg.querymany(qterms, scopes, ...)}:}{ Return the batch query result. }
  \item{\code{mg.metadata()}:}{ Get metadata for MyGene.info services. }
    
  }
}
