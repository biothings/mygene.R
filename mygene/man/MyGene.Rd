\name{MyGene-class}
\Rdversion{1.1}
\docType{class}
\alias{MyGene-class}


\title{Class \code{"MyGene"}}
\description{
R Client to access MyGene.Info services
}

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
    \item{\code{base.url}:}{http://mygene.info/v2 }
    \item{\code{delay}:}{sleep time }
    \item{\code{step}:}{batch limit }
    \item{\code{verbose}:}{print out information }

  }
}
\section{Methods}{
  \describe{
  \item{\code{getGene(geneid, fields, ...)}:}{ Return the gene object for the given geneid }
  \item{\code{getGenes(geneids, fields, ...)}:}{ Return the list of gene object for the given list of geneids. }
  \item{\code{query(q, ...)}:}{ Return  the query result. }
  \item{\code{queryMany(qterms, scopes, ...)}:}{ Return the batch query result. }
  \item{\code{metadata()}:}{ Get metadata for MyGene.info services. }
    
  }
}
