\name{MyGene-class}
\Rdversion{1.1}
\docType{class}
\alias{MyGene}

\title{Class \code{"MyGene"}}
\description{
R Client to access MyGene.Info annotation services
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{MyGene(...)}.

}
\section{Slots}{
  \describe{
    \item{\code{base.url}:}{"http://mygene.info/v2". Object of class \code{"character"}}
    \item{\code{delay}:}{Sleep time between batch retrieval. Object of class \code{"numeric"}}
    \item{\code{step}:}{Batch limit. Object of class \code{"numeric"}}
    \item{\code{version}:}{httr package version. Object of class \code{"character"}}
    \item{\code{verbose}:}{Object of class \code{"logical"}}
    \item{\code{debug}:}{Object of class \code{"logical"}}
  }
}
\section{Methods}{
  \describe{
  \item{\code{getGene(geneid, fields, ...)}:}{ Return the gene object for the given geneid }
  \item{\code{getGenes(geneids, fields, ...)}:}{ Return the list of gene object for the given list of geneids. }
  \item{\code{query(q, ...)}:}{ Return  the query result. }
  \item{\code{queryMany(qterms, scopes, ...)}:}{ Return the batch query result. }
  \item{\code{metaData()}:}{ Get metadata for MyGene.info services. }
    
  }
}
\references{
Wu C, MacLeod I, Su AI (2013) BioGPS and MyGene.info: organizing online, gene-centric information. Nucl. Acids Res. 41(D1): D561-D565.
}
\author{
Adam Mark, Chunlei Wu, Ryan Thompson
}

\examples{
showClass("MyGene")
}
\keyword{classes}
\usage{
mygene<-MyGene(base.url="http://mygene.info/v2", delay=1, step=1000,
        version=version, verbose=TRUE, debug=FALSE)
}

