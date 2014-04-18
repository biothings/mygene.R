\name{MyGene-class}
\Rdversion{1.1}
\docType{class}
\alias{MyGene-class}
\alias{MyGene}

\title{Class \code{"MyGene"}}
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
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
  \item{\code{getgene(geneid, fields, ...)}:}{ Return the gene object for the given geneid }
  \item{\code{getgenes(geneids, fields, ...)}:}{ Return the list of gene object for the given list of geneids. }
  \item{\code{query(q, ...)}:}{ Return  the query result. }
  \item{\code{querymany(qterms, scopes, ...)}:}{ Return the batch query result. }
  \item{\code{metadata()}:}{ Get metadata for MyGene.info services. }
    
  }
}
