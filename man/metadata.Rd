\name{metadata}
\docType{methods}
\alias{metadata}
\alias{metadata,MyGene-method}

\title{
metadata
}
\description{
Get metadata for MyGene.info services.
}
\usage{
metadata(x, ...)
}
\arguments{
  \item{x}{
MyGene object
}
  \item{...}{
MyGene object slot parameters
}}

\value{returns the metadata including available fields, genome assemblies, sources, statistics, taxonomy, and timestamp}

\references{
http://mygene.info/v2/metadata
}

\examples{
## Get metadata 
mygene<-MyGene()
metadata(mygene)

## get available fields to return
metadata(mygene)$available_fields
}
