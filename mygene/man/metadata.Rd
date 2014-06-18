\name{metadata}

\title{
metadata
}
\description{
Get metadata for MyGene.info services.
}
\usage{
metadata(mygene)
}

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

