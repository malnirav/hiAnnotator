\name{makeChunks}
\alias{makeChunks}
\title{Breaks two GRanges/RangedData objects into chunks of N size.}
\usage{
  makeChunks(sites.rd, features.rd, chunkSize = NULL)
}
\arguments{
  \item{sites.rd}{a RangedData/GRanges object.}

  \item{features.rd}{a RangedData/GRanges object.}

  \item{chunkSize}{number of rows to use per chunk of
  query. Default to length(sites.rd)/detectCores() or
  length(query)/getDoParWorkers() depending on parallel
  backend registered.}
}
\value{
  a list of RangedData/GRanges objects where each element
  is of length 2 representing query & subject chunks.
}
\description{
  Given a query and subject GRanges/RangedData objects, the
  function breaks query into chunks of N size where each
  chunk has a respective subject object filtered by
  seqnames/space present in the query chunk. This is a
  helper function used by one of the annotation function in
  'See Also' section where each chunk is sent to a parallel
  node for processing.
}
\examples{
data(sites)
data(genes)
sites <- makeGRanges(sites,soloStart=TRUE)
genes <- makeGRanges(genes)
makeChunks(sites, genes)
}
\seealso{
  \code{\link{makeGRanges}}, \code{\link{doAnnotation}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}},
  \code{\link{getSitesFeatureCounts}}.
}

