\name{makeGRanges}
\alias{makeGRanges}
\title{Make a sorted GRanges object from a dataframe.}
\usage{
  makeGRanges(x, freeze = NULL, ...)
}
\arguments{
  \item{x}{dataframe to be converted into a GRanges object}

  \item{freeze}{UCSC genome version of the data in x.
  Default is NULL. This parameter is generally used to
  populate \code{\link{seqinfo}} slot of GRanges objects.}

  \item{...}{parameters for \code{\link{makeRangedData}}.}
}
\value{
  a GRanges object converted from x.
}
\description{
  The function converts a dataframe into a GRanges object
  without too much hassle of renaming column names. The
  function finds column names that sound like seqname,
  chromosome, start, stop, position, etc and puts them in
  respective slots to facilitate the conversion of a
  dataframe to GRanges object. If more than one column that
  sounds like start, stop, or position is present, the
  function will use the first match as the representative.
  It is recommended to run this function before utilizing
  any other annotation functions since it will sort the
  object by chromosome and position for copying annotations
  back to their respective rows accurately. This function
  wraps around \code{\link{makeRangedData}}.
}
\examples{
# Convert a dataframe to GRanges object
data(genes)

makeGRanges(genes,soloStart=TRUE)
makeGRanges(genes, freeze="hg18", soloStart=TRUE)
makeGRanges(genes)
makeGRanges(genes, freeze="hg18")
}
\seealso{
  \code{\link{makeRangedData}},
  \code{\link{getNearestFeature}},
  \code{\link{getFeatureCounts}},
  \code{\link{getSitesInFeature}}.
}

