\name{makeRangedData}
\alias{makeRangedData}
\title{Make a sorted RangedData object from a dataframe.}
\usage{
  makeRangedData(x, positionsOnly = FALSE,
    soloStart = FALSE, chromCol = NULL, strandCol = NULL,
    startCol = NULL, stopCol = NULL)
}
\arguments{
  \item{x}{data frame to be converted into RangedData
  object}

  \item{positionsOnly}{boolean flag indicating to return
  only position based data or everything from the data
  frame. Defaults to FALSE.}

  \item{soloStart}{flag denoting whether only one position
  based column is available. In other words, only starts
  are present and no stops. Default=FALSE.}

  \item{chromCol}{use the defined column name for
  space/chromosome based data from the data frame. Defaults
  to NULL.}

  \item{strandCol}{use the defined column name for strand
  or orientation from the data frame. Defaults to NULL.}

  \item{startCol}{use the defined column name for start
  coordinate from the data frame. Defaults to NULL.}

  \item{stopCol}{use the defined column name for stop
  coordinate from the data frame. Defaults to NULL and not
  required if soloStart=TRUE.}
}
\value{
  a RangedData object converted from x.
}
\description{
  The function converts a dataframe into a RangedData
  object without too much hassle of renaming column names.
  The function finds column names that sound like space,
  chromosome, start, stop, position, etc and puts them in
  respective slots to facilitate the conversion of a
  dataframe to RangedData object. If more than one column
  that sounds like start, stop, or position is present, the
  function will use the first match as the representative.
  It is recommended to run this function before utilizing
  any other annotation functions since it will sort the
  object by chromosome and position for copying annotations
  back to their respective rows accurately.
}
\examples{
# Convert a data frame to RangedData object
data(sites)
head(sites)

makeRangedData(sites,soloStart=TRUE)
makeRangedData(sites,soloStart=TRUE,positionsOnly=TRUE)
makeRangedData(sites)

data(genes)
head(genes)

makeRangedData(genes,soloStart=TRUE)
makeRangedData(genes)
}
\seealso{
  \code{\link{getNearestFeature}},
  \code{\link{getFeatureCounts}},
  \code{\link{getSitesInFeature}}.
}

