\name{resizeRangedData}
\alias{resizeRangedData}
\title{Resize a RangedData object.}
\usage{
  resizeRangedData(rd, width = NULL, boundary = "center",
    spaceSizes = NULL, spaceMin = 1, limitLess = FALSE)
}
\arguments{
  \item{rd}{a RangedData object}

  \item{width}{the width of the resized ranges.}

  \item{boundary}{same as fix parameter in
  \code{\link{resize}}. One of "start", "end", and
  "center". Default is "center".}

  \item{spaceSizes}{named vector of chromosome/space sizes
  to be used for fixing ranges off the limits.}

  \item{spaceMin}{lowest value to be allowed as the start
  position. Default is 1.}

  \item{limitLess}{whether to limit the resized ranges by
  spaceSizes and spaceMin. Default is FALSE. In genomic
  context, this means starts can be <=0 and ends can be >
  chromosome size depending on the width.}
}
\value{
  a RangedData object with ranges sized to width and/or
  truncated to spaceSizes and spaceMin unless
  limitLess=TRUE.
}
\description{
  Function to resize a RangedData object by the given
  width, max space/chromosome size, and the boundary. This
  is one of the helper function used in
  \code{\link{getFeatureCountsBig}}.
}
\examples{
# Convert a dataframe to RangedData object
data(sites)
alldata.rd <- makeRangedData(sites,soloStart=TRUE)
library(BSgenome.Hsapiens.UCSC.hg18)
resizeRangedData(alldata.rd,width=10000,spaceSizes=seqlengths(Hsapiens))
resizeRangedData(alldata.rd,width=10000,limitLess=T)
}
\seealso{
  \code{\link{getFeatureCountsBig}},
  \code{\link{makeRangedData}}, \code{\link{makeGRanges}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}}.
}

