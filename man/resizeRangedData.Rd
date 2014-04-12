\name{resizeRangedData}
\alias{resizeRangedData}
\title{Resize a RangedData object.}
\usage{
resizeRangedData(rd, width = NULL, boundary = "center", spaceSizes = NULL,
  spaceMin = 1, limitLess = FALSE)
}
\arguments{
  \item{rd}{a RangedData object}

  \item{width}{the width of the resized ranges.}

  \item{boundary}{same as fix parameter in
  \code{\link[IRanges]{resize}}. One of "start", "end", and
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
truncated to spaceSizes and spaceMin unless limitLess=TRUE.
}
\description{
Function to resize a RangedData object by the given width,
max space/chromosome size, and the boundary. This is one of
the helper function used in
\code{\link{getFeatureCountsBig}}.
}
\examples{
# Convert a dataframe to RangedData object
data(sites)
alldata.rd <- makeGRanges(sites,soloStart=TRUE)
seqLengths <- structure(c(247249719L, 242951149L, 199501827L, 191273063L, 180857866L,170899992L, 158821424L, 146274826L, 140273252L, 135374737L, 134452384L, 132349534L, 114142980L, 106368585L, 100338915L, 88827254L, 78774742L, 76117153L, 63811651L, 62435964L, 46944323L, 49691432L, 154913754L, 57772954L, 16571L), .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"))
resizeRangedData(alldata.rd,width=10000,spaceSizes=seqLengths)
resizeRangedData(alldata.rd,width=10000,limitLess=TRUE)
}
\seealso{
\code{\link{getFeatureCountsBig}},
\code{\link{makeRangedData}}, \code{\link{makeGRanges}},
\code{\link{getNearestFeature}},
\code{\link{getSitesInFeature}}.
}

