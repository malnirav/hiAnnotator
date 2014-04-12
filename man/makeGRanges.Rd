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

  \item{...}{parameters for \code{\link{makeRangedData}}
  except for 'asGRanges'.}
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
dataframe to a GRanges object. If more than one column that
sounds like start, stop, or position is present, the
function will use the first match as the representative. It
is recommended to run this function before utilizing any
other annotation functions since it will sort the object by
chromosome and position for copying annotations back to
their respective rows confidently. This function wraps
around \code{\link{makeRangedData}} with asGRanges
parameter enabled. It is recommended to use a
\code{\link[GenomicRanges]{GRanges}} object instead of a
\code{\link[IRanges]{RangedData}} object if number of
distinct chromosomes or targets is greater than 50.
}
\examples{
# Convert a dataframe to GRanges object
data(genes)

makeGRanges(genes, soloStart=TRUE)
makeGRanges(genes)
#makeGRanges(genes, freeze="hg18", soloStart=TRUE)
#makeGRanges(genes, freeze="hg18")
}
\seealso{
\code{\link{makeRangedData}},
\code{\link{getNearestFeature}},
\code{\link{getFeatureCounts}},
\code{\link{getSitesInFeature}}.
}

