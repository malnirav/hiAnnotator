\name{doAnnotation}
\alias{doAnnotation}
\title{Annotate a RangedData object using one of annotation functions.}
\usage{
  doAnnotation(annotType = NULL, ...)
}
\arguments{
  \item{annotType}{one of following: within, nearest,
  counts.}

  \item{...}{Additional parameters to be passed to the
  respective annotation function.}
}
\value{
  a RangedData object with new annotation columns appended
  at the end of sites.rd.
}
\description{
  This is a wrapper function which calls one of following
  functions depending on annotType parameter:
  \code{\link{getFeatureCounts}},
  \code{\link{getNearestFeature}},
  code{\link{getSitesInFeature}}
}
\examples{
# Convert a dataframe to RangedData object
data(sites)
head(sites)
alldata.rd <- makeRangedData(sites,soloStart=TRUE)
alldata.rd

data(genes)
head(genes)
genes.rd <- makeRangedData(genes)
genes.rd

library(BSgenome.Hsapiens.UCSC.hg18)

doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene")
doAnnotation(annotType="counts",alldata.rd,genes.rd,"NumOfGene",seqlengths(Hsapiens))
doAnnotation(annotType="nearest",alldata.rd,genes.rd,"NearestGene")
}
\seealso{
  \code{\link{makeRangedData}},
  \code{\link{getFeatureCounts}},
  \code{\link{getNearestFeature}},
  code{\link{getSitesInFeature}}.
}

