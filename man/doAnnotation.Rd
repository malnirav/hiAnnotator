\name{doAnnotation}
\alias{doAnnotation}
\title{Annotate a RangedData/GRanges object using one of annotation functions.}
\usage{
  doAnnotation(annotType = NULL, ...,
    postProcessFun = NULL, postProcessFunArgs = list())
}
\arguments{
  \item{annotType}{one of following: within, nearest,
  twoNearest, counts, countsBig.}

  \item{...}{Additional parameters to be passed to the
  respective annotation function.}

  \item{postProcessFun}{function to call on the resulting
  object for any post processing steps.}

  \item{postProcessFunArgs}{additional arguments for
  postProcessFun as a list.}
}
\value{
  a RangedData/GRanges object with new annotation columns
  appended at the end of sites.rd.
}
\description{
  This is a wrapper function which calls one of following
  functions depending on annotType parameter:
  \code{\link{getFeatureCounts}},
  \code{\link{getFeatureCountsBig}},
  \code{\link{getNearestFeature}},
  \code{\link{get2NearestFeature}},
  code{\link{getSitesInFeature}}
}
\examples{
# Convert a dataframe to RangedData/GRanges object
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
doAnnotation(annotType="counts",alldata.rd,genes.rd,"NumOfGene")
doAnnotation(annotType="nearest",alldata.rd,genes.rd,"NearestGene")
doAnnotation(annotType="countsBig",alldata.rd,genes.rd,"ChipSeqCounts")
geneCheck <- function(x,wanted) { x$isWantedGene <- x$InGene \%in\% wanted; return(x) }
doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene", postProcessFun=geneCheck, postProcessFunArgs=list("wanted"=c("FOXJ3","SEPT9","RPTOR")) )
}
\seealso{
  \code{\link{makeRangedData}}, \code{\link{makeGRanges}},
  \code{\link{getFeatureCounts}},
  \code{\link{getFeatureCountsBig}},
  \code{\link{getNearestFeature}},
  \code{\link{get2NearestFeature}},
  code{\link{getSitesInFeature}}.
}

