\name{getFeatureCounts}
\alias{getFeatureCounts}
\title{Get counts of annotation within a defined window around each query range/position.}
\usage{
  getFeatureCounts(sites.rd, features.rd, colnam = NULL,
    chromSizes = NULL,
    widths = c(`1k` = 500, `10k` = 5000, `1M` = 5e+05),
    weightsColname = NULL, doInChunks = FALSE,
    chunkSize = 10000, parallel = FALSE, ...)
}
\arguments{
  \item{sites.rd}{RangedData object to be used as the
  query.}

  \item{features.rd}{RangedData obj to be used as the
  subject or the annotation table.}

  \item{colnam}{column name to be added to sites.rd for the
  newly calculated annotation...serves as a prefix to
  windows sizes!}

  \item{chromSizes}{named vector of chromosome/space sizes
  to be used for testing if a position is off the mappable
  region.}

  \item{widths}{named vector of window sizes to be used for
  casting a net around each position. Default: \code{c("1k"
  = 500,"10k" = 5000,"1M" = 500000)}}

  \item{weightsColname}{if defined, weigh each row from
  features.rd when performing the counts}

  \item{doInChunks}{break up sites.rd into small pieces of
  chunkSize to perform the calculations. Default to FALSE.
  Useful if you are expecting to find great deal of overlap
  between sites.rd and features.rd.}

  \item{chunkSize}{number of rows to use per chunk of
  sites.rd. Default to 10000. Only used if
  doInChunks=TRUE.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to
  FALSE. If no parallel backend is registered, then a
  serial version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}

  \item{...}{Additional parameters for
  \code{\link{findOverlaps}}.}
}
\value{
  a RangedData object with new annotation columns appended
  at the end of sites.rd.
}
\description{
  Given a query object and window size(s), the function
  first augments the ranges within the query by flanking
  starts and stops with window width. Therefore, a start of
  12 and end of 14 with width 10 will yield a range of
  2,24. This new range is then compared against the subject
  to find any overlapping ranges and then tallied up. If
  weights are assigned to each positions in the subject,
  then tallied counts are multiplied accordingly. If
  annotation object is large, spanning greater than 100
  million rows, then getFeatureCountsBig is used which
  drops any weight column if specified or additional
  parameters passed to \code{\link{findOverlaps}}.
}
\note{
  Try not to use this function for >50 spaces unless you
  have tons fo memory. If parallel=TRUE, then be sure to
  have a paralle backend registered before running the
  function. One can use any of the following libraries
  compatible with \code{\link{foreach}}: doMC, doSMP,
  doSNOW, doMPI. For example: library(doSMP); w <-
  startWorkers(2); registerDoSMP(w)
}
\examples{
# Convert a data frame to RangedData object
data(sites)
head(sites)
alldata.rd <- makeRangedData(sites,soloStart=TRUE)
alldata.rd

data(genes)
head(genes)
genes.rd <- makeRangedData(genes)
genes.rd

library(BSgenome.Hsapiens.UCSC.hg18)

geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",seqlengths(Hsapiens))
geneCounts
# Parallel version of getFeatureCounts
geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",seqlengths(Hsapiens), parallel=T)
geneCounts
# For large annotations, use getFeatureCountsBig
}
\seealso{
  \code{\link{makeRangedData}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}}.
}

