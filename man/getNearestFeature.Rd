\name{getNearestFeature}
\alias{getNearestFeature}
\title{Get nearest annotation boundary for a position range.}
\usage{
  getNearestFeature(sites.rd, features.rd, colnam = NULL,
    side = "either", feature.colnam = NULL,
    strand.colnam = NULL, parallel = FALSE)
}
\arguments{
  \item{sites.rd}{RangedData/GRanges object to be used as
  the query.}

  \item{features.rd}{RangedData/GRanges object to be used
  as the subject or the annotation table.}

  \item{colnam}{column name to be added to sites.rd for the
  newly calculated annotation...serves a core!}

  \item{side}{boundary of annotation to use to calculate
  the nearest distance. Options are '5p','3p', or the
  default 'either'.}

  \item{feature.colnam}{column name from features.rd to be
  used for retrieving the nearest feature name. By default
  this is NULL assuming that features.rd has a column that
  includes the word 'name' somewhere in it.}

  \item{strand.colnam}{column name from features.rd to be
  used for retrieving the nearest feature's orientation. By
  default this is NULL assuming that features.rd has a
  column that includes the word 'strand' somewhere in it.
  If it doesn't the function will assume the supplied
  annotation is in '+' orientation (5' -> 3'). The same
  applies to strand column in sites.rd.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to
  FALSE. If no parallel backend is registered, then a
  serial version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}
}
\value{
  a RangedData/GRanges object with new annotation columns
  appended at the end of sites.rd.
}
\description{
  Given a query object, the function retrieves the nearest
  feature and its properties from a subject and then
  appends them as new columns within the query object. When
  used in genomic context, the function can be used to
  retrieve a nearest gene 5' or 3' end relative to genomic
  position of interest.
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
# Convert a dataframe to RangedData/GRanges object
data(sites)
head(sites)
alldata.rd <- makeRangedData(sites,soloStart=TRUE)
alldata.rd

data(genes)
head(genes)
genes.rd <- makeRangedData(genes)
genes.rd

nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene")
nearestGenes
nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="5p")
nearestGenes
nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="3p")
nearestGenes

# Parallel version of getNearestFeature
nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene", parallel=TRUE)
nearestGenes
}
\seealso{
  \code{\link{makeRangedData}}, \code{\link{makeGRanges}},
  \code{\link{getFeatureCounts}},
  \code{\link{getSitesInFeature}},
  \code{\link{get2NearestFeature}}.
}

