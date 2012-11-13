\name{get2NearestFeature}
\alias{get2NearestFeature}
\title{Get two nearest upstream and downstream annotation boundary for a position range.}
\usage{
  get2NearestFeature(sites.rd, features.rd, colnam = NULL,
    side = "either", feature.colnam = NULL,
    strand.colnam = NULL)
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
}
\value{
  a RangedData/GRanges object with new annotation columns
  appended at the end of sites.rd.
}
\description{
  Given a query object, the function retrieves the two
  nearest feature upstream and downstream along with their
  properties from a subject and then appends them as new
  columns within the query object. When used in genomic
  context, the function can be used to retrieve two nearest
  gene upstream and downstream of the genomic position of
  interest.
}
\note{
  For cases where a position is at the edge and there are
  no feature up/down stream since it would fall off the
  chromosome, the function simply returns the nearest
  feature. In addition, if there are multiple locations
  where a query falls into, the function arbitrarily
  chooses one to serve as the nearest feature, then reports
  2 upstream & downstream feature. That may occasionally
  yield features which are the same upstream and dowstream,
  which is commonly encountered when studying spliced genes
  or phenomena related to it.
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

nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene")
nearestGenes
nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",side="5p")
nearestGenes
nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",side="3p")
nearestGenes
}
\seealso{
  \code{\link{getNearestFeature}},
  \code{\link{makeRangedData}}, \code{\link{makeGRanges}},
  \code{\link{getFeatureCounts}},
  \code{\link{getSitesInFeature}}.
}

