\name{get2NearestFeature}
\alias{get2NearestFeature}
\title{Get two nearest upstream and downstream annotation boundary for a position range.}
\usage{
get2NearestFeature(sites.rd, features.rd, colnam = NULL, side = "either",
  feature.colnam = NULL, relativeTo = "subject")
}
\arguments{
  \item{sites.rd}{RangedData/GRanges object to be used as
  the query.}

  \item{features.rd}{RangedData/GRanges object to be used
  as the subject or the annotation table.}

  \item{colnam}{column name to be added to sites.rd for the
  newly calculated annotation...serves a core!}

  \item{side}{boundary of annotation to use to calculate
  the nearest distance. Options are '5p','3p',
  'either'(default), or 'midpoint'.}

  \item{feature.colnam}{column name from features.rd to be
  used for retrieving the nearest feature name. By default
  this is NULL assuming that features.rd has a column that
  includes the word 'name' somewhere in it.}

  \item{relativeTo}{calculate distance relative to query or
  subject. Default is 'subject'. See documentation of
  \code{\link{getNearestFeature}} for more information.}
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
\itemize{ \item When side='midpoint', the distance to
nearest feature is calculated by (start+stop)/2.  \item For
cases where a position is at the edge and there are no
feature up/down stream since it would fall off the
chromosome, the function simply returns 'NA'.  \item If
there are multiple locations where a query falls into, the
function arbitrarily chooses one to serve as the nearest
feature, then reports 2 upstream & downstream feature. That
may occasionally yield features which are the same upstream
and downstream, which is commonly encountered when studying
spliced genes or phenomena related to it.  \item Try not to
use this function for >50 spaces/seqnames/chromosomes
unless you have tons fo memory.  \item If strand
information doesn't exist, then everything is defaults to
'+' orientation (5' -> 3') \item If parallel=TRUE, then be
sure to have a parallel backend registered before running
the function. One can use any of the following libraries
compatible with \code{\link[foreach]{foreach}}: doMC,
doSMP, doSNOW, doMPI, doParallel. For example:
library(doMC); registerDoMC(2) }
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

