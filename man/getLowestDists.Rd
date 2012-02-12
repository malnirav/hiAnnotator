\name{getLowestDists}
\alias{getLowestDists}
\title{Get the lowest distance from the 5' or 3' boundaries of query and subject.}
\usage{
  getLowestDists(query = NULL, subject = NULL,
    subjectOrt = NULL, ok.chrs = NULL, res.nrst = NULL,
    side = "either", cores.use = 1)
}
\arguments{
  \item{query}{Compressed RangedList object to be used as
  the query.}

  \item{subject}{Compressed RangedList object to be used as
  the subject.}

  \item{subjectOrt}{CompressedCharacterList of
  orientation/strand of features in subject broken up by
  ok.chrs...generated using values(subject)[,"strand"]}

  \item{ok.chrs}{vector of spaces or chromosomes present in
  query/subject to do the calculations for.}

  \item{res.nrst}{nearest indices as returned by
  \code{\link{nearest}}.}

  \item{side}{boundary of subject/annotation to use to
  calculate the nearest distance. Options are '5p','3p', or
  the default 'either'.}

  \item{cores.use}{number of cores to utilize for the
  calculation. Default is 1. If >1, then a parallel backend
  must be registered to perform calculation with
  \code{\link{foreach}}.}
}
\value{
  a list of lowest distances named by ok.chrs
}
\description{
  Given a query and subject with indicies from
  \code{\link{nearest}}, calculate the shortest distance to
  either boundaries of the query and subject. This is a
  helper function utilized in
  \code{\link{getNearestFeature}},
  \code{\link{get2NearestFeature}}
}
\seealso{
  \code{\link{getNearestFeature}},
  \code{\link{get2NearestFeature}}.
}

