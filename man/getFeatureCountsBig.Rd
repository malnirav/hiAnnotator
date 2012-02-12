\name{getFeatureCountsBig}
\alias{getFeatureCountsBig}
\title{Get counts of annotation within a defined window around each query range/position for large annotation objects spanning greater than 100 million rows. This is still in beta phase.}
\usage{
  getFeatureCountsBig(sites.rd, features.rd, colnam = NULL,
    widths = c(1000, 10000, 1e+06))
}
\arguments{
  \item{sites.rd}{RangedData object to be used as the
  query.}

  \item{features.rd}{RangedData obj to be used as the
  subject or the annotation table.}

  \item{colnam}{column name to be added to sites.rd for the
  newly calculated annotation...serves as a prefix to
  windows sizes!}

  \item{widths}{named vector of window sizes to be used for
  casting a net around each position. Default:
  \code{c(1000,10000,1000000)}}
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
  8,17. This new range is then compared against the
  midpoints of subject to find any overlapping ranges and
  then tallied up. Note that here counting is doing using
  midpoint of the ranges in subject instead of start-stop
  boundaries.
}
\seealso{
  \code{\link{makeRangedData}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}},
  \code{\link{getFeatureCounts}}.
}

