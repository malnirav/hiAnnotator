\name{getLowestDists}
\alias{getLowestDists}
\title{Get the lowest biological distance from the 5' or 3' boundaries of query and subject.}
\usage{
  getLowestDists(query = NULL, subject = NULL,
    res.nrst = NULL, side = "either")
}
\arguments{
  \item{query}{GRanges object to be used as the query which
  holds data for 'queryHits' attribute of res.nrst.}

  \item{subject}{GRanges object to be used as the subject
  which holds data for 'subjectHits' attribute of
  res.nrst.}

  \item{res.nrst}{a dataframe of nearest indices as
  returned by \code{\link{nearest}}.}

  \item{side}{boundary of subject/annotation to use to
  calculate the nearest distance. Options are '5p','3p', or
  the default 'either'.}
}
\value{
  res.nrst with lowest distances appended at the end.
}
\description{
  Given a query and subject with indicies from
  \code{\link{nearest}}, calculate the shortest biological
  distance to either boundaries of the query and subject.
  This is a helper function utilized in
  \code{\link{getNearestFeature}},
  \code{\link{get2NearestFeature}}
}
\note{
  for cases where a query has multiple nearest neighbors or
  overlaps with >1 subjects, the function will choose the
  subject with the lowest absolute distance.
}
\seealso{
  \code{\link{getNearestFeature}},
  \code{\link{get2NearestFeature}}.
}

