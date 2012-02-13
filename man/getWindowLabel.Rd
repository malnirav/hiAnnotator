\name{getWindowLabel}
\alias{getWindowLabel}
\title{Generate a window size label.}
\usage{
  getWindowLabel(x)
}
\arguments{
  \item{x}{vector of integers to generate the labels for.}
}
\value{
  a character vector of length(x) which has x normalized
  and suffixed by bp, Kb, Mb, or Gb depending on respective
  interval sizes.
}
\description{
  Function to generate asthetically pleasing window size
  label given an integer. This is one of the helper
  function used in \code{\link{getFeatureCounts}} &
  \code{\link{getFeatureCountsBig}}.
}
\examples{
# Convert a data frame to RangedData object
getWindowLabel(c(0,1e7,1e3,1e6,2e9))
}
\seealso{
  \code{\link{getFeatureCounts}},
  \code{\link{makeRangedData}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}}.
}

