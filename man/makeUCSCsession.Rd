\name{makeUCSCsession}
\alias{makeUCSCsession}
\title{Initiate UCSC genome browser session given the freeze argument.}
\usage{
  makeUCSCsession(freeze = "hg18")
}
\arguments{
  \item{freeze}{one of following: hg18, mm8, rheM, etc.
  Default is hg18.}
}
\value{
  browser session object compatible with rtracklayer
  functions.
}
\description{
  Initiate UCSC genome browser session given the freeze
  argument.
}
\examples{
session<-makeUCSCsession()
genome(session)
session<-makeUCSCsession("mm8")
genome(session)
}
\seealso{
  \code{\link{getUCSCtable}}, \code{\link{makeRangedData}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}}.
}

