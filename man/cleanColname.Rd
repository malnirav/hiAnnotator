\name{cleanColname}
\alias{cleanColname}
\title{Clean the supplied string from punctuations and spaces.}
\usage{
  cleanColname(x, description = "colnam")
}
\arguments{
  \item{x}{string or a vector to be cleaned.}

  \item{description}{OPTIONAL string identifying the
  purpose supplied string in x to be displayed in the
  cleaning message.}
}
\value{
  cleaned string or a vector.
}
\description{
  Function to clean the supplied string from punctuations
  and spaces so it can be used as column headings.
}
\examples{
cleanColname("HIV-test")
cleanColname("HIV*test")
cleanColname("HIV-test","myAlias")
}
\seealso{
  \code{\link{getFeatureCounts}},
  \code{\link{makeRangedData}},
  \code{\link{getNearestFeature}},
  \code{\link{getSitesInFeature}}.
}

