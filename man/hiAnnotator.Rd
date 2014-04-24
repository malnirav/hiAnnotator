\docType{package}
\name{hiAnnotator}
\alias{hiAnnotator}
\alias{hiAnnotator-package}
\title{Annotating GRanges objects with hiAnnotator.}
\description{
hiAnnotator contains set of functions which allow users to
annotate a GRanges object with custom set of annotations.
The basic philosophy of this package is to take two GRanges
objects (query & subject) with common set of seqnames (i.e.
chromosomes) and return associated annotation per seqnames
and rows from the query matching seqnames and rows from the
subject (i.e. genes or cpg islands). The package comes with
three types of annotation functions which calculates if a
position from query is: within a feature, near a feature,
or count features in defined window sizes. Moreover, one
can utilize parallel backend for each annotation function
to utilize the foreach package. In addition, the package is
equipped with wrapper functions, which finds appropriate
columns needed to make a GRanges object from a common
dataframe.
}
\author{
Nirav V Malani
}

