% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtmsSample.R
\name{getSample}
\alias{getSample}
\title{Get an RTMS Sample}
\usage{
getSample(x, peaks, ...)
}
\arguments{
\item{x}{The object from which the sample should be retrieved}

\item{peaks}{A list of objects of class \code{rtmsPeak}}

\item{...}{Other possible arguments to specify a particular sample to be
retrieved}
}
\value{
A sample object of class \code{rtmsSample}
}
\description{
Fetches an object of class \code{rtmsSample} from the specified object.
}
