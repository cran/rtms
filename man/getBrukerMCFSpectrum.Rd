% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtmsBrukerMCFReader.R
\name{getSpectrum.rtmsBrukerMCFReader}
\alias{getSpectrum.rtmsBrukerMCFReader}
\alias{getBrukerMCFSpectrum}
\title{Extract a spectrum from a Bruker MCF directory}
\usage{
\method{getSpectrum}{rtmsBrukerMCFReader}(x, ...)

getBrukerMCFSpectrum(reader, index)
}
\arguments{
\item{x}{The MCF reader object}

\item{...}{Additional parameters}

\item{reader}{An RTMS reader object of class \code{rtmsBrukerMCFReader}}

\item{index}{A single numeric index specifying which acquisition should be
extracted}
}
\value{
An RTMS spectrum object of class \code{rtmsSpectrum}
}
\description{
Extracts an RTMS spectrum object (of class \code{rtmsSpectrum}) from a multi-
acquisition Bruker MCF directory opened using an RTMS reader object (of
class \code{rtmsBrukerMCFReader}).  A numeric index is used to identify which
spectrum should be extracted.
}
\section{Functions}{
\itemize{
\item \code{getSpectrum(rtmsBrukerMCFReader)}: The S3 method \code{getSpectrum} for objects
of class \code{rtmsBrukerMCFReader}; calls \code{getBrukerMCFSpectrum}

}}
