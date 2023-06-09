% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtmsBrukerBAFReader.R
\name{newBrukerBAFReader}
\alias{newBrukerBAFReader}
\title{Open a Bruker single-acquisition BAF directory}
\usage{
newBrukerBAFReader(bafdir)
}
\arguments{
\item{bafdir}{A directory (usually with the extension ".d") containing data
from a single Bruker acquisition. This directory will contain a file with
extension ".baf" that holds the primary raw data, as well as an index and
calibration file (see Details).}
}
\value{
An object of class \code{rtmsBrukerBAFReader} which can extract raw data
from the specified directory
}
\description{
Creates an RTMS reader object (of class \code{rtmsBrukerBAFReader}) which can
extract data from a Bruker single acquisition directory (extension ".d")
}
\details{
Currently, RTMS can create reader objects for two binary Bruker data formats,
BAF (presumably standing for "Bruker acquisition format") holding data from
a single spectrum acquisition, and MCF (probably "multiacquisition container
format") containing data from multiple spectra acquired in a single run.
Both formats hold data in a directory marked with the extension ".d". The
single acquisition BAF format directory contains three essential data files:
the main raw data file with extension ".baf", an index file with extension
".baf_idx", and a calibration data file with extension ".baf_xtr". This
function processes the index and calibration files so that raw data can be
extracted quickly on demand from the ".baf" file.

An important note: when a MCF multi-acquistion reader is created, it creates
an open connection to the raw data file, which allows for quicker processing
of many spectra in a single file.  However, because a BAF file contains only
a single spectrum, there is little advantage to maintaining an open
connection, so the connection is re-opened every time data is read.  Thus,
while it is important to close an MCF reader object when all data is
extracted, it is not necessary to close an object of class
\code{rtmsBrukerBAFReader}.
}
