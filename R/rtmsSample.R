# rtmsSample.R

#' Get an RTMS Sample
#'
#' Fetches an object of class `rtmsSample` from the specified object.
#'
#' @param x The object from which the sample should be retrieved
#'
#' @param peaks A list of objects of class `rtmsPeak`
#'
#' @param ... Other possible arguments to specify a particular sample to be
#' retrieved
#'
#' @returns A sample object of class `rtmsSample`
#'
#' @export
getSample <- function(x,peaks,...) {
	UseMethod("getSample")
}
#' @export
getSample.default <- function(x,peaks,...) {
	sample <- lapply(peaks,emptySubsample)
	new_rtmsSample(sample,peaks)
}

#' Get an RTMS Sample Set
#'
#' Fetches an object of class `rtmsSampleSet` from the specified object.
#'
#' @param x The object from which the sample set should be retrieved
#'
#' @param peaks A list of objects of class `rtmsPeak`
#'
#' @param ... Other possible arguments to specify a particular sample set to be
#' retrieved
#'
#' @returns A sample set object of class `rtmsSampleSet`
#'
#' @export
getSampleSet <- function(x,peaks,...) {
	UseMethod("getSampleSet")
}
#' @export
getSampleSet.default <- function(x,peaks,...) {
	new_rtmsSampleSet(list(),peaks)
}

#' Generate an empty RTMS sample set
#'
#' @description
#' Produces a sample set (of class `rtmsSampleSet`) with no samples but a
#' specific peaks attribute.  Useful for building a sample set one sample at
#' a time.
#'
#' @param peaks A list of objects of class `rtmsPeak`
#'
#' @returns An empty object of class `rtmsSampleSet`
#'
#' @export
emptySampleSet <- function(peaks) {
	new_rtmsSampleSet(list(),peaks)
}

#' @name sampleAndSampleSet
#' @title Functions for creating and manipulating samples and sample sets
#'
#' @description
#' Select a subset of a sample set (returns an `rtmsSampleSet`)
#'
#' @details
#' The sample (class `rtmsSample`) and sample set (class `rtmsSampleSet`)
#' objects are the core structures used to extract meaningful data from mass
#' spectographic data.  In general, samples and sample sets will be created
#' automatically from other RTMS objects (such as readers or spectra) but in
#' the event that one wishes to manipulate them directly, it is important to
#' understand several details about how they work.
#'
#' In terms of the data it contains, an object of class `rtmsSample` is just a
#' list of smaller objects (of class `rtmsSubsample`); however, each of these
#' subsamples corresponds to an `rtmsPeak` object that was used to extract it;
#' the `rtmsSample` object therefor has a "peaks" attribute, which is a list of
#' objects of class `rtmsPeak` corresponding to the subsamples in the
#' `rtmsSample` object.  This attribute is used to determine how measurements
#' of the sample are reported and how the sample is plotted.
#'
#' Similarly, the data contained in an object of class `rtmsSampleSet` is just
#' a list of `rtmsSample` objects but with an important difference.  If many
#' `rtmsSample` objects were arranged into a list, there would be no guarantee
#' that they contain measurements of the same peaks; such guarantees are
#' essential for plotting sample sets together or constructing extracted ion
#' chromatograms.  The `rtmsSampleSet` therefore strips the "peaks" attribute
#' from its individual members, and applies a single shared "peaks" attribute
#' to the entire sample set.  Further samples can only be added to the sample
#' set if their peaks attributes are deemed compatible.
#'
#' @param x An object of class `rtmsSampleSet`
#'
#' @param i An index or set of indices
#' @param ... Included for S3 compatibility
#'
#' @returns An object of class `rtmsSampleSet`
#'
#' @export
`[.rtmsSampleSet` <- function(x, i, ...) {
	new_rtmsSampleSet(NextMethod(),peaks=attr(x,"peaks"))
}
#' @param i A single numeric index of the sample set
#' @describeIn sampleAndSampleSet Select a single element of a sample set
#' (returns an `rtmsSample`)
#'
#' @returns An object of class `rtmsSample`
#' @export
`[[.rtmsSampleSet` <- function(x, i, ...) {
	new_rtmsSample(NextMethod(),peaks=attr(x,"peaks"))
}
#' @param value An object of class `rtmsSample`
#' @describeIn sampleAndSampleSet Insert a sample into a sample set
#'
#' @returns An object of class `rtmsSampleSet`
#' @export
`[[<-.rtmsSampleSet` <- function(x, i, value) {
	stopifnot(checkPeakSetsEqual(attr(x,"peaks"),attr(value,"peaks")))
	new_rtmsSampleSet(NextMethod(),peaks=attr(x,"peaks"))
}
#' @describeIn sampleAndSampleSet Repeat a sample set multiple times
#'
#' @returns An object of class `rtmsSampleSet`
#' @export
rep.rtmsSampleSet <- function(x, ...) {
	new_rtmsSampleSet(NextMethod(),peaks=attr(x,"peaks"))
}
#' @param x An object of class `rtmsSample`
#' @describeIn sampleAndSampleSet Create a sample set by repeating a single
#' sample (returns an `rtmsSampleSet`)
#'
#' @returns An object of class `rtmsSampleSet`
#' @export
rep.rtmsSample <- function(x,...) {
	new_rtmsSampleSet(rep(list(unclass(x)),...),peaks=attr(x,"peaks"))
}

rtmsSpectrumPiece <- function(mz,intensity) {
	mz <- as.double(mz)
	intensity <- as.double(intensity)
	intensity <- intensity[order(mz)]
	mz <- mz[order(mz)]
	validate_rtmsSpectrumPiece(new_rtmsSpectrumPiece(mz,intensity))
}

emptySubsample <- function(x) {
	peakPiece <- new_rtmsSpectrumPiece(numeric(),numeric())
	maxima <- new_rtmsSpectrumPiece(numeric(),numeric())
	new_rtmsSubsample(peakPiece,maxima)
}

new_rtmsSampleSet <- function(subsampleset,peaks) {
	stopifnot(all(vapply(peaks,function(p) inherits(p,"rtmsPeak"),FALSE)))
	stopifnot(all(vapply(subsampleset,function(ss)
		all(vapply(ss, function(s) inherits(s,"rtmsSubsample"),FALSE)),
		FALSE)))

	structure(subsampleset,
			  peaks=peaks,
			  class="rtmsSampleSet")
}

new_rtmsSample <- function(subsamples,peaks) {
	stopifnot(all(vapply(peaks,function(p) inherits(p,"rtmsPeak"),FALSE)))
	stopifnot(all(vapply(subsamples,function(s) inherits(s,"rtmsSubsample"),FALSE)))

	structure(subsamples,
			  peaks=peaks,
			  class="rtmsSample")
}
validate_rtmsSample <- function(sample) {
	if (length(attr(sample,"peaks"))!=length(sample)) {
		stop("'peaks' and 'samples' must be lists of the same length",call. = FALSE)
	}
}

new_rtmsSubsample <- function(peakPiece,maxima,windowPiece=NULL) {
	stopifnot(inherits(peakPiece,"rtmsSpectrumPiece"))
	stopifnot(inherits(maxima,"rtmsSpectrumPiece"))
	if (!is.null(windowPiece)) { stopifnot(inherits(peakPiece,"rtmsSpectrumPiece")) }

	structure(list(peakPiece=peakPiece,maxima=maxima,windowPiece=windowPiece),
			  class="rtmsSubsample")
}
validate_rtmsSpectrumPiece <- function(piece) {
	if (length(piece$mz)!=length(piece$intensity)) {
		stop("'mz' and 'intensity' must be numeric vectors of the same length",call. = FALSE)
	}
	piece
}

new_rtmsSpectrumPiece <- function(mz,intensity) {
	stopifnot(is.double(mz))
	stopifnot(is.double(intensity))

	structure(list(mz=mz,intensity=intensity),
			  class="rtmsSpectrumPiece")
}
