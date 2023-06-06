# rtmsSpectrum.R

#' Get an RTMS Spectrum
#'
#' Fetches an object of class `rtmsSpectrum` from the specified object.
#'
#' @param x The object from which the spectrum should be retrieved
#'
#' @param ... Other possible arguments to specify a particular spectrum to be
#' retrieved
#'
#' @returns An spectrum object of class `rtmsSpectrum`
#'
#' @export
getSpectrum <- function(x,...) {
	UseMethod("getSpectrum")
}

#' @export
getSpectrum.default <- function(x,...) {
	rtmsSpectrum(numeric(),numeric())
}

#' @param x The spectrum object
#' @param peaks A list of peak objects
#' @param ... Additional parameters
#'
#' @describeIn getSampleFromSpectrum The S3 method `getSample` for objects
#' of class `rtmsSpectrum`; calls `getSampleFromSpectrum`
#' @export
getSample.rtmsSpectrum <- function(x,peaks,...) {
	getSampleFromSpectrum(x,peaks,...)
}

#' Extract a sample from an RTMS spectrum object
#'
#' Extracts a sample object of class `rtmsSample` from a single spectrum object
#' of class `rtmsSample` usin g a list of peaks
#'
#' @param spectrum A full spectrum of class `rtmsSpectrum`
#'
#' @param peaks A list of peak objects of class `rtmsPeak`
#' @param freqSpacing If TRUE (the default), local maxima (estimated via
#' quadratic interpolation) are calculated in inverse m/z (or frequency) space,
#' as in FTMS spectra.  If FALSE, maxima are calculated directy in m/z space
#' @param threshold If NULL, all local maxima will be returned for each
#' subsample; if set to particular value, only those maxima above that threshold
#' will be returned.
#'
#' @returns An RTMS sample object of class `rtmsSample`
#'
#' @export
#' @examples
#' peaks <- rtmsPeakList(c(1516.83,1530.84),peakWidth=0.2,windowWidth = c(5,10))
#' names(peaks) <- c("Product","Substrate")
#' sample <- getSample(exampleSpectrum,peaks)
getSampleFromSpectrum <- function(spectrum,peaks,freqSpacing=TRUE,threshold=NULL) {
	subsamples <- lapply(peaks, function(p) getSubsampleFromSpectrum(spectrum,p,freqSpacing,threshold))
	new_rtmsSample(subsamples,peaks)
}

#' Create a new RTMS spectrum
#'
#' Generates an RTMS spectrum object (of class `rtmsSpectrum`) from a given
#' vector of m/z and intensity values.
#'
#' @param mz A numeric vector of m/z values
#'
#' @param intensity A numeric vector of intensity values
#'
#' @returns An object of class `rtmsSpectrum` with the given m/z and intensity
#' values
#'
#' @export
rtmsSpectrum <- function(mz,intensity) {
	mz <- as.double(mz)
	intensity <- as.double(intensity)
	intensity <- intensity[order(mz)]
	mz <- mz[order(mz)]
	validate_rtmsSpectrum(new_rtmsSpectrum(mz,intensity))
}

getSubsampleFromSpectrum <- function(spectrum,peak,freqSpacing=TRUE,threshold=NULL) {
	prel <- which(spectrum$mz >= peak$bounds[1] & spectrum$mz <= peak$bounds[2])
	peakPiece <- new_rtmsSpectrumPiece(spectrum$mz[prel],spectrum$intensity[prel])
	maxima <- rtms_pickQifftPeaks(peakPiece,freqSpacing,threshold)
	if (!is.null(peak$window)) {
		wrel <- which(spectrum$mz >= peak$window[1] & spectrum$mz <= peak$window[2])
		windowPiece <- new_rtmsSpectrumPiece(spectrum$mz[wrel],spectrum$intensity[wrel])
	} else { windowPiece <- NULL }
	new_rtmsSubsample(peakPiece,maxima,windowPiece)
}

rtms_pickQifftPeaks <- function(piece,freqSpacing=TRUE,threshold=NULL) {
	yvec <- piece$intensity
	if (length(yvec)<=2) { return(new_rtmsSpectrumPiece(numeric(),numeric())) }
	xvec <- piece$mz
	if (freqSpacing) { xvec <- 1/xvec }
	yvec <- yvec[order(xvec)]
	xvec <- xvec[order(xvec)]
	dxm1 <- c(0,diff(xvec))
	dxp1 <- c(diff(xvec),0)
	dym1 <- c(-1,diff(yvec))
	dyp1 <- c(-diff(yvec),-1)
	rel <- which(dym1>0 & dyp1>0)
	if (length(rel)==0) { return(new_rtmsSpectrumPiece(numeric(),numeric())) }

	um1 <- dxm1[rel]
	up1 <- dxp1[rel]
	rm1 <- dym1[rel]/um1
	rp1 <- dyp1[rel]/up1

	a <- rm1+rp1
	b <- rm1*up1-rp1*um1
	c <- um1+up1

	mz <- xvec[rel]+b/(2*a)
	intensity <- yvec[rel]+(b^2)/(4*a*c)
	if (freqSpacing) { mz <- 1/mz }
	intensity <- intensity[order(mz)]
	mz <- mz[order(mz)]
	if (!is.null(threshold)) {
		mz <- mz[intensity>threshold]
		intensity <- intensity[intensity>threshold]
	}
	return(new_rtmsSpectrumPiece(mz,intensity))
}

validate_rtmsSpectrum <- function(piece) {
	if (length(piece$mz)!=length(piece$intensity)) {
		stop("'mz' and 'intensity' must be numeric vectors of the same length",call. = FALSE)
	}
	piece
}

new_rtmsSpectrum <- function(mz,intensity) {
	stopifnot(is.double(mz))
	stopifnot(is.double(intensity))

	structure(list(mz=mz,intensity=intensity),
			  class="rtmsSpectrum")
}
