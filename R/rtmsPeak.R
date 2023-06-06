# rtmsPeak.R

#' Create an RTMS m/z Peak Object
#'
#' Generates an object of class `rtmsPeak` which contains the m/z values
#' bounding a spectrometric peak to be measured.  A peak object specifies not
#' only the m/z value at the cetner of the peak, but the upper and lower bounds
#' within which the peak is to be quantified; it also may optionally include
#' wider upper and lower bounds used to plot the peak in a wider context of the
#' spectrum.
#'
#' @param value The m/z value that the peak is intended to measure
#'
#' @param peakWidth The width of the peak centered on `value`. If a single
#' numeric value, the lower bound of the peak will lie `peakWidth/2` below
#' `value`,and the upper bound will lie `peakWidth/2` above `value`.  If a
#' vector of two numeric values, the first value specifies how far below `value`
#' the lower bound lies, and the second value specifies how far above `value`
#' the upper bound lies.  If parameter `bounds` is not null, this parameter will
#' be ignored.
#' @param windowWidth The width of the optional wider window around `value`
#' used to show the peak in context.  Operates by the same principles as
#' `peakWidth` with a single value split evenly between lower and upper bounds,
#' and two values specifying how far below and above `value` each bound lies. If
#' parameter `window` is not null, this parameter will be ignored.
#' @param bounds If not null, a two-value numeric vector specifying the lower
#' and upper m/z bounds of the measured peak. One of `bounds` or `peakWidth`
#' must be not null, and if `bounds` is not null, then `peakWidth` will be
#' ignored.
#' @param window If not null, a two-value numeric vector specifying the lower
#' and upper m/z bounds of the wider context window of the peak. If `window` is
#' not null, then `windowWidth` will be ignored.
#'
#' @returns An object of class `rtmsPeak`, used by RTMS functions to extract
#' and measure peaks from mass spectra.
#'
#' @export
#' @examples
#' peaks <- rtmsPeak(1516.83,peakWidth=0.2,windowWidth = c(5,10))
rtmsPeak <- function(value,peakWidth=0.1,windowWidth=NULL,bounds=NULL,window=NULL) {
	value <- as.double(value)
	if (length(value)!=1) { stop("Argument 'value' must be a single numeric value.") }
	if (!is.null(bounds)) {
		bounds <- as.double(bounds)
		if (length(bounds)!=2 || bounds[1]>value || bounds[2]<value) {
			stop("Argument 'bounds', if specified, must be a pair of numeric values lying on either side of 'value'.")
		}
	} else if (!is.null(peakWidth)) {
		peakWidth <- as.double(peakWidth)
		if (length(peakWidth)==1) { peakWidth <- rep(peakWidth,2)/2 }
		if (length(peakWidth)!=2 || any(peakWidth<0)) {
			stop("Argument 'peakWidth' must be a single non-negative numeric value, or a pair of non-negative numeric values.")
		}
		bounds <- value+(c(-1,1)*peakWidth)
	} else {
		stop("Function call must specifiy either 'peakWidth' or 'bounds'.")
	}
	if (!is.null(window)) {
		window <- as.double(window)
		if (length(window)!=2 || window[1]>value || window[2]<value) {
			stop("Argument 'window', if specified, must be a pair of numeric values lying on either side of 'value'.")
		}
	} else if (!is.null(windowWidth)) {
		windowWidth <- as.double(windowWidth)
		if (length(windowWidth)==1) { windowWidth <- rep(windowWidth,2)/2 }
		if (length(windowWidth)!=2 || any(windowWidth<0)) {
			stop("Argument 'windowWidth' must be a single non-negative numeric value, or a pair of non-negative numeric values.")
		}
		window <- value+(c(-1,1)*windowWidth)
	}

	return(new_rtmsPeak(value,bounds,window))
}

#' Create a list of RTMS m/z peak objects
#'
#' Generates a list of objects of class `rtmsPeak` which can be used to extract
#' a sample or sample set from other RTMS objects.
#'
#' @param values The m/z values that the peaks is intended to measure
#'
#' @param peakWidth The width of each peak centered on `values`. If a single
#' numeric value, the lower bound of each peak will lie `peakWidth/2` below
#' the given m/z value ,and the upper bound will lie `peakWidth/2` above it.
#' If a vectorof two numeric values, the first value specifies how far below
#' each given m/z value the lower bounds lie, and the second value specifies how
#' far above each value the upper bounds lie.
#' @param windowWidth The width of each optional wider window around the m/z
#' values used to show the peaks in context.  Operates by the same principles as
#' `peakWidth` with a single value split evenly between lower and upper bounds,
#' and two values specifying how far below and above the m/z values each bound
#' lies.
#'
#' @returns A list of objects of class `rtmsPeak`
#'
#' @export
#' @examples
#' peaks <- rtmsPeakList(c(1516.83,1530.84),peakWidth=0.2,windowWidth = c(5,10))
#' names(peaks) <- c("Product","Substrate")
rtmsPeakList <- function(values,peakWidth=0.1,windowWidth=NULL) {
	lapply(values,function(p) rtmsPeak(p, peakWidth, windowWidth, bounds=NULL, window=NULL))
}

new_rtmsPeak <- function(value,bounds,window=NULL) {
	stopifnot(is.double(value))
	stopifnot(is.double(bounds))
	if (!is.null(window)) { stopifnot(is.double(window)) }

	structure(list(value=value,bounds=bounds,window=window),
			  class="rtmsPeak")
}

checkPeaksEqual <- function(peak1,peak2) {
	if ((peak1$value != peak2$value) || any(peak1$bounds!=peak2$bounds)) {
		return(FALSE)
	}
	if ((is.null(peak1$window)!=is.null(peak2$window)) ||
		(!is.null(peak1$window) && any(peak1$window!=peak2$window))) {
		return(FALSE)
	}
	return(TRUE)
}

checkPeakSetsEqual <- function(peakset1,peakset2) {
	if (length(peakset1)!=length(peakset2)) { return(FALSE) }
	for (index in seq_along(peakset1)) {
		if (!checkPeaksEqual(peakset1[[index]],peakset2[[index]])) {
			return(FALSE)
		}
	}
	return(TRUE)
}
