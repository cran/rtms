# rtmsMeasurement.R

#' Measure peaks and samples in an RTMS sample set
#'
#' `measureSampleSet()` extracts one or more measurements for every peak in
#' every sample in an RTMS sample set object (of class `rtmsSampleSet`).
#'
#' @param sampleset An object of class `rtmsSampleSet`
#'
#' @param measure A character vector of named measurements, or a list of
#' custom measurement functions.  Supported measurement names are
#' "PeakIntensity", which takes the total of any local maxima within the peak
#' width, "PeakArea", which takes the area under the intensity curve within the
#' peak width, and "NumPeaks", which counts the local maxima in the peak window.
#' If `measure` is a list of functions, each function must take an object of
#' class `rtmsSubsample`, and return a single numeric value.  If the functions
#' are named, those names will be returned in the "measure" column of the
#' resulting data frame; otherwise they will be identified as "Measure1",
#' "Measure2", etc.
#'
#' @returns A data frame with one row for each sample, peak, and measurement.
#' The data.frame will have a character column named "sample", containing either
#' the name of the sample (if the samples in `sampleset` are named) or the
#' index of the sample if they are not (but it will always be a character
#' column); a column named "peakName" with the name of the relevant peak  (if
#' the "peaks" attribute of `sampleset` is a named list); a column named
#' "peakValue" containing the m/z value at the center of the relevant peak; a
#' column named "measure" containing the name of the relevant measure; and a
#' column named "value" containing the numeric value of the particular measure
#' for that sample and peak.
#'
#' @export
#' @examples
#' peaks <- rtmsPeakList(c(1516.83,1530.84),peakWidth=0.2,windowWidth = c(5,10))
#' names(peaks) <- c("Product","Substrate")
#' sample <- getSample(exampleSpectrum,peaks)
#' sampleSet <- rep(sample,3)
#' names(sampleSet) <- c("A","B","C")
#'
#' measures <- measureSampleSet(sampleSet)
measureSampleSet <- function(sampleset,measure="PeakIntensity") {
	stopifnot(inherits(sampleset,"rtmsSampleSet"))

	if (is.null(names(sampleset))) { sampleNames <- as.character(seq_along(sampleset)) }
	else { sampleNames <- names(sampleset) }
	measuredf <- data.frame()
	for (i in seq_along(sampleset)) {
		curdf <- measureSample(sampleset[[i]],measure)
		curdf$sample <- sampleNames[i]
		measuredf <- rbind(measuredf,curdf)
	}
	measuredf
}

#' Measure peaks in an RTMS sample
#'
#' `measureSample()` extracts one or more measurements for every peak in
#' an RTMS sample object (of class `rtmsSample`).
#'
#' @param sample An object of class `rtmsSample`
#'
#' @param measure A character vector of named measurements, or a list of
#' custom measurement functions.  Supported measurement names are
#' "PeakIntensity", which takes the total of any local maxima within the peak
#' width, "PeakArea", which takes the area under the intensity curve within the
#' peak width, and "NumPeaks", which counts the local maxima in the peak window.
#' If `measure` is a list of functions, each function must take an object of
#' class `rtmsSubsample`, and return a single numeric value.  If the functions
#' are named, those names will be returned in the "measure" column of the
#' resulting data frame; otherwise they will be identified as "Measure1",
#' "Measure2", etc.
#'
#' @returns A data frame with one row for each peak and measurement in the
#' sample. The data.frame will have a column named "peakName" with the name of
#' the relevant peak  (if  the "peaks" attribute of `sample` is a named list);
#' a column named "peakValue" containing the m/z value at the center of the
#' relevant peak; a column named "measure" containing the name of the relevant
#' measure; and a column named "value" containing the numeric value of the
#' particular measure for that peak.
#'
#' @export
#' @examples
#' peaks <- rtmsPeakList(c(1516.83,1530.84),peakWidth=0.2,windowWidth = c(5,10))
#' names(peaks) <- c("Product","Substrate")
#' sample <- getSample(exampleSpectrum,peaks)
#'
#' measure <- measureSample(sample,c("PeakArea","PeakIntensity"))
#'
#' myFunctions <- list(PeakRawIntensity = function(s) max(s$peakPiece$intensity))
#' myMeasures <- measureSample(sample,myFunctions)
measureSample <- function(sample,measure="PeakIntensity") {
	measure <- getRtmsMeasures(measure)

	peaks <- attr(sample,"peaks")
	peakValues <- vapply(peaks,function(p) p$value,0.0)
	measuredf <- data.frame()
	basedf <- data.frame(peakValue=peakValues,peakIndex=seq_along(peakValues))
	if (!is.null(names(peaks))) {
		basedf$peakName <- names(peaks)
	}
	for (i in seq_along(measure)) {
		curMeasure <- measure[[i]]
		measureValue <- vapply(sample,curMeasure,0)
		currentdf <- basedf
		currentdf$measureIndex <- i
		currentdf$measure <- names(measure)[[i]]
		currentdf$value <- measureValue
		measuredf <- rbind(measuredf,currentdf)
	}
	measuredf <- measuredf[order(measuredf$peakIndex,measuredf$measureIndex),]
	measuredf$peakIndex <- NULL
	measuredf$measureIndex <- NULL
	measuredf

}

getRtmsMeasures <- function(measure) {
	if (is.character(measure)) {
		measure <- unique(measure)
		newMeasure <- list()
		validMeasure <- c()
		for (i in seq_along(measure)) {
			validMeasure[i] <- TRUE
			if (measure[[i]]=="PeakIntensity") {
				newMeasure[[i]] <- rtmsMeasure_PeakIntensity
			} else if (measure[[i]]=="PeakArea") {
				newMeasure[[i]] <- rtmsMeasure_PeakArea
			} else if (measure[[i]]=="NumPeaks") {
				newMeasure[[i]] <- rtmsMeasure_NumPeaks
			} else {
				validMeasure[[i]] <- FALSE
			}
		}
		names(newMeasure) <- measure
		newMeasure <- newMeasure[validMeasure]
		if (length(newMeasure)==0) {
			stop("No valid measure names were included, Supported named",
				 " measures are 'PeakIntensity', 'PeakArea', and 'NumPeaks'.",
				 " Any other measures must be passed as functions.")
		} else if (any(!validMeasure)) {
			warning("The following measure names are not supported and will be ",
					sprintf("ignored: %s",paste(measure[!validMeasure],collapse=", ")))
		}
	} else if (is.function(measure)) {
		newMeasure <- list(Measure1=measure)
	} else if (is.list(measure) && length(measure)>0 && all(vapply(measure,is.function,FALSE))) {
		newMeasure <- measure
		if (is.null(names(measure))) {
			names(newMeasure) <- paste0("Measure",seq_along(measure))
		}
	} else { stop("Invalid measure parameter.") }
	return(newMeasure)
}
rtmsMeasure_PeakIntensity <- function(subsample) {
	stopifnot(inherits(subsample,"rtmsSubsample"))
	return(sum(subsample$maxima$intensity))
}
rtmsMeasure_PeakArea <- function(subsample) {
	stopifnot(inherits(subsample,"rtmsSubsample"))
	return(diff(range(subsample$peakPiece$mz))*mean(subsample$peakPiece$intensity))
}
rtmsMeasure_NumPeaks <- function(subsample) {
	stopifnot(inherits(subsample,"rtmsSubsample"))
	return(length(subsample$maxima$mz))
}
