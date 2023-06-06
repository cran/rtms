# rtmsPlotting.R

#' Plot an RTMS sample object
#'
#' `plotRtmsSample()` takes an RTMS sample object and produces a `ggplot` object
#' depicting all extracted peaks, and their context windows if included.
#'
#' @param sample An object of class `rtmsSample`.
#'
#' @param usePeakNames If the list of peaks used to create the sample was a
#' named list, then setting this to TRUE (the default) will use those names to
#' label the facets of the plotted sample. If set to FALSE, the facets will be
#' labelled with the m/z values of each peak.  This parameter will be ignored if
#' the peaks are unnamed.
#' @param freey If TRUE (the default) the y-axes of each peak's facet will be
#' allowed to vary freely, so different peaks will be plotted on different
#' scales.  Setting this to FALSE will fix all peaks with in a sample on the
#' same y-axis scale.
#' @returns A `ggplot` object depicting the RTMS sample.
#'
#' @export
#' @examples
#' peaks <- rtmsPeakList(c(1516.83,1530.84),peakWidth=0.2,windowWidth = c(5,10))
#' names(peaks) <- c("Product","Substrate")
#' sample <- getSample(exampleSpectrum,peaks)
#'
#' plot1 <- plotRtmsSample(sample)
#' plot2 <- plotRtmsSample(sample,freey=FALSE)
plotRtmsSample <- function(sample,usePeakNames=TRUE,freey=TRUE) {
	stopifnot(inherits(sample,"rtmsSample"))

	mz <- NULL
	intensity <- NULL
	facet <- NULL

	peaks <- attr(sample,"peaks")

	spectrum_df <- data.frame()
	peak_df <- data.frame()
	maxima_df <- data.frame(peak=c(),mz=c(),intensity=c(),type=c())
	for (i in seq_along(peaks)) {
		if (usePeakNames && !is.null(names(peaks)) &&  names(peaks)[[i]]!="") {
			peakName <- names(peaks)[[i]]
		} else {
			peakName <- paste0(peaks[[i]]$value," m/z")
		}
		curpiece_df <- data.frame(peak=peakName,
								  mz=sample[[i]]$peakPiece$mz,
								  intensity=sample[[i]]$peakPiece$intensity,
								  type="Peak")
		peak_df <- rbind(peak_df,data.frame(peak=peakName,
											mz=peaks[[i]]$value,
											intensity=0))
		if (!is.null(sample[[i]]$windowPiece)) {
			curwindow_df <- data.frame(peak=peakName,
									   mz=sample[[i]]$windowPiece$mz,
									   intensity=sample[[i]]$windowPiece$intensity,
									   type="Context")
			curpiece_df <- rbind(curpiece_df,curwindow_df)
		}
		spectrum_df <- rbind(spectrum_df,curpiece_df)

		if (length(sample[[i]]$maxima$mz)>0) {
			maxima_df <- rbind(maxima_df,data.frame(peak=peakName,
													mz=sample[[i]]$maxima$mz,
													intensity=sample[[i]]$maxima$intensity,
													type="Peak"))
		}
	}

	facet_df <- data.frame(peak=rep(peak_df$peak,each=length(unique(spectrum_df$type))),
						   type=rep(unique(spectrum_df$type,times=length(peaks))))
	facet_df$facet <- paste0(facet_df$peak," (",facet_df$type,")")

	spectrum_df <- merge(spectrum_df,facet_df,by=c("peak","type"))
	maxima_df <- merge(maxima_df,facet_df,by=c("peak","type"))
	peak_df <- merge(peak_df,facet_df,by="peak")
	spectrum_df$facet <- factor(spectrum_df$facet,facet_df$facet)
	maxima_df$facet <- factor(maxima_df$facet,facet_df$facet)
	peak_df$facet <- factor(peak_df$facet,facet_df$facet)

	# Generate facetted ggplot
	if (freey) { fscale <- "free" }
	else { fscale <- "free_x" }
	gplot <- ggplot2::ggplot(spectrum_df,ggplot2::aes(x=mz,y=intensity))+
		ggplot2::geom_hline(yintercept=0,colour="black")+
		ggplot2::geom_vline(ggplot2::aes(xintercept=mz),colour="red",data=peak_df)
	if (nrow(maxima_df)>0) {
		gplot <- gplot+
			ggplot2::geom_vline(ggplot2::aes(xintercept=mz),colour="blue",linetype=2,data=maxima_df)
	}
	gplot <- gplot+
		ggplot2::geom_line()+
		ggplot2::scale_y_continuous(limits=c(0,NA),labels=function(v) sprintf("%sM",v/1e6))+
		ggplot2::labs(x="m/z",y="Intensity")+
		ggplot2::facet_wrap(ggplot2::vars(facet),ncol=length(unique(facet_df$type)),scale=fscale,drop=FALSE)
	gplot
}

#' Plot an RTMS sample set object
#'
#' `plotRtmsSampleSet()` takes an RTMS sample set object and produces a `ggplot`
#' object depicting all extracted peaks, and their context windows if included.
#'
#' @param sampleset An object of class `rtmsSampleSet`.
#'
#' @param usePeakNames If the list of peaks used to create the sample set was a
#' named list, then setting this to TRUE (the default) will use those names to
#' label the facets of the plotted sample set. If set to FALSE, the facets will
#' be labelled with the m/z values of each peak.  This parameter will be ignored
#' if the peaks are unnamed.
#' @param freey If TRUE (the default) the y-axes of each sample and peak's facet
#' will beallowed to vary freely, so different facets will be plotted on
#' different scales.  Setting this to FALSE will fix all peaks and samples
#' on the same y-axis scale.
#' @returns A `ggplot` object depicting the RTMS sample set.
#'
#' @export
#' @examples
#' peaks <- rtmsPeakList(c(1516.83,1530.84),peakWidth=0.2,windowWidth = c(5,10))
#' names(peaks) <- c("Product","Substrate")
#' sample <- getSample(exampleSpectrum,peaks)
#' sampleSet <- rep(sample,3)
#' names(sampleSet) <- c("A","B","C")
#'
#' plot1 <- plotRtmsSampleSet(sampleSet)
#' plot2 <- plotRtmsSampleSet(sampleSet,freey=FALSE) + ggplot2::theme_bw()
plotRtmsSampleSet <- function(sampleset,usePeakNames=TRUE,freey=TRUE) {
	stopifnot(inherits(sampleset,"rtmsSampleSet"))

	mz <- NULL
	intensity <- NULL
	facet <- NULL

	peaks <- attr(sampleset,"peaks")

	if (is.null(names(sampleset))) {
		names(sampleset) <- as.character(1:length(sampleset))
	}
	setNames <- names(sampleset)

	# Generate a data frame with the peak piece and (if applicable) window piece data
	# stored and marked by peak and type.  Also generate a frame containing all
	# local maxima.
	spectrum_df <- data.frame()
	peak_df <- data.frame()
	maxima_df <- data.frame(sample=c(),peak=c(),mz=c(),intensity=c(),type=c())
	for (j in seq_along(sampleset)) {
		cursample <- sampleset[[j]]
		for (i in seq_along(peaks)) {
			if (usePeakNames && !is.null(names(peaks)) &&  names(peaks)[[i]]!="") {
				peakName <- names(peaks)[[i]]
			} else {
				peakName <- paste0(peaks[[i]]$value," m/z")
			}
			curpiece_df <- data.frame(sample=setNames[j],
									  peak=peakName,
									  mz=cursample[[i]]$peakPiece$mz,
									  intensity=cursample[[i]]$peakPiece$intensity,
									  type="Peak")

			peak_df <- rbind(peak_df,data.frame(sample=setNames[j],
												peak=peakName,
												mz=peaks[[i]]$value,
												intensity=0))
			if (!is.null(cursample[[i]]$windowPiece)) {
				curwindow_df <- data.frame(sample=setNames[j],
										   peak=peakName,
										   mz=cursample[[i]]$windowPiece$mz,
										   intensity=cursample[[i]]$windowPiece$intensity,
										   type="Context")
				curpiece_df <- rbind(curpiece_df,curwindow_df)
			}
			spectrum_df <- rbind(spectrum_df,curpiece_df)

			if (length(cursample[[i]]$maxima$mz)>0) {
				maxima_df <- rbind(maxima_df,
								   data.frame(sample=setNames[j],
								   		   peak=peakName,
								   		   mz=cursample[[i]]$maxima$mz,
								   		   intensity=cursample[[i]]$maxima$intensity,
								   		   type="Peak"))
			}
		}
	}

	facet_df <- data.frame(sample=rep(peak_df$sample,each=length(unique(spectrum_df$type))),
						   peak=rep(peak_df$peak,each=length(unique(spectrum_df$type))),
						   type=rep(unique(spectrum_df$type,times=length(sampleset)*length(peaks))))
	facet_df$facet <- paste0(facet_df$sample,": ",facet_df$peak," (",facet_df$type,")")

	spectrum_df <- merge(spectrum_df,facet_df,by=c("sample","peak","type"))
	maxima_df <- merge(maxima_df,facet_df,by=c("sample","peak","type"))
	peak_df <- merge(peak_df,facet_df,by=c("sample","peak"))
	spectrum_df$facet <- factor(spectrum_df$facet,facet_df$facet)
	maxima_df$facet <- factor(maxima_df$facet,facet_df$facet)
	peak_df$facet <- factor(peak_df$facet,facet_df$facet)

	# Generate facetted ggplot
	if (freey) { fscale <- "free" }
	else { fscale <- "free_x" }
	gplot <- ggplot2::ggplot(spectrum_df,ggplot2::aes(x=mz,y=intensity))+
		ggplot2::geom_hline(yintercept=0,colour="black")+
		ggplot2::geom_vline(ggplot2::aes(xintercept=mz),colour="red",data=peak_df)
	if (nrow(maxima_df)>0) {
		gplot <- gplot+
			ggplot2::geom_vline(ggplot2::aes(xintercept=mz),colour="blue",linetype=2,data=maxima_df)
	}
	gplot <- gplot+
		ggplot2::geom_line()+
		ggplot2::scale_y_continuous(limits=c(0,NA),labels=function(v) sprintf("%sM",v/1e6))+
		ggplot2::labs(x="m/z",y="Intensity")+
		ggplot2::facet_wrap(ggplot2::vars(facet),ncol=length(peaks)*length(unique(facet_df$type)),scale=fscale,drop=FALSE)
	gplot
}

#' Plot an RTMS spectrum object
#'
#' `plotRtmsSpectrum()` takes an RTMS spectrum object and produces a `ggplot`
#' object depicting the spectrum
#'
#' @details
#' Unlike a sample object, an RTMS spectrum is actually quite simple; just a
#' vector of m/z values and vector of intensities.  Ordinarily, this could be
#' done using standard `ggplot2` functions, such as geom_line.  However, mass
#' spectra can often be quite large (on the order of millions of measurements),
#' and sending all that data to be plotted can be computationally intractable.
#' `plotRtmsSpectrum()` therefore selects a subset of up to 10000 m/z-intensity
#' pairs from the original spectrum to produce a representative plot without
#' rendering millions of points.  Any points that are sufficiently larger than
#' their local surroundings (including all relevant peaks) will be included in
#' this subset, as well as a random sampling of points closer to the baseline.
#' This ensures that the peaks plotted will always be present.  However, there
#' will be slight differences from one plot to the next in terms of baseline
#' points plotted.  This can be eliminated by fixing the random seed using
#' `set.seed` before plotting.
#'
#' We also strongly discourage using `xlim` or seetting the x-coordinate
#' boundaries using standard `ggplot2` methods, as these will only be applied
#' after the data has been down-sampled.  If you would like to plot a particular
#' subset of the spectrum, it is recommended that you use the `limits` parameter
#' of this function instead.
#'
#' @param spectrum An object of class `rtmsSpectrum`.
#'
#' @param limits An optional parameter to control the bounds of the m/z x axis.
#' If set to NULL, the default, the full spectrum will be plotted.  Otherwise,
#' `limits` should be a two element numeric vector; if one element is `NA`, then
#' only the other boundary will be enforced on the x-axis.
#'
#' @returns A `ggplot` object depicting the RTMS spectrum.
#'
#' @export
#' @examples
#' plot1 <- plotRtmsSpectrum(exampleSpectrum)
#' plot2 <- plotRtmsSpectrum(exampleSpectrum,limits=c(1500,1550)) +
#' 				ggplot2::geom_vline(xintercept=c(1516.83,1530.84),
#' 									colour="red",linetype=2)
plotRtmsSpectrum <- function(spectrum,limits=NULL) {
	mz <- NULL
	intensity <- NULL

	numPts <- 10000
	upperProp <- 0.125
	padding <- 2^16

	if (!inherits(spectrum,"rtmsSpectrum")) {
		stop("Parameter 'spectrum' must be an object of class 'rtmsSpectrum'.")
	}

	if (is.null(limits)) {
		rel_intensity <- spectrum$intensity
		rel_mz <- spectrum$mz
	} else {
		if (!is.numeric(limits) || length(limits)!=2) {
			stop("Parameter 'limits' must be a numeric vector of length 2.")
		}
		if (is.na(limits[1])) {
			if (is.na(limits[2])) {
				rel_intensity <- spectrum$intensity
				rel_mz <- spectrum$mz
			} else {
				rel_intensity <- spectrum$intensity[spectrum$mz<=limits[2]]
				rel_mz <- spectrum$mz[spectrum$mz<=limits[2]]
			}
		} else {
			if (is.na(limits[2])) {
				rel_intensity <- spectrum$intensity[spectrum$mz>=limits[1]]
				rel_mz <- spectrum$mz[spectrum$mz>=limits[1]]
			} else {
				if (limits[1]>=limits[2]) {
					stop("First limit value must be lower than second limit value.")
				}
				rel_intensity <- spectrum$intensity[spectrum$mz>=limits[1] & spectrum$mz<=limits[2]]
				rel_mz <- spectrum$mz[spectrum$mz>=limits[1] & spectrum$mz<=limits[2]]
			}
		}
	}
	if (length(rel_mz)==0) {
		stop("No measurements were found within the given limits.")
	}

	rel_intensity <- rel_intensity[order(rel_mz)]
	rel_mz <- rel_mz[order(rel_mz)]

	if (length(rel_mz)<(1.5*numPts)) {
		pinds <- seq_along(rel_mz)
	} else {
		mzdiff <- diff(rel_mz)
		mzdiff <- c(mzdiff[1],mzdiff)

		numOtherPts <- numPts*(1-upperProp)
		pv <- pmin(numOtherPts*mzdiff/sum(mzdiff),1)

		if (length(rel_mz)<=padding) {
			baseline <- mean(rel_intensity)
		} else {
			baseline <- localMean(rel_intensity,padding)
		}
		ratio <- rel_intensity/baseline
		threshold <- stats::quantile(ratio,1-(numPts*upperProp)/length(rel_mz))

		pinds <- which(ratio>threshold | stats::runif(length(rel_mz))<pv)
	}

	plot_mz <- rel_mz[pinds]
	plot_intensity <- rel_intensity[pinds]
	plot_df <- data.frame(mz=plot_mz,intensity=plot_intensity)
	ggplot2::ggplot(plot_df,ggplot2::aes(mz,intensity))+
		ggplot2::geom_line()
}

localMean <- function(vec,pad=2^16) {
	hpad <- floor(pad/2)
	pad <- 2*hpad
	fvec <- c(rev(vec[1:hpad]),
			  vec,
			  rev(vec)[1:hpad])
	zpad <- rep(0,pad)
	guess <- (cumsum(c(fvec,zpad))-cumsum(c(zpad,fvec)))
	guess <- guess[(pad+hpad+1):(length(vec)+pad+hpad)]
	guess/pad
}
