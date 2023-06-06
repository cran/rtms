# rtmsBrukerBAFReader.R

#' Open a Bruker single-acquisition BAF directory
#'
#' @description
#' Creates an RTMS reader object (of class `rtmsBrukerBAFReader`) which can
#' extract data from a Bruker single acquisition directory (extension ".d")
#'
#' @details
#' Currently, RTMS can create reader objects for two binary Bruker data formats,
#' BAF (presumably standing for "Bruker acquisition format") holding data from
#' a single spectrum acquisition, and MCF (probably "multiacquisition container
#' format") containing data from multiple spectra acquired in a single run.
#' Both formats hold data in a directory marked with the extension ".d". The
#' single acquisition BAF format directory contains three essential data files:
#' the main raw data file with extension ".baf", an index file with extension
#' ".baf_idx", and a calibration data file with extension ".baf_xtr". This
#' function processes the index and calibration files so that raw data can be
#' extracted quickly on demand from the ".baf" file.
#'
#' An important note: when a MCF multi-acquistion reader is created, it creates
#' an open connection to the raw data file, which allows for quicker processing
#' of many spectra in a single file.  However, because a BAF file contains only
#' a single spectrum, there is little advantage to maintaining an open
#' connection, so the connection is re-opened every time data is read.  Thus,
#' while it is important to close an MCF reader object when all data is
#' extracted, it is not necessary to close an object of class
#' `rtmsBrukerBAFReader`.
#'
#' @param bafdir A directory (usually with the extension ".d") containing data
#' from a single Bruker acquisition. This directory will contain a file with
#' extension ".baf" that holds the primary raw data, as well as an index and
#' calibration file (see Details).
#'
#' @returns An object of class `rtmsBrukerBAFReader` which can extract raw data
#' from the specified directory
#'
#' @export
newBrukerBAFReader <- function(bafdir) {
	# Locate the relevant data files in the specified Bruker directory:
	#
	# - A main data file, with extension .baf
	# - The corresponding BAF index file, with the same file name but
	#   extension .baf_idx
	# - A calibration data file, with an identical file name and extension
	#   .baf_xtr
	mainFile <- dir(bafdir,"\\.baf$")
	if (length(mainFile)==0) { stop("Main data file not found in directory.") }
	else { mainFile <- mainFile[1] }
	mainIndex <- paste0(mainFile,"_idx")
	if (!file.exists(paste0(bafdir,"/",mainIndex))) {
		stop("Main index file not found in directory.")
	}
	calibFile <- paste0(mainFile,"_xtr")
	if (!file.exists(paste0(bafdir,"/",calibFile))) {
		stop("Calibration data file not found in directory.")
	}
	files <- list(main=mainFile,index=mainIndex,calibration=calibFile)

	offsets <- readBafIndexFile(paste0(bafdir,"/",mainIndex))
	calibration <- readBafExtraFile(paste0(bafdir,"/",calibFile))
	parameters <- baf_readFileHeader(paste0(bafdir,"/",mainFile),offsets$header)

	new_rtmsBrukerBAFReader(dir=bafdir,files=files,
							calibration=calibration,offsets=offsets,
							parameters=parameters)
}

#' @param x The BAF reader object
#' @param ... Additional parameters
#'
#' @describeIn getBrukerBAFSpectrum The S3 method `getSpectrum` for objects
#' of class `rtmsBrukerBAFReader`; calls `getBrukerBAFSpectrum`
#' @export
getSpectrum.rtmsBrukerBAFReader <- function(x,...) {
	getBrukerBAFSpectrum(x,...)
}

#' @param x The BAF reader object
#' @param ... Additional parameters
#'
#' @describeIn getBrukerBAFSample The S3 method `getSample` for objects
#' of class `rtmsBrukerBAFReader`; calls `getBrukerBAFSample`
#' @export
getSample.rtmsBrukerBAFReader <- function(x,peaks,...) {
	getBrukerBAFSample(x,peaks,...)
}


#' Extract a spectrum from a Bruker BAF directory
#'
#' Extracts an RTMS spectrum object (of class `rtmsSpectrum`) from a single
#' acquisition Bruker BAF directory opened using an RTMS reader object (of
#' class `rtmsBrukerBAFReader`).  Because a BAF directory only contains one
#' spectrum, no additional parameters are needed to specify the spectrum to be
#' extracted.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerBAFReader`
#'
#' @returns An RTMS spectrum object of class `rtmsSpectrum`
#'
#' @export
getBrukerBAFSpectrum <- function(reader) {
	if (!inherits(reader,"rtmsBrukerBAFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerBAFReader'.") }

	fhigh <- reader$calibration$frequencyHigh
	fwidth <- reader$calibration$frequencyWidth
	fsize <- reader$calibration$size
	alpha <- reader$calibration$alpha
	beta <- reader$calibration$beta

	massToIndex <- function(p) { return(fsize*(fwidth-(fhigh/(p-alpha/fhigh) + beta))/fwidth) }
	indexToMass <- function(i) { return(fhigh/((fwidth*(fsize-i)/fsize)-beta) + alpha/fhigh) }

	fcon <- file(paste0(reader$dir,"/",reader$files$main),"rb")
	on.exit(close(fcon),add=TRUE)

	baf_goToOffset(fcon,reader$offsets$spectrum)
	sblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	sdelim <- readBin(fcon,"raw",4)
	sheaderSize <- readBin(fcon,"integer",1,size=4,endian="little")

	spectrumLength <- (sblockSize-sheaderSize)/4
	seek(fcon,sheaderSize-12,"current")
	spectrum <- readBin(fcon,"double",spectrumLength,size=4,endian="little")

	rtmsSpectrum(indexToMass(0:(spectrumLength-1)),spectrum)
}

#' Extract a sample from a Bruker BAF directory
#'
#' Extracts an RTMS sample object (of class `rtmsSample`) from a single
#' acquisition Bruker BAF directory opened using an RTMS reader object (of
#' class `rtmsBrukerBAFReader`).  Because a BAF directory only contains one
#' spectrum, no additional parameters are needed to specify the spectrum from
#' which to take the sample.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerBAFReader`
#' @param peaks A list of peak objects of class `rtmsPeak`
#'
#' @returns An RTMS sample object of class `rtmsSample`
#'
#' @export
getBrukerBAFSample <- function(reader,peaks) {
	if (!inherits(reader,"rtmsBrukerBAFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerBAFReader'.") }

	fhigh <- reader$calibration$frequencyHigh
	fwidth <- reader$calibration$frequencyWidth
	fsize <- reader$calibration$size
	alpha <- reader$calibration$alpha
	beta <- reader$calibration$beta

	massToIndex <- function(p) { return(fsize*(fwidth-(fhigh/(p-alpha/fhigh) + beta))/fwidth) }
	indexToMass <- function(i) { return(fhigh/((fwidth*(fsize-i)/fsize)-beta) + alpha/fhigh) }

	fcon <- file(paste0(reader$dir,"/",reader$files$main),"rb")
	on.exit(close(fcon),add=TRUE)

	# Move to extracted peak blob.  The sum of peaks within the
	# transformed index window is the current assay measure
	baf_goToOffset(fcon,reader$offsets$maxima)
	pblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	pdelim <- readBin(fcon,"raw",4)
	pheaderSize <- readBin(fcon,"integer",1,size=4,endian="little")

	seek(fcon,pheaderSize-12,"current")
	numPeaks <- readBin(fcon,"integer",1,size=4,endian="little")
	peakOffset <- readBin(fcon,"integer",1,size=4,endian="little")
	seek(fcon,peakOffset-8,"current")
	peakIndexes <- readBin(fcon,"double",numPeaks,size=8,endian="little")
	peakMasses <- indexToMass(peakIndexes)
	peakIntensities <- readBin(fcon,"double",numPeaks,size=8,endian="little")

	# Move to raw data blob.  From here we will extract all measured spectra
	# values within a given width of the desired peaks, and take their average
	baf_goToOffset(fcon,reader$offsets$spectrum)
	sblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	sdelim <- readBin(fcon,"raw",4)
	sheaderSize <- readBin(fcon,"integer",1,size=4,endian="little")

	numValues <- (sblockSize-sheaderSize)/4
	seek(fcon,sheaderSize-12,"current")

	# For each desired peak, extract the relevant peak piece, maxima and, (if
	# applicable) window piece of the overall spectrum. This can be read directly
	# from the file
	subsamples <- list()
	for (i in seq_along(peaks)) {
		currentPeak <- peaks[[i]]
		indexBounds <- massToIndex(currentPeak$bounds)
		indexBoundsD <- shrinkBounds(indexBounds)+1
		if (indexBoundsD[1]>numValues || indexBoundsD[2]<1) {
			peakPiece <- rtmsSpectrumPiece(c(),c())
		} else {
			indexboundsD <- pmin(pmax(indexBoundsD,1),numValues)
			currIndexes <- indexBoundsD[1]:indexBoundsD[2]
			seek(fcon,4*(indexBoundsD[1]-1),origin="current")
			currValues <- readBin(fcon,"double",diff(indexBoundsD)+1,size=4,endian="little")
			seek(fcon,-4*indexBoundsD[2],origin="current")
			peakPiece <- rtmsSpectrumPiece(indexToMass(currIndexes-1),currValues)
		}
		if (!is.null(currentPeak$window)) {
			indexBounds <- massToIndex(currentPeak$window)
			indexBoundsD <- shrinkBounds(indexBounds)+1
			if (indexBoundsD[1]>numValues || indexBoundsD[2]<1) {
				windowPiece <- rtmsSpectrumPiece(c(),c())
			} else {
				indexboundsD <- pmin(pmax(indexBoundsD,1),numValues)
				currIndexes <- indexBoundsD[1]:indexBoundsD[2]
				seek(fcon,4*(indexBoundsD[1]-1),origin="current")
				currValues <- readBin(fcon,"double",diff(indexBoundsD)+1,size=4,endian="little")
				seek(fcon,-4*indexBoundsD[2],origin="current")
				windowPiece <- rtmsSpectrumPiece(indexToMass(currIndexes-1),currValues)
			}
		} else { windowPiece <- NULL }

		relmax <- which(peakMasses>=currentPeak$bounds[1] & peakMasses<=currentPeak$bounds[2])
		maxima <- rtmsSpectrumPiece(peakMasses[relmax],peakIntensities[relmax])
		subsamples[[i]] <- new_rtmsSubsample(peakPiece,maxima,windowPiece)
	}

	new_rtmsSample(subsamples,peaks)
}

#' Retrieve specific metadata values from a Bruker BAF file
#'
#' Retrieves a list of specific metadata values (including instrument data,
#' acquisition parameters, processing and analysis directives, etc.) from a
#' Bruker single acquisition BAF directory (represented by an
#' `rtmsBrukerBAFReader` object).
#'
#' @param reader An RTMS reader object of class `rtmsBrukerBAFReader`
#'
#' @param names A character vector of metadata names
#'
#' @returns A named list of values corresponding to the metadata values
#' specified.  All values will be returned as a string, including numeric
#' quantities (with units if appropriate).
#'
#' @export
getBrukerBAFMetadata <- function(reader,names) {
	metadata <- getBrukerBAFAllMetadata(reader)

	mtab <- merge(data.frame(displayName=names),metadata)
	output <- list()
	for (iter in seq_len(nrow(mtab))) {
		output[[mtab$displayName[[iter]]]] <- mtab$stringValue[[iter]]
	}
	output
}

#' Retrieve all metadata values from a Bruker BAF file
#'
#' Retrieves a table of all metadata values (including instrument data,
#' acquisition parameters, processing and analysis directives, etc.) from a
#' Bruker single acquisition BAF directory (represented by an
#' `rtmsBrukerBAFReader` object).
#'
#' @param reader An RTMS reader object of class `rtmsBrukerBAFReader`
#'
#' @returns A data frame with all metadata parameters for the acquisition. The
#' data frame will have five columns: `rowIndex`, a numeric index for each
#' metadata value; `parameterName`, the internal identifier of the parameter in
#' Bruker software systems; `parameterGroup`, the group of parameters that each
#' value belongs to; `displayName`, the string used to specify the parameter to
#' users (i.e. how the parameter would be labelled in a user interface); and
#' `stringValue`, a character column containing the value of each metadata
#' parameter.  Numeric quantities will also be returned as strings, with units
#' if appropriate.
#'
#' @export
getBrukerBAFAllMetadata <- function(reader) {
	fcon <- file(paste0(reader$dir,"/",reader$files$main),"rb")
	on.exit(close(fcon),add=TRUE)

	baf_goToOffset(fcon,reader$offsets$metadata)
	mblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	mdelim <- readBin(fcon,"raw",4)
	mheaderSize <- readBin(fcon,"integer",1,size=4,endian="little")
	seek(fcon,mheaderSize-12,origin="current")

	nvblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	nvdelim <- readBin(fcon,"raw",4)
	nvheaderSize <- readBin(fcon,"integer",1,size=4,endian="little")
	nnumValues <- readBin(fcon,"integer",1,size=4,endian="little")

	ntable <- data.frame(rowIndex=rep(NA,nnumValues),numValue=rep(NA,nnumValues))
	for (niter in seq_len(nnumValues)) {
		ntable$rowIndex[[niter]] <- readBin(fcon,"integer",1,size=2,endian="little",signed=FALSE)
		ntable$numValue[[niter]] <- readBin(fcon,"double",1,size=8,endian="little")
	}
	ntable <- merge(ntable,reader$parameters,by="rowIndex",all.x=TRUE)
	ntable$stringValue <- ""
	for (niter in seq_len(nnumValues)) {
		if (is.na(ntable$values[[niter]])) {
			if (ntable$format[[niter]]=="%d") {
				ntable$stringValue[[niter]] <- sprintf(ntable$format[[niter]],round(ntable$numValue[[niter]]))
			} else {
				ntable$stringValue[[niter]] <- sprintf(ntable$format[[niter]],ntable$numValue[[niter]])
			}
		} else {
			ntable$stringValue[[niter]] <- baf_recode(round(ntable$numValue[[niter]]),ntable$values[[niter]])
		}
		if (!is.na(ntable$unit[[niter]])) {
			ntable$stringValue[[niter]] <- paste(ntable$stringValue[[niter]],ntable$unit[[niter]])
		}
	}
	ntable <- ntable[,c("rowIndex","parameterName","parameterGroup","displayName","stringValue")]

	svblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	svdelim <- readBin(fcon,"raw",4)
	svheaderSize <- readBin(fcon,"integer",1,size=4,endian="little")
	snumValues <- readBin(fcon,"integer",1,size=4,endian="little")

	stable <- data.frame(rowIndex=rep(NA,snumValues),stringValue=rep(NA,snumValues))
	for (siter in seq_len(snumValues)) {
		curRow <- readBin(fcon,"raw",2)
		rowSize <- readBin(curRow,"integer",1,size=2,endian="little",signed=FALSE)
		curRow <- c(curRow,readBin(fcon,"raw",rowSize-2))
		stable$rowIndex[[siter]] <- readBin(curRow[3:4],"integer",1,size=2,endian="little",signed=FALSE)
		stringLength <- readBin(curRow[5:8],"integer",1,size=4,endian="little")-1
		stable$stringValue[[siter]] <- rawToChar(curRow[9:(8+stringLength)])
	}
	stable <- merge(stable,reader$parameters,by="rowIndex",all.x=TRUE)
	stable <- stable[,c("rowIndex","parameterName","parameterGroup","displayName","stringValue")]

	mtable <- rbind(ntable,stable)
	mtable <- mtable[order(mtable$rowIndex),]
	return(mtable)
}

new_rtmsBrukerBAFReader <- function(dir,files,calibration,offsets,parameters) {
	stopifnot(is.character(dir))
	stopifnot(is.list(files))
	stopifnot(is.list(calibration))
	stopifnot(is.list(offsets))
	stopifnot(is.data.frame(parameters))

	structure(list(dir=dir, files=files, calibration=calibration,
				   offsets=offsets, parameters=parameters),
			  class="rtmsBrukerBAFReader")
}

readBafIndexFile <- function(filename) {
	rawVectors <- bafidx_getRawVectors(filename)

	outdf <- data.frame()
	for (rv in rawVectors) {
		if (rv[2]==as.raw(0)) { next }
		outdf <- rbind(outdf,data.frame(idA=readBin(rv[1],"integer",1,size=1,signed=FALSE),
										idB=readBin(rv[2],"integer",1,size=1,signed=FALSE),
										flag=readBin(rv[9:12],"integer",1,size=4,endian="little"),
										offset=readBin(rv[17:20],"integer",1,size=4,endian="little"),
										offsetPage=readBin(rv[21:24],"integer",1,size=4,endian="little")))
	}

	offsetLookup <- data.frame(name=c("header","metadata","spectrum","maxima"),
							   idA=c(1,2,1,2),idB=c(48,32,16,16))
	offsetTable <- merge(offsetLookup,outdf,by=c("idA","idB"),all.x=TRUE)
	offsets <- list()
	for (i in seq_len(nrow(offsetTable))) {
		offsets[[offsetTable$name[[i]]]] <- baf_reshiftOffset(c(offsetTable$offset[[i]],offsetTable$offsetPage[[i]]))
	}
	offsets
}
readBafExtraFile <- function(filename) {
	rawVectors <- bafidx_getRawVectors(filename)

	for (rv in rawVectors) {
		if (rv[1]!=as.raw(4) || rv[2]!=as.raw(1)) { next }
		calrow <- list()
		offset <- 4
		headerSize <- readBin(rv[(offset+1):(offset+4)],"integer",1,size=4,endian="little")
		offset <- offset+headerSize
		blockSize <- readBin(rv[(offset+1):(offset+4)],"integer",1,size=4,endian="little")
		offset <- offset+blockSize-headerSize
		calibMode <- readBin(rv[(offset+9):(offset+12)],"integer",1,size=4,endian="little")

		if (calibMode==5) { calrow$alpha <- readBin(rv[(offset+13):(offset+20)],"double",1,size=8,endian="little") }
		else { calrow$alpha <- 0 }
		calrow$beta <- readBin(rv[(offset+21):(offset+28)],"double",1,size=8,endian="little")
		calrow$frequencyHigh <- readBin(rv[(offset+29):(offset+36)],"double",1,size=8,endian="little")
		if (calibMode==4) { calrow$beta <- -calrow$beta }

		calrow$size <- readBin(rv[(offset+41):(offset+44)],"integer",1,size=4,endian="little")

		calrow$frequencyLow <- readBin(rv[(offset+49):(offset+56)],"double",1,size=8,endian="little")
		calrow$frequencyWidth <- readBin(rv[(offset+57):(offset+64)],"double",1,size=8,endian="little")
		return(calrow)
	}

	return(NULL)
}

bafidx_getRawVectors <- function(filename) {
	fcon <- file(filename,"rb")
	on.exit(close(fcon),add=TRUE)

	rawvecs <- list()
	index <- 1
	while (index>0) {
		temp <- readBin(fcon,"integer",1,size=4,endian="little")
		delim <- readBin(fcon,"raw",4)
		blockSize <- readBin(fcon,"integer",1,size=4,endian="little")
		if (blockSize==0) {
			index = 0
		} else {
			rawvecs[[index]] <- readBin(fcon,"raw",blockSize-24)
			closer <- readBin(fcon,"raw",12)
			index <- index+1
		}
	}

	rawvecs
}

baf_reshiftOffset <- function(offset) {
	if (length(offset)==1) { offset <- c(offset,0) }
	offset[2] <- offset[2]*2
	if (offset[1]<0) {
		offset[2] <- offset[2]+1
		offset[1] <- offset[1]+(2^31)
	}
	offset
}
baf_goToOffset <- function(fcon,offset) {
	seek(fcon,0)
	if (length(offset)>1 && offset[2]>0) {
		for (piter in seq_len(offset[2])) { seek(fcon,2^31,origin="current") }
	}
	seek(fcon,offset[1],origin="current")
}
baf_readFileHeader <- function(filename,headerOffset) {
	fcon <- file(filename,"rb")
	on.exit(close(fcon),add=TRUE)

	baf_goToOffset(fcon,headerOffset)
	blockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	delim <- readBin(fcon,"raw",4)
	subdelim <- readBin(fcon,"raw",4)
	headerSize <- readBin(fcon,"integer",1,size=4,endian="little")
	nextOffset <- readBin(fcon,"integer",1,size=4,endian="little")
	seek(fcon,nextOffset-20,origin="current")
	mblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	seek(fcon,mblockSize-4,origin="current")
	hblockSize <- readBin(fcon,"integer",1,size=4,endian="little")
	hdelim <- readBin(fcon,"raw",4)
	hheaderSize <- readBin(fcon,"integer",1,size=2,endian="little",signed=FALSE)
	hnumRows <- readBin(fcon,"integer",1,size=2,endian="little",signed=FALSE)

	headerRows <- data.frame()
	for (hni in 1:hnumRows) {
		nextRow <- baf_readHeaderRow(fcon)
		headerRows <- rbind(headerRows,nextRow)
	}
	return(headerRows)
}
baf_readHeaderRow <- function(fcon) {
	curRow <- readBin(fcon,"raw",2)
	rowSize <- readBin(curRow,"integer",1,size=2,endian="little",signed=FALSE)
	curRow <- c(curRow,readBin(fcon,"raw",rowSize-2))
	rowIndex <- readBin(curRow[3:4],"integer",1,size=2,endian="little",signed=FALSE)
	rowType <- readBin(curRow[5:8],"integer",1,size=4,endian="little")

	fieldOffsets <- readBin(curRow[9:20],"integer",6,size=2,endian="little",signed=FALSE)
	fieldLengths <- c(fieldOffsets,rowSize)
	for (foi in 6:1) {
		if (fieldLengths[foi]==0) { fieldLengths[foi] <- fieldLengths[foi+1] }
	}
	fieldLengths <- diff(fieldLengths)
	fields <- rep(NA,6)
	for (foi in 1:6) {
		if (fieldOffsets[foi]==0) { next }
		fields[foi] <- rawToChar(curRow[(fieldOffsets[foi]+1):(fieldOffsets[foi]+fieldLengths[foi])])
	}

	return(list(rowIndex=rowIndex,rowType=rowType,parameterName=fields[1],parameterGroup=fields[2],
				displayName=fields[3],values=fields[4],format=fields[5],unit=fields[6]))
}
baf_recode <- function(value,valueList) {
	valuePairs <- strsplit(valueList,";")[[1]]
	valueTable <- as.data.frame(do.call(rbind,strsplit(valuePairs,":")))
	names(valueTable) <- c("code","stringValue")

	return(valueTable$stringValue[[which(valueTable$code==value)[[1]]]])
}
