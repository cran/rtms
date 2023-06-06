# rtmsBrukerMCFReader.R

#' Open a Bruker multi-acquisition MCF directory
#'
#' @description
#' Creates an RTMS reader object (of class `rtmsBrukerMCFReader`) which can
#' extract data from a Bruker multi-acquisition directory (extension ".d")
#'
#' @details
#' Currently, RTMS can create reader objects for two binary Bruker data formats,
#' BAF (presumably standing for "Bruker acquisition format") holding data from
#' a single spectrum acquisition, and MCF (probably "multi-acquisition container
#' format") containing data from multiple spectra acquired in a single run.
#' Both formats hold data in a directory marked with the extension ".d". The
#' multi-acquisition MCF format directory contains four essential data files:
#' the main raw data file ending in "_1" with extension ".mcf", a matching
#' index file with extension ".mcf_idx", and a calibration data file ending in
#' "_2" with  extension ".mcf", and the matching calibration index file with
#' extension ".mcf_idx". This function preprocesses the main index file and
#' calibration files so that raw data can be extracted quickly on demand from
#' the main ".mcf" file.
#'
#' An important note: when a MCF multi-acquisition reader is created, it creates
#' an open connection to the raw data file, which allows for quicker processing
#' of many spectra in a single file.  This connection will remain open until
#' the reader is closed with the `close` function.
#'
#' @param mcfdir A directory (usually with the extension ".d") containing data
#' from a Bruker multi-acquisition run. This directory will contain a files with
#' extenssion ".mcf" and matching index files with extension ".mcf_idx" (see
#' Details).
#'
#' @returns A reader object of class `rtmsBrukerMCFReader` with an open
#' connection to the main ".mcf" data file.
#'
#' @export
newBrukerMCFReader <- function(mcfdir) {
	# Locate the relevant data files in the specified Bruker directory:
	#
	# - A main data file, with filename ending in "_1" and extension .mcf
	# - The corresponding SQLite index file, with the same file name but
	#   extension .mcf_idx
	# - A calibration datafile, with an identical file name to the main
	#   data file but ending in "_2", with extension .mcf
	# - The SQLite calibration data index file, with extension .mcf_idx
	mainFile <- dir(mcfdir,"_1.mcf$")
	if (length(mainFile)==0) { stop("Main data file not found in directory.") }
	else { mainFile <- mainFile[1] }
	mainIndex <- paste0(mainFile,"_idx")
	if (!file.exists(paste0(mcfdir,"/",mainIndex))) {
		stop("Main index file not found in directory.")
	}
	calibFile <- sub("_1\\.mcf","_2.mcf",mainFile)
	if (!file.exists(paste0(mcfdir,"/",calibFile))) {
		stop("Calibration data file not found in directory.")
	}
	calibIndex <- paste0(calibFile,"_idx")
	if (!file.exists(paste0(mcfdir,"/",calibIndex))) {
		stop("Calibration index file not found in directory.")
	}
	files <- list(main=mainFile,index=mainIndex,calibration=calibFile,calIndex=calibIndex)

	# Initialize the processing of the multispectrum file by reading in the
	# location of binary data blocks in the main data file and calibration
	# datafile from their corresponding SQLite index files.
	offsetTable <- readBrukerMCFIndexFile(paste0(mcfdir,"/",mainIndex),calibration=FALSE)
	calOffsets <- readBrukerMCFIndexFile(paste0(mcfdir,"/",calibIndex),calibration=TRUE)
	calOffsets$index <- 1:nrow(calOffsets)

	# Locate the relevant calibration data blocks in the calibration data
	# file for each spot measured, and store the calibration data
	spotCalOffsets <- merge(calOffsets,offsetTable[offsetTable$BlobResType==258,"id",drop=FALSE],by.x="toId",by.y="id",sort=FALSE)
	spotCalOffsets <- spotCalOffsets[order(spotCalOffsets$index),]
	spotCalData <- readMCFCalibration(paste0(mcfdir,"/",calibFile),spotCalOffsets)
	spotCalData$Index <- 1:nrow(spotCalData)

	# Use spot calibration data as a spot table
	spotTable <- spotCalData[order(spotCalData$Index),]
	names(spotTable)[names(spotTable)=="toId"] <- "id"

	# Open the main data file
	mcfcon <- file(paste0(mcfdir,"/",mainFile),"rb")

	# Retrieve a full representation of the method metadata, including
	# parameter keys, and lookup table values
	paramBlobOffset <- offsetTable$Offset[3]
	metadata <- retrieveMCFMetadata(mcfcon,paramBlobOffset,0)

	# Returned value is a object with the following fields:
	# dir : ".d" Bruker data directory
	# files : Relevant data files
	# spotTable : Table of spots with well ids, calibration parameters
	# offsetTable : Table of binary offsets in main data file
	# metadata : A two-element list containing metadata declarations and values
	# con : Open connection to main Bruker mcf data file
	return(new_rtmsBrukerMCFReader(dir=mcfdir,
								   files=files,
								   spotTable=spotTable,
								   offsetTable=offsetTable,
								   metadata=metadata,
								   con=mcfcon))
}

#' Reopen a reader object
#'
#' An S3 generic function for reopening reader objects that have been close
#'
#' @param x The object to be re-opened
#'
#' @param ... Other possible arguments for specific object types
#'
#' @returns An object of the same type as `x`
#'
#' @export
reopen <- function(x,...) {
	UseMethod("reopen")
}
#' @export
reopen.default <- function(x,...) {
	invisible(x)
}
#' @param x The reader object to reopen
#' @param ... Included for S3 compatibility
#'
#' @describeIn closeBrukerMCFReader The S3 method `close` for objects of class
#' `rtmsBrukerMCFReader`; calls `closeBrukerMCFReader`
#' @export
reopen.rtmsBrukerMCFReader <- function(x,...) {
	reopenBrukerMCFReader(x)
}
#' @param con The reader object to be closed
#'
#' @describeIn closeBrukerMCFReader The S3 method `reopen` for objects of class
#' `rtmsBrukerMCFReader`; calls `reopenBrukerMCFReader`
#' @export
close.rtmsBrukerMCFReader <- function(con,...) {
	closeBrukerMCFReader(con)
}

#' @param x The MCF reader object
#' @param ... Additional parameters
#'
#' @describeIn getBrukerMCFSpectrum The S3 method `getSpectrum` for objects
#' of class `rtmsBrukerMCFReader`; calls `getBrukerMCFSpectrum`
#' @export
getSpectrum.rtmsBrukerMCFReader <- function(x,...) {
	getBrukerMCFSpectrum(x,...)
}
#' @param x The MCF reader object
#' @param ... Additional parameters
#'
#' @describeIn getBrukerMCFSample The S3 method `getSample` for objects
#' of class `rtmsBrukerMCFReader`; calls `getBrukerMCFSample`
#' @export
getSample.rtmsBrukerMCFReader <- function(x,peaks,...) {
	getBrukerMCFSample(x,peaks,...)
}
#' @param x The MCF reader object
#' @param ... Additional parameters
#'
#' @describeIn getBrukerMCFSampleSet The S3 method `getSample` for objects
#' of class `rtmsBrukerMCFReader`; calls `getBrukerMCFSampleSet`
#' @export
getSampleSet.rtmsBrukerMCFReader <- function(x,peaks,...) {
	getBrukerMCFSampleSet(x,peaks,...)
}

#' Manage an MCF reader file connection
#'
#' Closes the open connection to the main data file in a Bruker MCF reader
#' object.
#'
#' @details
#' Because Bruker MCF directories contain a potentially large number of
#' spectra, reopening a connection to the main data file when reading many
#' spectra or samples from it is inefficient and slow, especially if the file
#' is being accessed over a network connection.  The `rtmsBrukerMCFReader`
#' object therefore maintains an open connection to the main binary data file
#' until it is closed by the user. Of course, the reader object still maintains
#' all the index and calibration data, making it possible to reopen a
#' connection to the MCF directory without all the preprocessing required when
#' first opening.  Unfortunately, taking advantage of this fact is a little
#' tricky.
#'
#' In most cases, R is a functional language, with limited side effects on R
#' objects; so it is difficult to alter the state of reader object without
#' returning it explicitly.  However, one of the few cases where side-effects
#' are quite important is R's management of open file connections, with
#' functions like `close`.  Thus, calling `closeBrukerMCFReader` (and the S3
#' function `close` which wraps it) will in fact close the connection, but will
#' not return an altered copy of the reader object reflecting that it is
#' closed.  So if the user wishes to close a reader object with the possibility
#' of reopening it, they must close the reader AND assign the returned reader
#' object to the relevant name.  This will store the fact the connection has
#' been closed, and allow the `reopen` function to operate correctly.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @returns The same reader object with a closed connection
#'
#' @export
closeBrukerMCFReader <- function(reader) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	if (!is.null(reader$con)) {
		close(reader$con)
		reader$con <- NULL
	}
	invisible(reader)
}
#' @describeIn closeBrukerMCFReader Reopens a file connection to the main
#' binary data file in a Bruker MCF directgry so that data can be extracted
#' @export
reopenBrukerMCFReader <- function(reader) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	if (is.null(reader$con) || !inherits(reader$con,"file")) {
		reader$con <- file(paste0(reader$dir,"/",reader$files$main),"rb")
	}
	invisible(reader)
}

#' Get spot names and indices from a Bruker MCF file
#'
#' Assembles a table of all acquisitions in a Bruker MCF file; Bruker
#' measurements are often identified by the metadata parameter "Spot Number",
#' so this function extracts that specific metadata value and joins it with the
#' indices used to pick out spectra in other functions.  Also retrieves the
#' timestamp at which acquisition was taken, if acquisitions must be identified
#' by order.
#'
#' @param reader An openRTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @returns A data.frame with an `Index` column containing the indices of each
#' acquisition (used by other functions such as `getSpectrum` or `getSample`),
#' a `SpotNumber` column containing the "Spot Number" metadata value for each
#' acquisition, and a `Timestampe` column containing the time at which each
#' acquisition was collected by the instrument.
#'
#' @export
getBrukerMCFSpots <- function(reader) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	mcfdir <- reader$dir
	mainIndex <- reader$files$index

	# Open index file for decoding
	scon <- file(paste0(mcfdir,"/",mainIndex),"rb")
	on.exit(close(scon),add=TRUE)

	# Locate all B-tree table leaf pages for the main table (which lists the
	# titles, column names, and root pages of all tables and index in the
	# SQLite file.
	mainHandler <- function(values,cellIndex) {
		return(list(name=values[[2]],value=values[[4]]))
	}
	mainTable <- sqlt_parseBTreeTable(scon,1,mainHandler,TRUE)$output

	metaRoot <- as.numeric(mainTable$value[which(mainTable$name=="MetaDataString")])
	metastringHandler <- function(values,cellIndex) {
		valueList <- values
		names(valueList) <- c("GuidA","GuidB","MetadataId","Text")
		# Guid values are 8-byte integers in the file, and hence return two integers
		# each, which must be pasted together to be a single column in a data table
		valueList$id <- paste(c(valueList$GuidA,valueList$GuidB),collapse=" ")
		valueList$GuidA <- NULL
		valueList$GuidB <- NULL
		return(as.data.frame(valueList))
	}
	metasout <- sqlt_parseBTreeTable(scon,metaRoot,metastringHandler,first=FALSE)$output

	wmeta1 <- metasout[which(metasout$MetadataId==64),c("id","Text"),drop=FALSE]
	names(wmeta1)[[2]] <- "SpotNumber"
	wmeta2 <- metasout[which(metasout$MetadataId==34),c("id","Text"),drop=FALSE]
	names(wmeta2)[[2]] <- "Timestamp"
	wmeta <- merge(wmeta1,wmeta2,by="id",all=TRUE)

	spots <- merge(reader$spotTable[c("id","Index")],wmeta,by="id",all.x=TRUE)
	spots <- spots[order(spots$Index),-1,drop=FALSE]
	spots
}

#' @describeIn getBrukerMCFSpots Retrieves a vector of all the indices
#' (beginning with zero) of the acquisitions in the MCF file.  Faster than
#' `getBrukerMCFSpots` but contains no metadata or spot names
#' @export
getBrukerMCFIndices <- function(reader) {
	return(reader$spotTable$Index)
}

#' Extract a spectrum from a Bruker MCF directory
#'
#' Extracts an RTMS spectrum object (of class `rtmsSpectrum`) from a multi-
#' acquisition Bruker MCF directory opened using an RTMS reader object (of
#' class `rtmsBrukerMCFReader`).  A numeric index is used to identify which
#' spectrum should be extracted.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @param index A single numeric index specifying which acquisition should be
#' extracted
#'
#' @returns An RTMS spectrum object of class `rtmsSpectrum`
#'
#' @export
getBrukerMCFSpectrum <- function(reader,index) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	if (index <=0 || index>nrow(reader$spotTable)) {
		stop("Parameter 'index' is out of bounds.")
	}

	fcon <- reader$con
	spotId <- reader$spotTable$id[index]
	fhigh <- reader$spotTable$frequencyHigh[index]
	fwidth <- reader$spotTable$frequencyWidth[index]
	fsize <- reader$spotTable$size[index]
	alpha <- reader$spotTable$alpha[index]
	beta <- reader$spotTable$beta[index]

	massToIndex <- function(p) { return(fsize*(fwidth-(fhigh/(p-alpha/fhigh) + beta))/fwidth) }
	indexToMass <- function(i) { return(fhigh/((fwidth*(fsize-i)/fsize)-beta) + alpha/fhigh) }

	# Move to raw data blob.  From here we will extract all measured spectra
	# values within a given width of the desired peaks, and take their average
	rawdataBlob <- which(reader$offsetTable$id==spotId & reader$offsetTable$BlobResType==258)
	mcf_blobSeek(fcon,reader$offsetTable$Offset[rawdataBlob],reader$offsetTable$OffsetPage[rawdataBlob])

	mcf_checkBlobCodeWord(fcon)
	bin_checkBytes(fcon,charToRaw("\x01\x01"))
	seek(fcon,16,origin="current")
	numNamedEls <- bin_readVarInt(fcon)

	name <- mcf_readNamedName(fcon)
	if (name!="Intensities") { stop("First named element of raw data blob should be 'Intensities'") }
	# Intensities should be an array of 4-byte floats
	bin_checkBytes(fcon,as.raw(3),"Raw spectra must be an array.")
	typeByte <- readBin(fcon,"raw",1)
	if (typeByte==as.raw(32)) { # x20 indicates floats are uncompressed
		compressedSpectrum <- FALSE
		numValues <- bin_readVarInt(fcon)
		spectrum <- readBin(fcon,"double",numValues,size=4,endian="little")
	} else if (typeByte==as.raw(34)) { #x22 indicates floats have been compressed using gzip
		compressedSpectrum <- TRUE
		numBytes <- bin_readVarInt(fcon)
		gzippedBytes <- readBin(fcon,"raw",numBytes)
		unzippedBytes <- memDecompress(gzippedBytes,type="gzip")
		numValues <- length(unzippedBytes)/4
		spectrum <- readBin(unzippedBytes,"double",numValues,size=4,endian="little")
	} else { stop("Raw spectra must be an array of 32-bit floats.") }

	rtmsSpectrum(indexToMass(0:(numValues-1)),spectrum)
}

#' Extract a sample from a Bruker MCF directory
#'
#' Extracts an RTMS sample object (of class `rtmsSample`) from a multi-
#' acquisition Bruker MCF directory opened using an RTMS reader object (of
#' class `rtmsBrukerMCFReader`). A numeric index is used to identify which
#' acquisition the sample should be extracted from.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @param peaks A list of peak objects of class `rtmsPeak`
#' @param index A single numeric index specifying which acquisition the sample
#' set should be extracted from
#'
#' @returns An RTMS sample object of class `rtmsSample`
#'
#' @export
getBrukerMCFSample <- function(reader,peaks,index) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	if (!all(vapply(peaks,function(p) inherits(p,"rtmsPeak"),FALSE))) {
		stop("Parameter 'peaks' must be a list of objects of class 'rtmsPeak'.")
	}
	if (index <=0 || index>nrow(reader$spotTable)) {
		stop("Parameter 'index' is out of bounds.")
	}

	fcon <- reader$con
	spotId <- reader$spotTable$id[index]
	fhigh <- reader$spotTable$frequencyHigh[index]
	fwidth <- reader$spotTable$frequencyWidth[index]
	fsize <- reader$spotTable$size[index]
	alpha <- reader$spotTable$alpha[index]
	beta <- reader$spotTable$beta[index]

	massToIndex <- function(p) { return(fsize*(fwidth-(fhigh/(p-alpha/fhigh) + beta))/fwidth) }
	indexToMass <- function(i) { return(fhigh/((fwidth*(fsize-i)/fsize)-beta) + alpha/fhigh) }

	# Move to extracted peak blob.  The sum of peaks within the
	# tranformed index window is the current assay measure
	peakdataBlob <- which(reader$offsetTable$id==spotId & reader$offsetTable$BlobResType==257)
	mcf_blobSeek(fcon,reader$offsetTable$Offset[peakdataBlob],reader$offsetTable$OffsetPage[peakdataBlob])

	mcf_checkBlobCodeWord(fcon)
	bin_checkBytes(fcon,charToRaw("\x01\x01"))
	seek(fcon,16,origin="current")
	numNamedEls <- bin_readVarInt(fcon)
	seek(fcon,25,origin="current") # Skip over "PeakFinderType" element

	name <- mcf_readNamedName(fcon)
	if (name!="Indexes") { stop("Second named element of peak data blob should be 'Indexes") }
	bin_checkBytes(fcon,charToRaw("\x03\x1F")) # Indexes should be an array of 8-byte floats
	numIndexes <- bin_readVarInt(fcon)
	peakIndexes <- readBin(fcon,"double",numIndexes,size=8,endian="little")
	peakMasses <- indexToMass(peakIndexes)

	name <- mcf_readNamedName(fcon)
	if (name!="Intensities") { stop("Third named element of peak data blob should be 'Intensities") }
	bin_checkBytes(fcon,charToRaw("\x03\x20")) # Intensities should be an array of 4-byte floats
	numIntensities <- bin_readVarInt(fcon)
	if (numIndexes!=numIntensities) { stop("Number of peak indicies should match number of peak intensities") }
	peakIntensities <- readBin(fcon,"double",numIndexes,size=4,endian="little")

	# Move to raw data blob.  From here we will extract all measured spectra
	# values within a given width of the desired peaks, and take their average
	rawdataBlob <- which(reader$offsetTable$id==spotId & reader$offsetTable$BlobResType==258)
	mcf_blobSeek(fcon,reader$offsetTable$Offset[rawdataBlob],reader$offsetTable$OffsetPage[rawdataBlob])

	mcf_checkBlobCodeWord(fcon)
	bin_checkBytes(fcon,charToRaw("\x01\x01"))
	seek(fcon,16,origin="current")
	numNamedEls <- bin_readVarInt(fcon)

	name <- mcf_readNamedName(fcon)
	if (name!="Intensities") { stop("First named element of raw data blob should be 'Intensities'") }
	# Intensities should be an array of 4-byte floats
	bin_checkBytes(fcon,as.raw(3),"Raw spectra must be an array.")
	typeByte <- readBin(fcon,"raw",1)
	if (typeByte==as.raw(32)) { # x20 indicates floats are uncompressed
		compressedSpectrum <- FALSE
		numValues <- bin_readVarInt(fcon)
	} else if (typeByte==as.raw(34)) { #x22 indicates floats have been compressed using gzip
		compressedSpectrum <- TRUE
		numBytes <- bin_readVarInt(fcon)
		gzippedBytes <- readBin(fcon,"raw",numBytes)
		unzippedBytes <- memDecompress(gzippedBytes,type="gzip")
		numValues <- length(unzippedBytes)/4
		spectrum <- readBin(unzippedBytes,"double",numValues,size=4,endian="little")
	} else { stop("Raw spectra must be an array of 32-bit floats.") }

	# For each desired peak, extract the relevant peak piece, maxima and, (if
	# applicable) window piece of the overall spectrum. This can be read directly
	# from the file if the data is uncompressed; otherwise it can be taken from
	# the extracted and decompressed spectrum.
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
			if (compressedSpectrum) { currValues <- spectrum[currIndexes] }
			else {
				seek(fcon,4*(indexBoundsD[1]-1),origin="current")
				currValues <- readBin(fcon,"double",diff(indexBoundsD)+1,size=4,endian="little")
				seek(fcon,-4*indexBoundsD[2],origin="current")
			}
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
				if (compressedSpectrum) { currValues <- spectrum[currIndexes] }
				else {
					seek(fcon,4*(indexBoundsD[1]-1),origin="current")
					currValues <- readBin(fcon,"double",diff(indexBoundsD)+1,size=4,endian="little")
					seek(fcon,-4*indexBoundsD[2],origin="current")
				}
				windowPiece <- rtmsSpectrumPiece(indexToMass(currIndexes-1),currValues)
			}
		} else { windowPiece <- NULL }

		relmax <- which(peakMasses>=currentPeak$bounds[1] & peakMasses<=currentPeak$bounds[2])
		maxima <- rtmsSpectrumPiece(peakMasses[relmax],peakIntensities[relmax])
		subsamples[[i]] <- new_rtmsSubsample(peakPiece,maxima,windowPiece)
	}

	new_rtmsSample(subsamples,peaks)
}

#' Extract a sample set from a Bruker MCF directory
#'
#' Extracts an RTMS sample object (of class `rtmsSampleSet`) from a multi-
#' acquisition Bruker MCF directory opened using an RTMS reader object (of
#' class `rtmsBrukerMCFReader`). A vector  numeric indices is used to identify
#' which acquisitions the sample set should be extracted from.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @param peaks A list of peak objects of class `rtmsPeak`
#' @param indices A vector of numeric indices specifying which acquisitions
#' the sample set should be extracted from
#'
#' @returns An RTMS sample set object of class `rtmsSampleSet`
#'
#' @export
getBrukerMCFSampleSet <- function(reader,peaks,indices) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	if (!all(vapply(peaks,function(p) inherits(p,"rtmsPeak"),FALSE))) {
		stop("Parameter 'peaks' must be a list of objects of class 'rtmsPeak'.")
	}
	if (any(indices<=0) || any(indices>nrow(reader$spotTable))) {
		stop("Parameter 'indices' contains values that are out of bounds.")
	}

	fcon <- reader$con

	subsampleset <- vector("list",length(indices))
	for (index in indices) {
		spotId <- reader$spotTable$id[index]
		fhigh <- reader$spotTable$frequencyHigh[index]
		fwidth <- reader$spotTable$frequencyWidth[index]
		fsize <- reader$spotTable$size[index]
		alpha <- reader$spotTable$alpha[index]
		beta <- reader$spotTable$beta[index]

		massToIndex <- function(p) { return(fsize*(fwidth-(fhigh/(p-alpha/fhigh) + beta))/fwidth) }
		indexToMass <- function(i) { return(fhigh/((fwidth*(fsize-i)/fsize)-beta) + alpha/fhigh) }

		# Move to extracted peak blob.  The sum of peaks within the
		# tranformed index window is the current assay measure
		peakdataBlob <- which(reader$offsetTable$id==spotId & reader$offsetTable$BlobResType==257)
		mcf_blobSeek(fcon,reader$offsetTable$Offset[peakdataBlob],reader$offsetTable$OffsetPage[peakdataBlob])

		mcf_checkBlobCodeWord(fcon)
		bin_checkBytes(fcon,charToRaw("\x01\x01"))
		seek(fcon,16,origin="current")
		numNamedEls <- bin_readVarInt(fcon)
		seek(fcon,25,origin="current") # Skip over "PeakFinderType" element

		name <- mcf_readNamedName(fcon)
		if (name!="Indexes") { stop("Second named element of peak data blob should be 'Indexes") }
		bin_checkBytes(fcon,charToRaw("\x03\x1F")) # Indexes should be an array of 8-byte floats
		numIndexes <- bin_readVarInt(fcon)
		peakIndexes <- readBin(fcon,"double",numIndexes,size=8,endian="little")
		peakMasses <- indexToMass(peakIndexes)

		name <- mcf_readNamedName(fcon)
		if (name!="Intensities") { stop("Third named element of peak data blob should be 'Intensities") }
		bin_checkBytes(fcon,charToRaw("\x03\x20")) # Intensities should be an array of 4-byte floats
		numIntensities <- bin_readVarInt(fcon)
		if (numIndexes!=numIntensities) { stop("Number of peak indicies should match number of peak intensities") }
		peakIntensities <- readBin(fcon,"double",numIndexes,size=4,endian="little")

		# Move to raw data blob.  From here we will extract all measured spectra
		# values within a given width of the desired peaks, and take their average
		rawdataBlob <- which(reader$offsetTable$id==spotId & reader$offsetTable$BlobResType==258)
		mcf_blobSeek(fcon,reader$offsetTable$Offset[rawdataBlob],reader$offsetTable$OffsetPage[rawdataBlob])

		mcf_checkBlobCodeWord(fcon)
		bin_checkBytes(fcon,charToRaw("\x01\x01"))
		seek(fcon,16,origin="current")
		numNamedEls <- bin_readVarInt(fcon)

		name <- mcf_readNamedName(fcon)
		if (name!="Intensities") { stop("First named element of raw data blob should be 'Intensities'") }
		# Intensities should be an array of 4-byte floats
		bin_checkBytes(fcon,as.raw(3),"Raw spectra must be an array.")
		typeByte <- readBin(fcon,"raw",1)
		if (typeByte==as.raw(32)) { # x20 indicates floats are uncompressed
			compressedSpectrum <- FALSE
			numValues <- bin_readVarInt(fcon)
		} else if (typeByte==as.raw(34)) { #x22 indicates floats have been compressed using gzip
			compressedSpectrum <- TRUE
			numBytes <- bin_readVarInt(fcon)
			gzippedBytes <- readBin(fcon,"raw",numBytes)
			unzippedBytes <- memDecompress(gzippedBytes,type="gzip")
			numValues <- length(unzippedBytes)/4
			spectrum <- readBin(unzippedBytes,"double",numValues,size=4,endian="little")
		} else { stop("Raw spectra must be an array of 32-bit floats.") }

		# For each desired peak, extract the relevant peak piece, maxima and, (if
		# applicable) window piece of the overall spectrum. This can be read directly
		# from the file if the data is uncompressed; otherwise it can be taken from
		# the extracted and decompressed spectrum.
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
				if (compressedSpectrum) { currValues <- spectrum[currIndexes] }
				else {
					seek(fcon,4*(indexBoundsD[1]-1),origin="current")
					currValues <- readBin(fcon,"double",diff(indexBoundsD)+1,size=4,endian="little")
					seek(fcon,-4*indexBoundsD[2],origin="current")
				}
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
					if (compressedSpectrum) { currValues <- spectrum[currIndexes] }
					else {
						seek(fcon,4*(indexBoundsD[1]-1),origin="current")
						currValues <- readBin(fcon,"double",diff(indexBoundsD)+1,size=4,endian="little")
						seek(fcon,-4*indexBoundsD[2],origin="current")
					}
					windowPiece <- rtmsSpectrumPiece(indexToMass(currIndexes-1),currValues)
				}
			} else { windowPiece <- NULL }

			relmax <- which(peakMasses>=currentPeak$bounds[1] & peakMasses<=currentPeak$bounds[2])
			maxima <- rtmsSpectrumPiece(peakMasses[relmax],peakIntensities[relmax])
			subsamples[[i]] <- new_rtmsSubsample(peakPiece,maxima,windowPiece)
		}
		subsampleset[[which(indices==index)]] <- subsamples
	}
	if (!is.null(names(indices))) { names(subsampleset) <- names(indices) }

	new_rtmsSampleSet(subsampleset,peaks)
}

#' Retrieve peak intensities directly from an MCF file
#'
#' The size of mass spectrometric data in general, and Bruker MCF directories
#' specifically, makes the extraction of data a resource intensive and time
#' consuming process. `rtms` as a package is designed to reduce this burden,
#' but pulling a sample set from an MCF file can (in the event of compressed
#' spectra) requires reading nearly all data out of a file, which could take an
#' extremely long time over a network connection.  Since peak intensity
#' (calculated as the sum of local intensity maxima within a given peak width)
#' is one of the most common measurements used in evaluating spectra, and
#' because this measure can be extracted without extracting the full spectra,
#' this function aims to avoid expensive reading time by skipping the creation
#' of a sample set object and calculating peak intensity directly.  Other
#' measurements are not possible, but full spectra do not have to be read, even
#' when spectra are compressed, as local maxima are pre-processed and stored
#' separately in a Bruker MCF file.
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @param peaks A list of peak objects of class `rtmsPeak`
#' @param indices A vector of numeric indices specifying which acquisitions
#' the measurements should be taken from
#'
#' @returns A data frame containing columns specifying the index of each
#' acquisition, the name of each acquistion (if `indices` is a named vector),
#' the peak value of the peak measure, the peak name (if `peaks` is a named
#' list), the measure name (which will always be "PeakIntensity"), and the
#' measured value for that sample and peak. Format matches the output of
#' `measureSampleSet`
#'
#' @export
getBrukerMCFIntensities <- function(reader,peaks,indices) {
	if (!inherits(reader,"rtmsBrukerMCFReader")) { stop("Parameter 'reader' must be of class 'rtmsBrukerMCFReader'.") }
	if (!all(vapply(peaks,function(p) inherits(p,"rtmsPeak"),FALSE))) {
		stop("Parameter 'peaks' must be a list of objects of class 'rtmsPeak'.")
	}
	if (any(indices<=0) || any(indices>nrow(reader$spotTable))) {
		stop("Parameter 'indices' contains values that are out of bounds.")
	}

	fcon <- reader$con

	intensityTable <- data.frame()
	# subsampleset <- vector("list",length(indices))
	for (iter in seq_along(indices)) {
		index <- indices[[iter]]
		spotId <- reader$spotTable$id[index]
		fhigh <- reader$spotTable$frequencyHigh[index]
		fwidth <- reader$spotTable$frequencyWidth[index]
		fsize <- reader$spotTable$size[index]
		alpha <- reader$spotTable$alpha[index]
		beta <- reader$spotTable$beta[index]

		massToIndex <- function(p) { return(fsize*(fwidth-(fhigh/(p-alpha/fhigh) + beta))/fwidth) }
		indexToMass <- function(i) { return(fhigh/((fwidth*(fsize-i)/fsize)-beta) + alpha/fhigh) }

		# Move to extracted peak blob.  The sum of peaks within the
		# tranformed index window is the current assay measure
		peakdataBlob <- which(reader$offsetTable$id==spotId & reader$offsetTable$BlobResType==257)
		mcf_blobSeek(fcon,reader$offsetTable$Offset[peakdataBlob],reader$offsetTable$OffsetPage[peakdataBlob])

		mcf_checkBlobCodeWord(fcon)
		bin_checkBytes(fcon,charToRaw("\x01\x01"))
		seek(fcon,16,origin="current")
		numNamedEls <- bin_readVarInt(fcon)
		seek(fcon,25,origin="current") # Skip over "PeakFinderType" element

		name <- mcf_readNamedName(fcon)
		if (name!="Indexes") { stop("Second named element of peak data blob should be 'Indexes") }
		bin_checkBytes(fcon,charToRaw("\x03\x1F")) # Indexes should be an array of 8-byte floats
		numIndexes <- bin_readVarInt(fcon)
		peakIndexes <- readBin(fcon,"double",numIndexes,size=8,endian="little")
		peakMasses <- indexToMass(peakIndexes)

		name <- mcf_readNamedName(fcon)
		if (name!="Intensities") { stop("Third named element of peak data blob should be 'Intensities") }
		bin_checkBytes(fcon,charToRaw("\x03\x20")) # Intensities should be an array of 4-byte floats
		numIntensities <- bin_readVarInt(fcon)
		if (numIndexes!=numIntensities) { stop("Number of peak indicies should match number of peak intensities") }
		peakIntensities <- readBin(fcon,"double",numIndexes,size=4,endian="little")

		# For each desired peak, extract the relevant maxima of the overall
		# spectrum.
		thisTable <- data.frame(index=index,
								peakValue=vapply(peaks,function(p) p$value,0),
								measure="PeakIntensity",
								value=0)
		if (!is.null(names(indices))) { thisTable$sample <- names(indices)[[iter]] }
		else { thisTable$sample <- as.character(index) }
		if (!is.null(names(peaks))) { thisTable$peakName <- names(peaks) }
		for (i in seq_along(peaks)) {
			currentPeak <- peaks[[i]]
			relmax <- which(peakMasses>=currentPeak$bounds[1] & peakMasses<=currentPeak$bounds[2])
			thisTable$value[[i]] <- sum(peakIntensities[relmax])
		}
		intensityTable <- rbind(intensityTable,thisTable)

	}

	intensityTable
}

#' Retrieve specific metadata values from a Bruker BAF file
#'
#' Retrieves a list of specific metadata values (including instrument data,
#' acquisition parameters, processing and analysis directives, etc.) for a
#' specific acquisition from from a Bruker multi-acquisition BAF directory
#' (represented by an `rtmsBrukerBAFReader` object).
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @param names A character vector of metadata names
#' @param index A single numeric index specifying which acquisition the sample
#' set should be extracted from
#'
#' @returns A named list of values corresponding to the metadata values
#' specified.  All values will be returned as a string, including numeric
#' quantities (with units if appropriate).
#'
#' @export
getBrukerMCFMetadata <- function(reader,names,index) {
	meta <- getBrukerMCFAllMetadata(reader,index)
	filtered <- merge(data.frame(DisplayName=names),meta,sort=FALSE)
	values <- as.list(filtered$Value)
	names(values) <- filtered$DisplayName
	values
}

#' Retrieve all metadata values from a Bruker MCF file
#'
#' Retrieves a table of all metadata values (including instrument data,
#' acquisition parameters, processing and analysis directives, etc.) for a
#' specific acquisition from a Bruker multi-acquisition MCF directory
#' (represented by an `rtmsBrukerMCFReader` object).
#'
#' @param reader An RTMS reader object of class `rtmsBrukerMCFReader`
#'
#' @param index A single numeric index specifying which acquisition the sample
#' set should be extracted from
#'
#' @returns A data frame with all metadata parameters for the acquisition. The
#' data frame will have five columns: `Index`, the numeric index of the
#' acquisition; `PermanentName`, the internal identifier of the parameter in
#' Bruker software systems; `GroupName`, the group of parameters that each
#' value belongs to; `DisplayName`, the string used to specify the parameter to
#' users (i.e. how the parameter would be labelled in a user interface); and
#' `Value`, a character column containing the value of each metadata parameter.
#' Numeric quantities will also be returned as strings, with units if
#' appropriate.
#'
#' @export
getBrukerMCFAllMetadata <- function(reader,index) {
	blobRow <- which(reader$offsetTable$id==reader$spotTable$id[index] &
					 	reader$offsetTable$BlobResType==259)
	dectable <- reader$metadata$declarations
	reptable <- reader$metadata$replacements

	fcon <- reader$con
	metaBlobOffset <- reader$offsetTable$Offset[blobRow]
	metaBlobOffsetPage <- reader$offsetTable$OffsetPage[blobRow]
	mcf_blobSeek(fcon,metaBlobOffset,metaBlobOffsetPage)

	mcf_checkBlobCodeWord(fcon)
	bin_checkBytes(fcon,charToRaw("\x01\x01"))
	seek(fcon,16,origin="current")
	numNamedEls <- bin_readVarInt(fcon)

	thisMetaTable <- data.frame()
	for (neiter in seq_len(numNamedEls)) {
		name <- mcf_readNamedName(fcon)
		if (readBin(fcon,"raw",1)==as.raw(0)) { next }
		else { seek(fcon,-1,origin="current") }

		keyTable <- mcf_readKeyValueTable(fcon)

		if (name=="IntValues") {
			names(keyTable)[names(keyTable)=="Value"] <- "Code"
			keyTable <- merge(keyTable,reptable,all.x=TRUE)
			keyTable$Value[is.na(keyTable$Value)] <- as.character(keyTable$Code[is.na(keyTable$Value)])
			keyTable$Code <- NULL
			keyTable <- merge(dectable,keyTable,all.y=TRUE)
		} else if (name=="StringValues") {
			keyTable <- merge(dectable,keyTable,all.y=TRUE)
		} else if (name=="DoubleValues") {
			keyTable <- merge(dectable,keyTable,all.y=TRUE)
			keyTable$Value <- sprintf(keyTable$DisplayFormat,keyTable$Value)
		} else {
			stop("Other metadata values types not supported.")
		}
		thisMetaTable <- rbind(thisMetaTable,keyTable)
	}
	thisMetaTable$Index <- index
	rel <- thisMetaTable$Unit != ""
	thisMetaTable$Value[rel] <- paste(thisMetaTable$Value[rel],thisMetaTable$Unit[rel])
	thisMetaTable[,c("Index","PermanentName","DisplayName","GroupName","Value")]
}

new_rtmsBrukerMCFReader <- function(dir,files,spotTable,offsetTable,metadata,con=NULL) {
	stopifnot(is.character(dir))
	stopifnot(is.list(files))
	stopifnot(is.data.frame(spotTable))
	stopifnot(is.data.frame(offsetTable))
	stopifnot(is.list(metadata))
	if (!is.null(con)) { stopifnot(inherits(con,"file")) }

	structure(list(dir=dir, files=files, spotTable=spotTable,
				   offsetTable=offsetTable, metadata=metadata, con=con),
			  class="rtmsBrukerMCFReader")
}

readBrukerMCFIndexFile <- function(path,calibration=FALSE) {
	# Extract a table of blob information from an mcf_idx file in the SQLite format

	# Open index file for decoding
	scon <- file(path,"rb")
	on.exit(close(scon),add=TRUE)

	# Locate all B-tree table leaf pages for the main table (which lists the
	# titles, column names, and root pages of all tables and index in the
	# SQLite file.
	mainHandler <- function(valueList,cellIndex) {
		if (valueList[[2]]=="ContainerIndex") {
			if (wholenumber(valueList[[4]])) {
				return(data.frame(name="blobRootPage",value=valueList[[4]]))
			} else { stop("Page record column should be an integer") }
		} else if (valueList[[2]]=="Relations") {
			if (wholenumber(valueList[[4]])) {
				return(data.frame(name="relRootPage",value=valueList[[4]]))
			} else { stop("Page record column should be an integer") }
		} else { return(c()) }
	}
	mainTable <- sqlt_parseBTreeTable(scon,1,mainHandler,TRUE)$output

	containerRoot <- mainTable$value[which(mainTable$name=="blobRootPage")]
	if (length(containerRoot)==0) { stop("ContainerIndex table not found.") }
	containerHandler <- function(values,cellIndex) {
		valueList <- values
		names(valueList) <- c("GuidA","GuidB","BlobResType","Offset","BlobSize")
		valueList[["Index"]] <- cellIndex
		# Guid values are 8-byte integers in the file, and hence return two integers
		# each, which must be pasted together to be a single column in a data table
		valueList$id <- paste(c(valueList$GuidA,valueList$GuidB),collapse=" ")
		valueList$GuidA <- NULL
		valueList$GuidB <- NULL
		valueList[["OffsetPage"]] <- 0
		if (length(valueList$Offset)>1) { valueList$OffsetPage <- valueList$Offset[2]*2 }
		if (valueList$Offset[1]<0) {
			valueList$OffsetPage <- valueList$OffsetPage+1
			valueList$Offset <- valueList$Offset[1]+(2^31)
		} else { valueList$Offset <- valueList$Offset[1] }
		return(as.data.frame(valueList))
	}
	containerdf <- sqlt_parseBTreeTable(scon,containerRoot,containerHandler,FALSE)$output

	if (calibration) {
		relationRoot <- mainTable$value[which(mainTable$name=="relRootPage")]
		relationHandler <- function(values,cellIndex) {
			valueList <- values
			names(valueList) <- c("GuidA","GuidB","ToGuidA","ToGuidB","RelationType")
			# Guid values are 8-byte integers in the file, and hence return two integers
			# each, which must be pasted together to be a single column in a data table
			valueList$id <- paste(c(valueList$GuidA,valueList$GuidB),collapse=" ")
			valueList$GuidA <- NULL
			valueList$GuidB <- NULL
			valueList$toId <- paste(c(valueList$ToGuidA,valueList$ToGuidB),collapse=" ")
			valueList$ToGuidA <- NULL
			valueList$ToGuidB <- NULL
			return(as.data.frame(valueList))
		}
		relationdf <- sqlt_parseBTreeTable(scon,relationRoot,relationHandler,FALSE)$output

		containerdf <- merge(containerdf,relationdf,all.x=TRUE)
		containerdf <- containerdf[order(containerdf$Index),c("Index","id","toId","RelationType","BlobResType","Offset","OffsetPage","BlobSize")]
	} else {
		containerdf <- containerdf[order(containerdf$Index),c("Index","id","BlobResType","Offset","OffsetPage","BlobSize")]
	}
	return(containerdf)
}

readMCFCalibration <- function(path,offsetTable) {
	# Open the full Bruker calibration data file for reading
	fsize <- file.size(path)
	fcon <- file(path,"rb")
	rawfile <- readBin(fcon,"raw",fsize+1000)
	close(fcon)

	fcon <- rawConnection(rawfile,"rb")
	on.exit(close(fcon),add=TRUE)

	caltab <- data.frame()
	for (iter in seq_len(nrow(offsetTable))) {
		seek(fcon,offsetTable$Offset[iter])

		# Process blob header
		mcf_checkBlobCodeWord(fcon)
		bin_checkBytes(fcon,charToRaw("\x01\x01"))
		seek(fcon,16,origin="current")
		numNamedEls <- bin_readVarInt(fcon)

		# First named element should be "Transformator" a Bruker list that contains the list of
		# calibration parameters
		name <- mcf_readNamedName(fcon)
		if (name!="Transformator") { stop("First named element of parameter blob should be 'Transformator'") }
		bin_checkBytes(fcon,charToRaw("\x01")) # Transformator should be a named list of parameters
		seek(fcon,16,origin="current")

		# Convert calibration mode 4 to calibration mode 5 for simplicity
		calrow <- mcf_readNamedTableRow(fcon)
		if (calrow$calibMode==4) {
			calrow$alpha <- 0
			calrow$beta <- -calrow$beta
		}
		calrow$toId <- offsetTable$toId[iter]
		caltab <- rbind(caltab,as.data.frame(calrow))
	}
	caltab
}

retrieveMCFMetadata <- function(fcon,offset,offsetPage=0) {
	mcf_blobSeek(fcon,offset,offsetPage)

	# Process blob header
	mcf_checkBlobCodeWord(fcon)
	bin_checkBytes(fcon,charToRaw("\x01\x01"))
	seek(fcon,16,origin="current")
	numNamedEls <- bin_readVarInt(fcon)

	# First named element should be "Declarations" a Bruker data table that includes
	# the integer keys of all metadata elements
	name <- mcf_readNamedName(fcon)
	if (name!="Declarations") { stop("First named element of parameter blob should be 'Declarations'") }
	outTable <- mcf_readTable(fcon)

	# First named element should be "Replacements" a Bruker link table that contains
	# lookup tables for enumerated parameter values
	name <- mcf_readNamedName(fcon)
	if (name!="Replacements") { stop("Second named element of parameter blob should be 'Replacements'") }
	bin_checkBytes(fcon,charToRaw("\x02\x01")) # Declarations should be a index of parameter replacements
	seek(fcon,16,origin="current")
	numReps <- bin_readVarInt(fcon)

	valTable <- data.frame()
	for (nriter in seq_len(numReps)) {
		bin_checkBytes(fcon,as.raw(1))
		seek(fcon,16,origin="current")
		bin_checkBytes(fcon,as.raw(3))
		if (nriter==1) {
			nextName <- mcf_readNamedName(fcon)
		} else { bin_checkBytes(fcon,as.raw(11)) }
		codes <- mcf_readPrimitiveArray(fcon)
		bin_checkBytes(fcon,as.raw(10))
		values <- mcf_readPrimitiveArray(fcon)
		bin_checkBytes(fcon,as.raw(9))
		bin_checkBytes(fcon,as.raw(37))
		key <-  readBin(fcon,"integer",1,size=4,endian="little")
		if (length(codes)!=length(values)) {
			stop("Replacement codes must be matched by the same number of values")
		}
		if (length(codes>0)) {
			valTable <- rbind(valTable,data.frame(Key=key,Code=codes,Value=as.character(values)))
		}
	}

	return(list(declarations=outTable,replacements=valTable))
}

# Helper function to convert a floating-point range of indexes to an integer range
shrinkBounds <- function(bounds) { return(c(ceiling(bounds[1]),floor(bounds[2]))) }

# Checks if a data file represents a whole number
wholenumber <- function(dat) { return(is.numeric(dat) && round(dat)==dat) }

bin_readVarInt <- function(fcon,endian="little") {
	# A varint represents an integer with anywhere from one to nine bytes,
	# with each byte representing a "digit" of the integer in base-128, and each byte
	# but the last with its highest bit set to 1.  For example, the number 13 would just
	# be represented by a single byte "0D", but the number 290 (being equal to 2*128 + 34)
	# would be represented by "22 82" ( 34 in hex, 2 + 128 in hex ) in little endian form
	# and "02 A2" ( 2 in hex, 34 + 128 in hex) in big endian form
	nums <- c()
	for (iter in 1:8) {
		nextByte <- readBin(fcon,"integer",1,size=1,signed=FALSE)
		if (nextByte>127) { nextNum <- nextByte-128 }
		else { nextNum <- nextByte }
		if (endian=="little") { nums <- c(nums,nextNum) }
		else { nums <- c(nextNum,nums) }
		if (nextByte<=127) { break }
	}
	return(sum(nums*(128^(0:(length(nums)-1)))))
}
bin_checkBytes <- function(fcon,bytes,message="Invalid bytes") {
	# Check that the next string of bytes in a file matches a given set
	if (!all(readBin(fcon,"raw",length(bytes))==bytes)) { stop(message) }
}

sqlt_getPage <- function(scon,page) {
	# SQLite files are broken into pages, which we assume here have length
	# 1024; this function jumps to the location of a given page, loads the page
	# into memory and returns a raw connection to that data; it optionally jumps
	# to a given offset within that page
	seek(scon,(page-1)*1024)
	rawpage <- readBin(scon,"raw",1024)
	rawcon <- rawConnection(rawpage,"rb")
	return(rawcon)
}
sqlt_parseBTreeInteriorPage <- function(pagecon,upper) {
	seek(pagecon,2,origin="current")
	numCells <- readBin(pagecon,"integer",1,size=2,signed=FALSE,endian="big")
	cellStart <- readBin(pagecon,"integer",1,size=2,signed=FALSE,endian="big")
	seek(pagecon,1,origin="current")
	finalPage <- readBin(pagecon,"integer",size=4,endian="big")
	cellPointers <- readBin(pagecon,"integer",numCells,size=2,signed=FALSE,endian="big")

	newpagedf <- data.frame(Page=rep(NA,numCells),Upper=rep(NA,numCells),Type=rep("None",numCells))
	for (citer in seq_len(numCells)) {
		seek(pagecon,cellPointers[citer])
		newpagedf$Page[[citer]] <- readBin(pagecon,"integer",size=4,endian="big")
		newpagedf$Upper[[citer]] <- bin_readVarInt(pagecon,"big")
	}
	return(rbind(newpagedf,data.frame(Page=finalPage,Upper=upper,Type="None")))
}
sqlt_parseBTreeLeafPage <- function(pagecon,handler) {
	if (is.null(handler)) { return(c()) }
	seek(pagecon,2,origin="current")
	numCells <- readBin(pagecon,"integer",1,size=2,signed=FALSE,endian="big")
	cellStart <- readBin(pagecon,"integer",1,size=2,signed=FALSE,endian="big")
	seek(pagecon,1,origin="current")
	cellPointers <- readBin(pagecon,"integer",numCells,size=2,signed=FALSE,endian="big")

	rows <- rep(list(c()),numCells)
	for (citer in seq_len(numCells)) {
		seek(pagecon,cellPointers[citer])
		cellSize <- bin_readVarInt(pagecon,"big")
		cellIndex <- bin_readVarInt(pagecon,"big")
		valueList <- sqlt_parseRecord(pagecon)
		rows[[citer]] <- handler(valueList,cellIndex)
	}
	return(as.data.frame(do.call(rbind,rows)))
}
sqlt_parseBTreeTable <- function(scon,rootpage,handler,first=FALSE) {

	pastTheFirst <- !first
	pagedf <- data.frame(Page=rootpage,Upper=Inf,Type="None")
	outputdf <- data.frame()
	while (TRUE) {
		if (!any(pagedf$Type=="None")) { break }
		curpageind <- which(pagedf$Type=="None")[1]
		curpage <- pagedf$Page[curpageind]
		pagecon <- sqlt_getPage(scon,curpage)
		if (!pastTheFirst) {
			seek(pagecon,100)
			pastTheFirst <- TRUE
		}
		pageType <- readBin(pagecon,"integer",1,size=1,signed=FALSE)
		if (pageType==5) {
			pagedf <- rbind(pagedf,sqlt_parseBTreeInteriorPage(pagecon,pagedf$Upper[curpageind]))
			pagedf <- pagedf[-curpageind,]
		} else if (pageType==13) {
			outputdf <- rbind(outputdf,sqlt_parseBTreeLeafPage(pagecon,handler))
			pagedf$Type[curpageind] <- "Leaf"
		} else if (pageType!=5) {
			close(pagecon)
			stop("Page must be b-tree leaf or interior page")
		}
		close(pagecon)
	}
	pagedf <- pagedf[order(pagedf$Upper),]
	return(list(pages=pagedf,output=outputdf))
}
sqlt_parseRecord <- function(scon) {
	# In an SQLite data file, a row of a data table is converted to a string of bytes.
	# The byte string begins with a header of varints; one varint specifying the total
	# size of the header in bytes, and the remaining varints encoding the data types of
	# the elements of the record.  Then each element of the record is included in order
	output <- list()
	headerSize <- bin_readVarInt(scon,"big")
	headerRemaining <- headerSize-1
	dataTypes <- c()
	while (headerRemaining>0) {
		nextCode <- bin_readVarInt(scon,"big")
		dataTypes <- c(dataTypes,nextCode)
		headerRemaining <- headerRemaining-max(ceiling(log(nextCode,base=128)),1)
	}

	for (dtiter in 1:length(dataTypes)) {
		output[[dtiter]] <- sqlt_parseRecordValue(scon,dataTypes[dtiter])
	}
	return(output)
}
sqlt_parseRecordValue <- function(scon,type) {
	# Each possible data type in an SQLite record stream is represented by an integer.
	# The most important data types are represented by
	#   - 1, 2, 3, 4, 5, and 6, representing one, two, three, four, six, and eight
	#     byte signed integers, respectively
	#   - 7, representing an 8-byte (64 bit) floating point number
	#   - Even numbers greater than 12, representing a string of raw binary data
	#   - Odd numbers greater than 13, representing a variable-length string
	# Because R cannot reliably represent 8-byte integer values, data types 5 and 6
	# (representing six and eight-byte integers) actually return two four-byte integers
	# which can be handled later
	if (type==0) { return(NA) }
	else if (type<0 || !wholenumber(type)) { stop("Unsupported record data type") }
	else if (type==8) { return(0) }
	else if (type==9) { return(1) }
	else if (type==12) { return(as.raw(c())) }
	else if (type==13) { return("") }
	else if (type==1) { return(readBin(scon,"integer",1,size=1,signed=TRUE)) }
	else if (type==2) { return(readBin(scon,"integer",1,size=2,signed=TRUE,endian="big")) }
	else if (type==3) {
		tempraw <- c(as.raw(0),readBin(scon,"raw",3))
		return(readBin(tempraw,"integer",1,size=4,endian="big"))
	} else if (type==4) { return(readBin(scon,"integer",1,size=4,endian="big")) }
	else if (type==5) {
		tempraw <- c(as.raw(c(0,0)),readBin(scon,"raw",6))
		return(rev(readBin(tempraw,"integer",2,size=4,endian="big")))
	} else if (type==6) { return(rev(readBin(scon,"integer",2,size=4,endian="big"))) }
	else if (type==7) { return(readBin(scon,"double",1,size=8,endian="big")) }
	else if (type==10 || type==11) { stop("Unsupported record data type") }
	else if ((type%%2)==0) { return(readBin(scon,"raw",(type-12)/2)) }
	else { return(rawToChar(readBin(scon,"raw",(type-13)/2))) }
}

mcf_blobSeek <- function(scon,offset,page) {
	# A hacked little survival function.  It turns out R does not consistently
	# work with integers larger than 2147483647 (2^31-1); most of the offsets in
	# the Bruker raw data file go into the billions, so passing those offsets to
	# seek is not possible. This function represents an offset in terms of "pages"
	# of size 2^31, and then smaller offsets within those pages, which allows us
	# to traverse a 10 GB file
	seek(scon,0)
	if (page>0) {
		for (piter in 1:page) { seek(scon,2^31,origin="current") }
	}
	seek(scon,offset,origin="current")
}
mcf_checkBlobCodeWord <- function(fcon) {
	# Check the beginning of a Bruker raw data file blob for the code word
	# "C0DEAFFE" or "code monkey" in German.
	bin_checkBytes(fcon,charToRaw("\xC0\xDE\xAF\xFE"),"Invalid code word.")
}
mcf_readNamedName <- function(fcon) {
	# In a Bruker raw data file, a named element is marked by four bytes, counting down
	# from "FF FF FF FF", so the first byte must be greater than 127, and the next three
	# bytes must be "FF FF FF" (unless a blob contains more than 127 elements, which we
	# assume is impossible).  The next byte must be "0F", marking the beginning of the
	# name string, followed by a byte representing the length of the name string.  The
	# string itself can be read in ASCII from the next set of bytes corresponding to the
	# string length.
	if (!(readBin(fcon,"integer",1,size=1,signed=FALSE)>127)) { stop("Invalid object name.") }
	bin_checkBytes(fcon,charToRaw("\xFF\xFF\xFF\x0F"),"Invalid object name.")
	namelen <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	return(rawToChar(readBin(fcon,"raw",namelen)))
}
mcf_readPrimitive <- function(fcon) {
	# Extracts primitive data types from the stream of a Bruker raw data file.
	# Primitive data types here mean integers, floats and strings. Each primitive type is
	# marked by a single byte:
	#   - x27 and x28 represent eight-byte integers (probably one signed and the other unsigned)
	#   - x25 and x26 represent four-byte integers (again probably one signed and the other unsigned)
	#   - x20 represents a four-byte (32 bit) floating point number
	#   - x1F represents an eight-byte (64 bit) floating point number
	#   - x2A represents an ASCII string, with a following varint specifying the number of characters
	# Other data type byte markers not handled by this function is x03, which represents an array of
	# items, which is then followed by bytes specifying that type and a varint specifying the number
	# of elements, and x01, which appears to allow specifying a custom data type, often a named list
	dataType <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	if (dataType %in% c(39,40)) { return(readBin(fcon,"integer",1,size=8,endian="little")) }
	else if (dataType %in% c(37,38)) { return(readBin(fcon,"integer",1,size=4,endian="little")) }
	else if (dataType==32) { return(readBin(fcon,"double",1,size=4,endian="little")) }
	else if (dataType==31) { return(readBin(fcon,"double",1,size=8,endian="little")) }
	else if (dataType==0) { return(NULL) }
	else if (dataType==42) {
		stringLength <- bin_readVarInt(fcon)
		return(rawToChar(readBin(fcon,"raw",stringLength)))
	} else { stop(sprintf("Invalid primitive data type %d.",dataType)) }
}
mcf_readPrimitiveArray <- function(fcon) {
	# Reads an array of primitive data values using the same type marking as
	# the previous function.  Arrays are marked by an array byte (x03), a type
	# byte, and a varint listing the number of values to be extracted
	bin_checkBytes(fcon,as.raw(3)) # Arrays are marked with x03
	dataType <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	numEls <- bin_readVarInt(fcon)
	if (numEls==0) { return(c()) }
	outVec <- rep(NA,numEls)
	for (viter in seq_along(outVec)) {
		if (dataType %in% c(39,40)) { outVec[viter] <- readBin(fcon,"integer",1,size=8,endian="little") }
		else if (dataType %in% c(37,38)) { outVec[viter] <- readBin(fcon,"integer",1,size=4,endian="little") }
		else if (dataType==32) { outVec[viter] <- readBin(fcon,"double",1,size=4,endian="little") }
		else if (dataType==31) { outVec[viter] <- readBin(fcon,"double",1,size=8,endian="little") }
		else if (dataType==0) {  }
		else if (dataType==42) {
			stringLength <- bin_readVarInt(fcon)
			outVec[viter] <- rawToChar(readBin(fcon,"raw",stringLength))
		} else { stop(sprintf("Invalid primitive data type %d.",dataType)) }
	}
	return(outVec)
}
mcf_readNamedTableRow <- function(fcon) {
	# In a data table in a Bruker raw data file, there is no separate row specifying column names
	# Instead, the first row of the table *may* contained named elements rather than bare primitives
	# so the names along with their values should be extracted as a list.
	#
	# The format of a table row is as follows:
	#   - A byte specifying the number (N) of columns in the row (its possible this is a varint and
	#     I've simply never seen a table with more than 127 columns, but we'll ignore that for now)
	#   - N elements:
	#      - Either:
	#        - a byte specifying the index of the element in the row (offset to begin at 2 for some
	#          reason, or
	#        - the name of the element
	#      - The primitive element itself
	#
	# The structure is not dissimilar to R lists, in which elements can be specified either by
	# their position, or by a name
	numRowCells <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	output <- list()
	for (rciter in 1:numRowCells) {
		nextByte <- readBin(fcon,"integer",1,size=1,signed=FALSE)
		if (nextByte<=127) {
			if (nextByte!=(rciter+1)) { stop("Invalid row cell index") }
			name <- as.character(nextByte)
		} else {
			seek(fcon,-1,origin="current")
			name <- mcf_readNamedName(fcon)
		}
		prim <- mcf_readPrimitive(fcon)
		if (nextByte<=127) { output[[rciter]] <- prim }
		else { output[[name]] <- prim }
	}
	return(output)
}
mcf_readTableRow <- function(fcon,namevec=NULL) {
	# All rows but the first row (see previous function) in a data table contain unnamed
	# primitive elements.  Ideally the name of each element will have been extracted from
	# the first row, and can be passed to this function to apply the same names to the
	# elements of this row
	numRowCells <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	if (is.null(namevec)) { namevec <- rep("",numRowCells) }
	else if (length(namevec)!=numRowCells) {
		stop("Name vector must contain the same number of elements as the table row")
	}
	output <- list()
	for (rciter in 1:numRowCells) {
		nextByte <- readBin(fcon,"integer",1,size=1,signed=FALSE)
		if (nextByte!=(rciter+1)) { stop("Invalid row cell index") }
		prim <- mcf_readPrimitive(fcon)
		if (namevec[rciter]=="") { output[[rciter]] <- prim }
		else { output[[namevec[rciter]]] <- prim }
	}
	return(output)
}
mcf_readTable <- function(fcon) {
	# Must be marked as an array of similar objects
	bin_checkBytes(fcon,charToRaw("\x03\x01"))
	seek(fcon,16,origin="current")
	numRows <- bin_readVarInt(fcon)

	outTable <- data.frame()
	for (riter in seq_len(numRows)) {
		if (riter==1) {
			curRow <- mcf_readNamedTableRow(fcon)
			namevec <- names(curRow)
			outTable <- rbind(outTable,as.data.frame(curRow))
		} else {
			curRow <- mcf_readTableRow(fcon,namevec)
			outTable <- rbind(outTable,as.data.frame(curRow))

		}
	}
	outTable
}
mcf_readNamedKeyValueRow <- function(fcon) {
	# A key-value table is similar to a generic Bruker raw data file table, with the
	# added constraint that each row must contain two elements, an integer named "Key"
	# and a primitive named "Value"
	bin_checkBytes(fcon,as.raw(2),"Key-value row must contain two elements")
	output <- list()
	for (rciter in 1:2) {
		nextByte <- readBin(fcon,"integer",1,size=1,signed=FALSE)
		if (nextByte<=127) {
			if (nextByte!=(rciter+1)) { stop("Invalid row cell index") }
			if (rciter==1) { name <- "Key" }
			else { name <- "Value" }
		} else {
			seek(fcon,-1,origin="current")
			name <- mcf_readNamedName(fcon)
			if ((rciter==1 && name!="Key") || (rciter==2 && name!="Value")) {
				stop("Names in a key-value row must be 'Key' and 'Value'")
			}
		}
		prim <- mcf_readPrimitive(fcon)
		if (rciter==1) {
			if (!is.integer(prim)) { stop("Key in key-value row must be an integer.") }
			output$Key <- prim
		} else { output$Value <- prim }
	}
	return(output)
}
mcf_readKeyValueRow <- function(fcon) {
	# Extracts an unnamed later row from a Bruker raw data file key-value table
	bin_checkBytes(fcon,as.raw(2),"Key-value row must contain two elements")
	output <- list()
	for (rciter in 1:2) {
		nextByte <- readBin(fcon,"integer",1,size=1,signed=FALSE)
		if (nextByte!=(rciter+1)) { stop("Invalid row cell index") }
		prim <- mcf_readPrimitive(fcon)
		if (rciter==1) {
			if (!is.integer(prim)) { stop("Key in key-value row must be an integer.") }
			output$Key <- prim
		} else { output$Value <- prim }
	}
	return(output)
}
mcf_readKeyValueTable <- function(fcon) {
	# Must be marked as an array of similar objects
	bin_checkBytes(fcon,charToRaw("\x03\x01"))
	seek(fcon,16,origin="current")
	numRows <- bin_readVarInt(fcon)

	keyTable <- data.frame()
	for (riter in 1:numRows) {
		if (riter==1) { curRow <- mcf_readNamedKeyValueRow(fcon) }
		else { curRow <- mcf_readKeyValueRow(fcon) }
		keyTable <- rbind(keyTable,as.data.frame(curRow))
	}
	keyTable
}
mcf_readNumericPrimitiveBlank <- function(fcon) {
	# This function calculates the correct number of bytes to read based on
	# the next data-type byte in a given file stream, useful for skipping
	# over a large stream of similar data elements
	# Obviously, as strings are variable in length, this function cannot be
	# applied to them
	dataType <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	if (dataType%in%c(31,39,40)) { return(8) }
	else if (dataType%in%c(32,37,38)) { return(4) }
	else { stop(sprintf("Invalid primitive data type %d.",dataType)) }
}
mcf_readNamedTableRowBlank <- function(fcon) {
	# Like the previous function, determines the number of bytes that make
	# up each later row in a Bruker raw data file, so that they can be
	# skipped quickly.  This function assumes that all table elements are
	# numeric primitives
	numRowCells <- readBin(fcon,"integer",1,size=1,signed=FALSE)
	count <- 1
	for (rciter in 1:numRowCells) {
		nextByte <- readBin(fcon,"integer",1,size=1,signed=FALSE)
		count <- count+1
		if (nextByte<=127) {
			if (nextByte!=(rciter+1)) { stop("Invalid row cell index") }
		} else {
			seek(fcon,-1,origin="current")
			mcf_readNamedName(fcon)
		}
		primSize <- mcf_readNumericPrimitiveBlank(fcon)
		count <- count+1+primSize
		readBin(fcon,"raw",primSize)
	}
	return(count)
}
