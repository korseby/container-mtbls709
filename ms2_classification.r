
# ---------- Import spectra from MSP ----------
importMs1Ms2data <- function(fileMs1Path, fileMs2Path, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, mzDeviationAbsolute_grouping, mzDeviationInPPM_grouping){
    # Process ms1 / ms2 data
    source("MetFamily/R_packages.R")
    source("MetFamily/FragmentMatrixFunctions.R")
    source("MetFamily/DataProcessing.R")
    
    # Fixed variables
    proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- 0.9
    mzDeviationAbsolute_mapping <- 0.01
    
    # Box parameters
    parameterSet <- list()
############### TODO: GALAXY PARAMETER
    parameterSet$minimumIntensityOfMaximalMS2peak                  <- minimumIntensityOfMaximalMS2peak
    parameterSet$minimumProportionOfMS2peaks                       <- minimumProportionOfMS2peaks
    parameterSet$mzDeviationAbsolute_grouping                      <- mzDeviationAbsolute_grouping
    parameterSet$mzDeviationInPPM_grouping                         <- mzDeviationInPPM_grouping
    parameterSet$doPrecursorDeisotoping                            <- TRUE
    parameterSet$mzDeviationAbsolute_precursorDeisotoping          <- 0.001
    parameterSet$mzDeviationInPPM_precursorDeisotoping             <- 10
    parameterSet$maximumRtDifference                               <- 0.02
    parameterSet$doMs2PeakGroupDeisotoping                         <- TRUE
    parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping       <- 0.01
    parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping          <- 10
    parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
    parameterSet$mzDeviationAbsolute_mapping                       <- mzDeviationAbsolute_mapping
    parameterSet$minimumNumberOfMS2PeaksPerGroup                   <- 1
    parameterSet$neutralLossesPrecursorToFragments                 <- TRUE
    parameterSet$neutralLossesFragmentsToFragments                 <- FALSE
############### TODO: GALAXY PARAMETER
    
    error <- NULL
    resultObj <- tryCatch(
        {
            convertToProjectFile(
                filePeakMatrix = fileMs1Path, 
                fileSpectra = fileMs2Path, 
                parameterSet = parameterSet, 
                progress = NA
            )
        }, error = function(e) {
            error <<- e
        }
    )
    
    # Errors?
    if(!is.null(error)){
        stop(paste(
            "There occurred an error while processing the input files. Please check the file format and content and try again.", "\n",
            "Occurred error: ", error, sep = ""
        ))
        return()
    }
    if(length(resultObj) == 1){
        if(resultObj == "Number of spectra is zero"){
            stop(paste("There are no MS/MS spectra which fulfill the given criteria. Please refine parameter 'Spectrum intensity' and try 'Import MS\u00B9 and MS/MS data' again."))
            return()
        }
    }
    
    lines <- sparseMatrixToString(matrixRows = resultObj$matrixRows, matrixCols = resultObj$matrixCols, matrixVals = resultObj$matrixVals, parameterSet = parameterSet)
    
    # Process project file
    error <- NULL
    dataList <- tryCatch(
        {
            readProjectData(fileLines = lines, progress = FALSE)
        }, error = function(e) {
            error <<- e
        }
    )
    
    if(!is.null(error)){
        stop(paste(
            "There occurred an error while processing the project file. Please check the file format and content and try again.", "\n",
            "Occurred error: ", error, sep = ""
        ))
        return()
    }
    
    return(dataList)
}



# ---------- Apply classifier on sample spectra ----------
applyClassifierMs2 <- function(classifierFile,
                               propertiesFile,
                               fileMs1Path = NULL,
                               fileMs2Path,
                               fileClasses,
                               minimumIntensityOfMaximalMS2peak = 2000,
                               minimumProportionOfMS2peaks = 0.05,
                               mzDeviationAbsolute_grouping = 0.01,
                               mzDeviationInPPM_grouping = 10) {
    # Classes of interest
    lines <- readLines(con = fileClasses)
    theseLines <- grepl(pattern = "^Organic compounds", x = lines)
    classesCanonical <- lines[which(theseLines)-1]
    classes <- lines[theseLines]
    classesCanonical <- gsub(x = classesCanonical, pattern = ":", replacement = "")
    classes <- gsub(x = classes, pattern = "/", replacement = "; ")

    # Process MS-MS data (msp)
    dataList <- importMs1Ms2data(fileMs1Path, fileMs2Path, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, mzDeviationAbsolute_grouping, mzDeviationInPPM_grouping)
    
    # Load and apply classifier
    source("MetFamily/Annotation.R")
    source("MetFamily/Classifiers.R")
    
    propertiesList <- getClassifierProperties(propertiesFile)
    
    # Load and apply classifier
    error <- NULL
    resultObj <- tryCatch(
        {
            doAnnotation(filePath = classifierFile, propertiesList = propertiesList, featureMatrix = dataList$featureMatrix, parameterSet = dataList$importParameterSet, classesWhiteList = classes, progress = FALSE)
        }, error = function(e) {
            error <<- e
        }
    )
    
    if(!is.null(error)){
        stop(paste(
            "There occurred an error while applying the classifiers. Please check the file format and content and try again.", "\n",
            "Occurred error: ", error, sep = ""
        ))
        return()
    }
    
    classToSpectra_class    <- resultObj$classToSpectra_class
    properties_class        <- resultObj$properties_class
    mappingSpectraToClassDf <- resultObj$mappingSpectraToClassDf
    
    # Box classifier results to data.frame
    annotationRows <- list()
    for(classIdx in seq_along(classToSpectra_class)){
        
        precursorIndeces <- as.integer(names(classToSpectra_class[[classIdx]]))
        pValues          <- unname(classToSpectra_class[[classIdx]])
        
        class <- names(classToSpectra_class)[[classIdx]]
        
        for(idx in seq_along(precursorIndeces)){
            precursorIndex     <- precursorIndeces[[idx]]
            precursorLabel     <- dataList$precursorLabels[[precursorIndex]]
            mz                 <- dataList$dataFrameInfos[[precursorIndex, "m/z"]]
            rt                 <- dataList$dataFrameInfos[[precursorIndex, "RT"]]
            metaboliteName     <- dataList$dataFrameInfos[[precursorIndex, "Metabolite name"]]
            
            presentAnnotations <- unlist(dataList$annoArrayOfLists[[precursorIndex]])
            if(length(presentAnnotations) == 0){
                presentAnnotations <- ""
            } else {
                presentAnnotations <- sort(unlist(presentAnnotations))
                presentAnnotations <- paste(presentAnnotations, collapse = "; ")
            }
            
            pValue <- pValues[[idx]]
            
            annotationRows[[length(annotationRows)+1]] <- c(
                "Index" = precursorIndex, 
                "Label" = precursorLabel, 
                "m/z"   = mz, 
                "RT"    = rt, 
                "Metabolite name" = metaboliteName, 
                "Annotation (present)" = presentAnnotations, 
                "Annotation (putative)" = class, 
                "pValue" = pValue
            )
        }
    }
    
    # Box
    head <- c(
        "Index", 
        "Label", 
        "m/z", 
        "RT", 
        "Metabolite name", 
        "Annotation (present)", 
        "Annotation (putative)", 
        "pValue"
    )
    
    if (is.null(unlist(annotationRows))==TRUE) {
        annotationDf <- data.frame(matrix("", ncol=length(head), nrow=5))
        colnames(annotationDf) <- head
    } else {
        annotationDf <- as.data.frame(t(matrix(data = unlist(annotationRows), nrow = length(head))))
        colnames(annotationDf) <- head
    }
    
    return(annotationDf)
}


