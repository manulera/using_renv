setClassUnion("GRanges_Or_NULL", c("GRanges", "NULL"))

validHDA <- function(object) {
  
  #Create a vector to store error messages
  error <- character(0)
  
  #Check that the library file is OK by initialising a HDALibrary object from the rowData
  libraryFile <- HDALibrary(object)
  
  #Check that the batch is indicated
  if(is.null(colData(object)$batch)) {
    msg <- "Batches must be indicated in the colData, even if there is only one batch."
    error <- c(error, msg)
  }
  
  #Check that the sample is indicated
  if(is.null(colData(object)$sample)) {
    msg <- "Samples must be indicated in the colData and match the column names of the HDA object."
    error <- c(error, msg)
  }
  
  #Make sure that the col names match the sample column of the colData
  if(identical(colnames(object), colData(object)$sample)==F) {
    msg <- "The sample column of the colData must be identical to the column names of the HDA object."
    error <- c(error, msg)
  }
  
  #If there is a control, do checks
  if(is.null(control(object))==F) {
    #Check that the control(s) are in the samples
    if(all(control(object) %in% colnames(object)) == F) {
      msg <- "The control(s) must share the names of the samples."
      error <- c(error, msg)
    }
    #Check that there is a control column in the colData
    if(is.null(colData(object)$control)) {
      msg <- "There must be a 'control' column in the colData."
      error <- c(error, msg)
    } else {
      if(is.logical(colData(object)$control) == F) {
        msg <- "The class of the 'control' column in the colData must be 'logical'."
        error <- c(msg, error)
      } else {
        #Check that the controls in the metadata match the controls in the colData
        metadataControl <- control(object)
        colDataControl <- rownames(colData(object))[colData(object)$control]
        metadataControl <- metadataControl[order(metadataControl)]
        colDataControl <- colDataControl[order(colDataControl)]
        if(identical(metadataControl, colDataControl) == F) {
          msg <- "The control(s) indicated in the metadata must match the control(s) indicated in the colData."
          error <- c(error, msg)
        }
        #Check that each batch has one control
        batchSums <- plyr::ddply(as.data.frame(colData(object)), plyr::.(batch), plyr::summarize, sum=sum(control))
        if(all(batchSums$sum==1)==F) {
          msg <- "There must be only one control per batch."
          error <- c(error, msg)
        }
      }
    }
  }
  
  #Make sure replicate is not one of the columns names of the colData
  if("replicate" %in% colnames(colData(object))) {
    msg <- "'replicate' is a reserved name for the columns names of the colData. "
    error <- c(error, msg)
  }
  
  if(length(assays(object))>0) {
    #Check that there is an active dataset
    if(is.null(activeDataset(object))) {
      msg <- "One of the assays must be marked as the activeDataset."
      error <- c(error, msg)
    } else {
      #Check that the active dataset is one of the assays
      if(activeDataset(object) %in% names(assays(object)) == F) {
        msg <- "The activeDataset must be one of the names of the assays."
        error <- c(error, msg)
      }
      #Check invalid assay names
      invalidAssayNames <- c("all", "footprints", "auto")
      if(any(names(assay(object)) %in% invalidAssayNames)) {
        msg <- paste0("Assay names cannot be '", paste(invalidAssayNames, collapse="', '"), "'.")
        error <- c(error, msg)
      }
    }
  }
  
  #Make sure that the column names are not shared with HDALibrary
  if(any(colnames(object) %in% colnames(HDALibrary()))) {
    msg <- "Column names cannot be the same as columns names for a 'HDALibrary' object."
    error <- c(error, msg)
  }
  
  #Check other invalid column names
  invalidColNames <- c("replicate", "all", "seqnames", "start", "end", "strand")
  if(any(colnames(object) %in% invalidColNames)) {
    msg <- paste0("Column names cannot be '", paste(invalidColNames, collapse="', '"), "'.")
    error <- c(error, msg)
  }
  
  #Make sure that anything listed as a quality control dataset is present in the assays
  if(is.null(qualityControl(object))==F) {
    if(any(qualityControl(object) %in% names(assays(object)) == F)) {
      msg <- "Anything marked as quality control must be one of the assays."
      error <- c(error, msg)
    }
  }
  
  #Make sure that the organism is supported
  if(is.null(organism(object))==F) {
    supportedOrganisms <- c("Saccharomyces cerevisiae", "Schizosaccharomyces pombe")
    if(organism(object) %in% supportedOrganisms == F) {
      msg <- paste0(
        "The organism is not currently supported in 'HDAr'. Currently supported organisms are:\n",
        paste0(supportedOrganisms, collapse = "\n")
      )
      error <- c(error, msg)
    }
  }
  
  #Make sure that any rowRanges are the same length as other rows and that they have the same name
  if(is.null(rowRanges(object))==F) {
    if(nrow(object) != length(rowRanges(object))) {
      msg <- "The length of the rowRanges must be the same as the number of rows of the HDA."
      error <- c(error, msg)
    }
    if(identical(rownames(object), names(rowRanges(object)))==F) {
      msg <- "The row names of the HDA must match the names of the rowRanges of the HDA."
      error <- c(error, msg)
    }
  }
  
  #Make sure that colRanges are the same length as the number of columns
  if(length(colRanges(object))!=ncol(object)) {
    msg <- "The length of the colRanges must be the same as the number of columns of the HDA."
    error <- c(error, msg)
  }
  
  #Make sure that the col names match the names of the colRanges
  if(identical(colnames(object), names(colRanges(object)))==F) {
    msg <- "The column names of the HDA must match the names of the colRanges of the HDA."
    error <- c(error, msg)
  }
  
  #Return any errors if necessary
  if(length(error)==0) {return(T)} else {return(error)}
  
}

#'@include HDALibrary.R
#'@include SampleIdentities.R
setClass(
  "HDA",
  contains = "SummarizedExperiment",
  prototype = list(colData=DataFrame(batch=logical(), sample=character()), elementMetadata=HDALibrary()),
  slots = c(
    rowRanges = "GRanges_Or_NULL",
    colRanges = "GRangesList"
  ),
  validity = validHDA
)

setMethod("HDALibrary", "HDA", function(libraryFile) {
  HDALibrary(rowData(libraryFile))
})

#Function to construct a matrix of raw colony sizes from files from image analysis software
constructRSMatrix <- function(files, directory, imageAnalysisSoftware) {
  oldDir <- getwd()
  on.exit({setwd(oldDir)})
  setwd(directory)
  
  #Create a list, with each element containing a vector of all colony sizes for a particular sample
  rawSizes <- lapply(1:length(files), function(x) {
    
    #For each sample
    #Get the sample name
    name <- names(files)[x]
    
    #Get the files
    x <- files[[x]]
    
    #Create a list, with each element a data.frame for a particular file
    fileList <- lapply(1:length(x), function(y) {
      if(imageAnalysisSoftware=="gitter") {
        tryCatch({
          file <- gitter::gitter.read(x[y])
        }, error=function(e) {
          stop(paste0("Could not read 'gitter' file: '",
                      x[y],
                      "'. Please make sure you have specified the correct file format."
          ))
        })
      }
      if(imageAnalysisSoftware=="spotsizer") {
        tryCatch({
          file <- read.csv(x[y], stringsAsFactors=F)
        }, error=function(e) {
          stop(paste0("Could not read 'spotsizer' file: '",
                      y,
                      "'. Please make sure you have specified the correct file format."
          ))
        })
      }
      return(file)
    })
    
    #Extract the colony size information and save as a vector
    fileList <- lapply(1:length(fileList), function(y) {
      name <- x[y]
      y <- fileList[[y]]
      tryCatch({
        if(imageAnalysisSoftware=="gitter") {
          y <- y$size
          type <- "size"
        }
        if(imageAnalysisSoftware=="spotsizer") {
          y <- y$area
          type <- "area"
        }
      }, error=function(e) {
        stop(paste0(
          "Could not find '",
          type,
          "' column in the '",
          imageAnalysisSoftware,
          "' file: '",
          name,
          "'. Please make sure you have specified the correct file format."
        ))
      })
      if(is.numeric(y)==F) {
        stop(paste0(
          "The '",
          type,
          "' column in the '",
          imageAnalysisSoftware,
          "' file: '",
          name,
          "' does not contain numeric data. Please make sure the file is correct."
        ))
      }
      return(y)
    })
    
    #Get the number of colonies in each file
    dims <- sapply(fileList, length)
    
    #If the number of colonies in each file is not equal, return an error
    if(length(unique(dims))!=1) {
      error <- paste0(
        "The number of rows in each file supplied for '",
        name,
        " 'is not equal.\n",
        paste0(unlist(lapply(1:length(files), function(i) {
          paste0("'", x[[i]], "' has ", dims[i], " rows.\n")
        })), collapse="")
      )
      stop(error)
    }
    
    #Combine colony sizes into a single vector and return
    sizeVector <- do.call("c", fileList)
    return(sizeVector)
  })
  
  #Get the number of colonies in each sample
  dims <- sapply(rawSizes, length)
  
  #If the number of colonies in each sample is not equal, return an error
  if(length(unique(dims))!=1) {
    error <- paste0(
      "The number of colonies for each sample is not equal.\n",
      paste0(unlist(lapply(1:length(rawSizes), function(i) {
        paste0("'", names(files)[i], "' has ", dims[i], " colonies.\n")
      })), collapse="")
    )
    stop(error)
  }
  
  #Create a matrix of colony sizes and return
  rawSizes <- do.call("cbind", rawSizes)
  colnames(rawSizes) <- names(files)
  return(rawSizes)
}

#Create a function to construct an object of class 'HDA' from output files from an image analysis software
#'@export
#'@title High-density Array
#'@description \code{HDA} is an class for the representation of data from high-throughput microbial screens such as synthetic genetic arrays or chemical screens.
HDA <- function(files, directory=getwd(), imageAnalysisSoftware="gitter", libraryFile, batches=1) {
  rawSizes <- constructRSMatrix(files, directory, imageAnalysisSoftware)
  colRanges <- GRangesList(rep(list(GRanges()), ncol(rawSizes)))
  names(colRanges) <- colnames(rawSizes)
  colData <- DataFrame(matrix(rep(1, ncol(rawSizes))))
  colnames(colData) <- "batch"
  colData$batch <- batches
  colData$sample <- colnames(rawSizes)
  new(
    "HDA",
    SummarizedExperiment(
      assays=list(RS=rawSizes),
      rowData=HDALibrary(libraryFile),
      colData=colData,
      metadata=list(activeDataset="RS", pipeline=SimpleList())
    ),
    colRanges=colRanges
  )
}

