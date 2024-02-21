#Function to check the validity of 'FileExtractor' object class
FileExtractorValidity <- function(object) {
  
  error <- character(0)
  
  nGroups <- ncol(stringr::str_match("", object@regularExpression))-1
  if(nGroups!=2) {
    msg <- paste0("There are ", nGroups, " capturing groups in the regular expression. There must be two.")
    error <- c(error, msg)
  }
  
  if(length(object@capturingGroups)!=2) {
    msg <- "The length of 'capturingGroups' is not equal to 2."
    error <- c(error, msg)
  }
  
  if("plate" %in% object@capturingGroups==F) {
    msg <- "'plate' has not been included in the 'capturingGroups'."
    error <- c(error, msg)
  }
  
  if("sample" %in% object@capturingGroups==F) {
    msg <- "'sample' has not been included in the 'capturingGroups'."
    error <- c(error, msg)
  }
  
  if(length(error)==0) {return(T)} else (return(error))
  
}

#Function to create an object of class 'FileExtractor'
#'@export
#'@description \code{FileExtractor} creates an object which represents the location and arrangement of files in a directory.
#'The default assumes the files for each sample are located in their own folder (with the folder names representing the sample names), and the file names consisting of a plate identified followed by the file extnsion '\code{.jpg.dat}' (i.e. \code{\link[gitter]{gitter}} output files).
#'@param regularExpression Character vector of length 1.
#'This should be a \code{\link[base]{regular expression}} representing relative file paths of the output files of image analysis softwares.
#'This serves two purposes. First, to identify the files which should be analysed.
#'Secondly, the regular expression must indicate which parts of the file paths contain information on sample identity and plate number.
#'This is done using capturing groups - i.e. '(...)'. See \code{\link[stringr]{str_match}}.
#'Therefore, there must be two capturing groups in the regular expression - one to extract the sample, and one to extract the plate.
#'@param capturingGroups Character vector of length 2, containing the strings \code{'sample'} and \code{'plate'}.
#'The order indicates which of the capturing groups indicates the sample, and which of the capturing groups indicates the plate.
#'@return \code{FileExtractor} returns an object of class \code{FileExtractor}.
FileExtractor <- function(regularExpression="(.*)/(.*)\\.jpg\\.dat", capturingGroups=c("sample", "plate")) {
  new("FileExtractor", regularExpression=regularExpression, capturingGroups=capturingGroups)
}

#'@rdname FileExtractor
setClass(
  "FileExtractor",
  slots=c(regularExpression="character", capturingGroups="character"),
  prototype=list(
    regularExpression="(.*)/(.*)\\.jpg\\.dat", capturingGroups=c("sample", "plate")
  ),
  validity=FileExtractorValidity
)

#Function to create a list of files which will go into an object of class 'HDA'
#'@export
#'@rdname FileExtractor
#'@title Specification of directory structure and identification of files
#'@description \code{createFileList} identifies all files in a directory which are matched by the specified \code{fileExtractor}.
#'In addition, numeric elements in the file names are identified and the files are ordered based on this in order to ensure that each plate is correctly mapped the mutants in the library.
#'@param directory Directory within which to search for files. Defaults to the current working directory.
#'@param FileExtractor An object of class \code{FileExtractor} which specifies the strucutre of the sub directories and file names.
#'@param samples2exclude Optional vector specifying any samples to not include in the list.
#'@param plates2exclude Optional vector specifying any plates to not include in the list.
#'@param verbose Logical. Should details of the files identified be printed. Defaults to \code{TRUE}.
#'@return \code{createFileList} returns a \code{\link[S4Vectors]{SimpleList}}. Each element of the list contains a character vector of the files which correspond to a particular sample, with names of each element of the list representing the sample names.
#'@examples
#'#Create a FileExtractor to identify files with the extension '.jpg.txt', where the files for each sample are contained within their own folder.
#'myFileExtractor <- FileExtractor(regularExpression="(.*)/(.*)\\.jpg\\.txt")
#'
#'#Get the path of a directory to be analysed
#'myDirectory <- system.file("extdata", "SGA_data", package="HDAr")
#'
#'#Show all files in the directory
#'list.files(myDirectory, recursive=TRUE)
#'
#'#Create a SimpleList containing all files which are matched by myFileExtractor in the specified directory.
#'myFileList <- createFileList(myDirectory, myFileExtractor)
#'
#'#Create a SimpleList containing all files which are matched by myFileExtractor in the specified directory, but excluding the last plate named 'plates_33-36.jpg.txt'.
#'myFileList <- createFileList(myDirectory, myFileExtractor, plates2exclude="plates_33-36")
createFileList <- function(directory=getwd(), FileExtractor=FileExtractor(), samples2exclude, plates2exclude, verbose=T) {
  
  #Import 'samples2exclude' and 'plates2exclude'
  if(missing(samples2exclude)==F) {
    samples2exclude <- samples2exclude
  } else {
    samples2exclude <- NULL
  }
  if(missing(plates2exclude)==F) {
    plates2exclude <- plates2exclude
  } else {
    plates2exclude <- NULL
  }
  
  #Get all files
  allFiles <- list.files(path=directory, recursive=T)
  
  #Get files which match pattern
  allFiles <- allFiles[grep(FileExtractor@regularExpression, allFiles)]
  
  if(length(allFiles)<1) {
    stop("No files have been found which match the specified pattern")
  }
  
  #Perform matches, get unique sample and plates
  matches <- stringr::str_match(allFiles, FileExtractor@regularExpression)
  samples <-  unique(matches[, which(FileExtractor@capturingGroups=="sample") + 1])
  plates <- unique(matches[, which(FileExtractor@capturingGroups=="plate") + 1])
  
  #Check the samples indicated to exclude are there and exclude if necessary
  if(is.null(samples2exclude)==F) {
    #Check that all 'samples2exclude' are in 'samples'
    if(all(samples2exclude %in% samples)==F) {
      stop(paste0(
        "Could not find all values of 'samples2exclude' in the identified samples:\n",
        paste0(samples, collapse=", "),
        "\n"
      ))
    }
    matches <- matches[-which(matches[, which(FileExtractor@capturingGroups=="sample") + 1] %in% samples2exclude),]
  }
  
  #Check the plates indicated to exclude are there. They will actually be excluded later however
  if(is.null(plates2exclude)==F) {
    #Check that all 'samples2exclude' are in 'samples'
    if(all(plates2exclude %in% plates)==F) {
      stop(paste0(
        "Could not find all values 'plates2exclude' in the identified plates:\n",
        paste0(plates, collapse=", "),
        "\n"
      ))
    }
    matches <- matches[-which(matches[, which(FileExtractor@capturingGroups=="plate") + 1] %in% plates2exclude),]
  }
  
  #Create a list of ordered files for each sample
  fileList <- lapply(unique(matches[, which(FileExtractor@capturingGroups=="sample") + 1]), function(x) {
    
    #Get all plates which match sample
    sampleMatches <- matches[which(matches[, which(FileExtractor@capturingGroups=="sample") + 1] == x),]
    
    #Sort the vector of plates
    plateIds <- sampleMatches[, which(FileExtractor@capturingGroups=="plate") + 1]
    plateIdsSorted <- sort(plateIds)
    
    #Search for the length of the longest common prefex between the first and last plate in the vector
    lcPrefixLength <- Biostrings::lcprefix(plateIdsSorted[1], plateIdsSorted[length(plateIdsSorted)])
    
    #Cut off any common prefix
    plateIds <- substr(plateIds, 1+lcPrefixLength, nchar(plateIds))
    
    #Get the locations and lengths of the first numeric elements in the plate
    numericElements <- stringr::str_match(plateIds, "\\d+")[,1]
    
    if(length(unique(numericElements))!=length(numericElements)) {
      stop(paste0(
        "There are not enough unique numeric elements in the file names for ",
        make.names(x),
        ". The identified elements are ",
        paste(numericElements[order(numericElements)], collapse=", "),
        ". Please edit the file names. Alternatively, check that you have correctly specified the regular expression and the identity of the capturing groups in the 'FileExtractor'"
      ))
    }
    
    #Order the vector of files based on the numeric elements
    fileVector <- sampleMatches[,1][order(as.numeric(numericElements))]
    
    #Pass a message saying how many files have been found
    if(verbose) {
      message(paste0(
        "Found ",
        length(fileVector),
        " files for ",
        make.names(x),
        " and ordered based on numeric elements from the file name:\n",
        paste(numericElements[order(as.numeric(numericElements))], collapse=", "),
        ".\n"
      ))
    }
    
    #Return
    return(fileVector)
    
  })
  
  #Name each element in the list and convert to class 'SimpleList'
  if(is.null(samples2exclude)==F) {
    samples <- samples[-which(samples %in% samples2exclude)]
  }
  names(fileList) <- make.names(samples)
  
  #Check that each vector has equal lengths and return a warning if not
  vectorLengths <- sapply(fileList, length)
  if(length(unique(vectorLengths))!=1) {
    warning("Different numbers of files were found for some samples")
  }
  
  #Return as a simple list
  SimpleList(fileList)
  
}
