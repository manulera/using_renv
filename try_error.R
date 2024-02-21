#A blank library
blankHDALibrary <- DataFrame(
  "library_mutant"=character(0),
  "replicate"=character(0),
  "plate"=integer(0),
  "row"=integer(0),
  "col"=integer(0)
)

#A function to check the validity of an object of class 'HDALibrary'
validHDALibrary <- function(object) {
  
  #Create a vector to store error messages
  error <- character(0)
  
  #Check the column names of the library
  col_names <- colnames(object)
  if(any(colnames(blankHDALibrary) %in% col_names == F)) {
    msg <- "Column names for the library file do not include 'library_mutant', 'replicate', 'plate', 'row', and 'col'."
    error <- c(error, msg)
  }
  
  #If the correct columns are present
  if(length(error)==0) {
    
    #Check the classes of each of the columns in the library
    if(is.character(object$library_mutant)==F) {
      msg <- "The class of the 'library_mutant' column is not 'character'."
      error <- c(error, msg)
    }
    if(is.character(object$replicate)==F) {
      msg <- "The class of the 'replicate' column is not 'character'."
      error <- c(error, msg)
    }
    if(is.integer(object$plate)==F) {
      msg <- "The class of the 'plate' column is not 'integer'."
      error <- c(error, msg)
    }
    if(is.integer(object$row)==F) {
      msg <- "The class of the 'row' column is not 'integer'."
      error <- c(error, msg)
    }
    if(is.integer(object$col)==F) {
      msg <- "The class of the 'col' column is not 'integer'."
      error <- c(error, msg)
    }
    
    #If the classes are correct
    if(length(error)==0) {
      
      #Make sure each position in the library is unique
      positions <- DataFrame("plate"=object$plate, "row"=object$row, "col"=object$col)
      if(nrow(unique(positions))!=nrow(positions)) {
        msg <- "Not every position in the library (combination of 'plate', 'row' and 'col') is unique."
        error <- c(error, msg)
      }
      
    }
    
  }
  
  #Return any errors if necessary
  if(length(error)==0) {return(T)} else {return(error)}
  
}

setClass(
  "HDALibrary",
  contains = "DataFrame",
  slots=c("library_mutant", "replicate", "plate", "row", "col"),
  prototype = blankHDALibrary,
  validity = validHDALibrary
)

#A function to fill in empty positions with NAs
fillEmpty <- function(libraryFile) {
  
  #Create an expanded grid containing every possible combination of plate, row and col
  newLibraryFile <- DataFrame(expand.grid(plate=1:max(libraryFile$plate), row=1:max(libraryFile$row), col=1:max(libraryFile$col)))
  
  #Merge the expaned grid with the original library file, filling missing positions with NA
  newLibraryFile <- base::merge(newLibraryFile, libraryFile, all=T)
  
  #Re-order the new library file by plate, row then col
  newLibraryFile <- newLibraryFile[order(newLibraryFile$plate, newLibraryFile$row, newLibraryFile$col),]
  
  #Rearrange the column order
  newLibraryFile <- HDALibrary(newLibraryFile[, c(4:ncol(newLibraryFile), 1:3)])
  
  #Return
  return(newLibraryFile)
  
}

#'@export
#'@rdname HDALibrary
#'@title Library File Manipulation
#'@description \code{HDALibrary} is a constructor function for objects of class \code{HDALibrary}.
#'This class is essentially a \code{\link[S4Vectors]{DataFrame}}, but must have five columns which describes the mutant library being used as follows:
#'\code{library_mutant} should contain the Systematic IDs of the mutants in the library.
#'\code{replicate} is used to identify and aggregate technical replicates. For this, any rows which have the same value for \code{replicate} will be treated as technical replicates.
#'In the standard case, this column will be the same as the \code{library_mutant} column, meaning should mutants appear multiple times in the library, they will be treated as technical replicates.
#'Should this behaviour not be desired, then each time a mutant appears in the library it should be marked with a different value in the \code{replicate} column.
#'\code{plate}, \code{row} and \code{col} should contain the positional information of each mutant.
#'@param libraryFile Argument to specify which library file to use.
#'In the case that this is class \code{character}, this will be matched against the libraries available in this version of the package, which can be listed using \code{listLibraryFiles}.
#'Should a unique match be found, this library will be retrieved from the internal package data.
#'Should no unique match be found, \code{libraryFile} will be assumed to represent the file path to a \code{csv} file containing a custom library.
#'Alternatively, an object of class \code{\link[base]{data.frame}} or \code{\link[S4Vectors]{DataFrame}} can be supplied, from which a library file will be constucted.
#'@return \code{HDALibrary} returns an object of class \code{HDALibrary}.
#'@seealso \code{\link{libraryFile}}
#'@examples
#'#Create a new library file from scratch of 1 plate in 96 well format
#'myHDALibrary <- expand.grid(plate=1L, row=1:8, col=1:12)
#'myMutants <- paste("mutant", 1:96, sep="_")
#'myHDALibrary <- DataFrame(library_mutant=myMutants, replicate=myMutants, myHDALibrary)
#'myHDALibrary <- HDALibrary(myHDALibrary)
#'
#'#Create a library file by reading from a csv
#'myHDALibrary <- system.file("extdata", "HDALibrary_example.csv", package="HDAr")
#'myHDALibrary <- HDALibrary(myHDALibrary)
#'
#'#List the available library files
#'availableLibraries <- listLibraryFiles()
#'
#'Get the Fission Yeast Bioneer V5 96 Well Format Library
#'myHDALibrary <- HDALibrary("Fission_Yeast-Bioneer_V5-96_Well_Format")
#'
#'#Convert into 384 Well Format by condesing from 36 plates to 9 plates
#'myHDALibrary <- rearrangeLibrary(myHDALibrary, c(16, 24), "condense")
#'
#'#Convert into 1536 (quadruplicate) Well Format by replicating
#'myHDALibrary <- rearrangeLibrary(myHDALibrary, c(32, 48), "replicate")
setGeneric("HDALibrary", function(libraryFile) {
  standardGeneric("HDALibrary")
})

setMethod("HDALibrary", "missing", function(libraryFile) {
  new("HDALibrary")
})

setMethod("HDALibrary", "data.frame", function(libraryFile) {
  libraryFile <- DataFrame(libraryFile)
  callGeneric()
})

setMethod("HDALibrary", "DataFrame", function(libraryFile) {
  new("HDALibrary", libraryFile)
})

setMethod("HDALibrary", "character", function(libraryFile) {
  libraryFile <- tryCatch({
    match.arg(libraryFile, listLibraryFiles())
    attributes(libraryFile)$error <- F
    libraryFile
  }, error=function(e) {
    attributes(libraryFile)$error <- T
    libraryFile
  })
  if(attributes(libraryFile)$error) {
    libraryFile <- read.csv(libraryFile, stringsAsFactors = F)
    callGeneric()
  }
  libraryFileSource <- librarySourceFile[which(listLibraryFiles()==libraryFile)]
  nTransformations <- sum(as.numeric(librarySourceFile[1:which(listLibraryFiles()==libraryFile)]==libraryFileSource))-1
  libraryFile <- libraryFileList[[libraryFileSource]]
  if(nTransformations>=1) {
    for(j in 1:nTransformations) {
      libraryFile <- rearrangeLibrary(libraryFile, c(max(libraryFile$row)*libraryFileTransformations[[libraryFileSource]][[j]][[2]][[1]], max(libraryFile$col)*libraryFileTransformations[[libraryFileSource]][[j]][[2]][[2]]), libraryFileTransformations[[libraryFileSource]][[j]][[1]])
    }
  }
  callGeneric()
})

#A function to take a library file and produce a new one in a different pinning format
#'@export
#'@rdname HDALibrary
#'@description \code{rearrangeLibrary} can be used to rearrange an existing library into a new pinning format.
#'@param outputFormat Length 2 \code{numeric} vector containing the row and columns lengths of the new pinning format.
#'@param method The method to be used, either "\code{condense}" or "\code{replicate}".
#'"\code{condense}" specifies that the library is being condensed into a smaller number of plates.
#'In this case, the factor by which the library is being condensed must be a factor of the original number of plates.
#'"\code{replicate}" specifies that technical replicates are being added by repeatedly pinning the library in adjacent positions.
#'@param byrow In the case that the library is being condensed, should positions in the new library be filled by row (as opposed to by column)?
#'Defaults to \code{TRUE}.
#'@return \code{rearrangeLibrary} returns a new object of class \code{HDALibrary} in the desired pinning format.
rearrangeLibrary <- function(libraryFile, outputFormat, method=c("condense", "replicate"), byrow=T) {
  
  #Check arguments
  if(any(method %in% c("condense", "replicate"))==F) {
    stop("Please make sure that the value of 'method' is either 'condense' or 'replicate'")
  }
  
  #Check that the dimensions of the two libraries are compatible
  if(outputFormat[1] %% max(libraryFile$row) != 0) {
    stop("Please make sure the number of rows in the outputFormat is a multiple of the number of rows in the library")
  }
  if(outputFormat[2] %% max(libraryFile$col) != 0) {
    stop("Please make sure the number of columns in the outputFormat is a multiple of the number of columns in the library")
  }
  
  #Fill in any missing positions
  libraryFile <- fillEmpty(libraryFile)
  
  #Calculate by what factor the nubmber of rows and columns are being increased
  rowExpansion <- outputFormat[1] / max(libraryFile$row)
  colExpansion <- outputFormat[2] / max(libraryFile$col)
  
  #Calculate the total increase in colonies on the plate
  totalExpansion <- rowExpansion * colExpansion
  
  #If the library is being condensed, check that the number of plates in the library is divisible by the factor by which the library is being condensed
  if(method=="condense") {
    if(max(libraryFile$plate) %% totalExpansion !=0) {
      stop("The number of plates in the original library is not a divisible by the factor by which the library is being condensed.")
    }
  }
  
  #Create a matrix to define the positions of the condensed/replicated plates
  x <- matrix(1:totalExpansion, nrow=rowExpansion, byrow=byrow)
  x <- x[rep(seq_len(nrow(x)), max(libraryFile$row)), rep(seq_len(ncol(x)), max(libraryFile$col))]
  
  #Get the number of new plates to create and a list showing how plates from the old library are to be placed into the new plates
  if(method=="condense") {
    nNewPlates <- max(libraryFile$plate) / totalExpansion
    plates2subset <- as.list(as.data.frame(matrix(1:max(libraryFile$plate), ncol=nNewPlates)))
  }
  if(method=="replicate") {
    nNewPlates <- max(libraryFile$plate)
    plates2subset <- lapply(1:max(libraryFile$plate), rep, totalExpansion)
  }
  
  #Create a new empty library file
  newLibraryFile <- HDALibrary()
  
  #Get column names to rearrange
  cols <- colnames(libraryFile)
  cols <- colnames(libraryFile)[-which(colnames(libraryFile) %in% c("plate", "row", "col"))]
  
  #For each new plate to be created
  for(j in 1:nNewPlates) {
    
    #Create a new plate to fill
    newPlate <- rep(list(x), length(cols))
    
    #For each new expanded position
    for(i in 1:totalExpansion) {
      
      #Subset the library file to get the right plate and make sure it is correctly ordered
      libraryFileSub <- libraryFile[libraryFile$plate==plates2subset[[j]][i],]
      libraryFileSub <- libraryFileSub[order(libraryFileSub$row, libraryFileSub$col),]
      
      #Create a matrix of the old plate
      oldPlate <- list()
      oldPlate <- lapply(cols, function(k) {
        matrix(DataFrame(libraryFileSub)[,k], nrow=max(libraryFile$row), byrow=T)
      })
      
      #Fill in the new plate at the appropriate position
      for(k in 1:length(cols)) {
        newPlate[[k]][x==i] <- oldPlate[[k]]
      }
      
    }
    
    #Construct a library file for this plate
    libraryPositions <- expand.grid(row=1:nrow(x), col=1:ncol(x))
    libraryPositions <- libraryPositions[order(libraryPositions$row, libraryPositions$col),]
    libraryPositions$plate <- as.integer(j)
    for(i in 1:length(cols)) {
      libraryPositions[,cols[i]] <- as.vector(t(newPlate[[i]]))
    }
    libraryPositions <- HDALibrary(libraryPositions[, c(4:ncol(libraryPositions), 3, 1, 2)])
    
    #Append the new plate to the library file
    newLibraryFile <- HDALibrary(rbind(newLibraryFile, libraryPositions))
    
  }
  
  return(newLibraryFile)
  
}

#Simple conversion of numeric to character describing how many replicates
nReplicates <- function(n) {
  if(n <= length(nReplicatesVector)) {return(nReplicatesVector[n])}
  paste(n, "replicates")
}

#A function to list the available library files
#'@export
#'@rdname HDALibrary
#'@description \code{listLibraryFiles} is a function which returns the default libraries available in this package version.
#'@return \code{listLibraryFiles} returns a vector containing the default libraries available in this package version.
listLibraryFiles <- function() {
  
  #Get the original library names
  originalNames <- names(libraryFileList)
  
  unlist(lapply(1:length(originalNames), function(i) {
    
    #Create a vector to contain the new library names which will be constructed from this single original
    newNames <- originalNames[i]
    
    #If there are transformation to be conducted
    if(length(libraryFileTransformations[[originalNames[i]]])>=1) {
      
      #For each transformation
      for(j in 1:length(libraryFileTransformations[[originalNames[i]]])) {
        
        #Create a new name (from the old)
        newName <- newNames[j]
        newName <- strsplit(newName, "-")
        
        #Calculate the new library density from the old library name and the expansion factor
        originalDensity <- newName[[1]][3]
        originalDensity <- as.numeric(regmatches(originalDensity, regexec("^\\d+", originalDensity))[[1]])
        newDensity <- originalDensity*prod(libraryFileTransformations[[originalNames[i]]][[j]][[2]])
        
        #Add a comment about replicates if neccesary
        if(libraryFileTransformations[[originalNames[i]]][[j]][[1]]=="replicate") {newDensity <- paste0(newDensity, "_(", nReplicates(newDensity/originalDensity), ")")}
        
        #Construct the new name
        newName[[1]][3] <- gsub("^\\d+", newDensity, newName[[1]][3])
        newName <- paste(newName[[1]], collapse = "-")
        
        #Append the new name to the vector of names
        newNames <- c(newNames, newName)
      }
    }
    
    #Return
    newNames
    
  }))
}