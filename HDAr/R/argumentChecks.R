checkHDA <- function(HDA) {
  
  #Check arguments
  if(class(HDA)!="HDA") {
    stop("Please supply an object of class 'HDA'.")
  }
  
}

checkHDAwhat <- function(HDA, what) {
  
  #Check arguments
  if(class(HDA)!="HDA") {
    stop("Please supply an object of class 'HDA'.")
  }
  if(class(what)!="character") {
    stop("Class of 'what' must be character.")
  }
  if(length(what)!=1) {
    stop("'what' must be a length 1 character vector.")
  }
  
  #Find the dataset to work on
  if(what=="auto") {
    what <- activeDataset(HDA)
  }
  if(what %in% names(assays(HDA))==F) {
    stop("'what' must be one of the assays of the HDA object.")
  }
  
  return(what)
  
}

checkExclude <- function(HDA, exclude, remove) {
  
  if(is.null(exclude)==F) {
    
    if(class(exclude)!="character") {
      stop("Class of 'exclude' must be character.")
    }
    
    if(identical(exclude, "all")) {
      exclude <- c(qualityControl(HDA), "footprints")
    }
    
    if(missing(remove)==F) {
      if(any(remove %in% exclude)) {
        exclude <- exclude[-which(exclude %in% remove)]
      }
    }
    
    if(all(exclude %in% c(qualityControl(HDA), "footprints")) == F) {
      stop("'exclude' must represent datasets which are marked as quality control in the HDA object. See '?qualityControl'.")
    }
    
  }
  
  return(exclude)
  
}

checkRowRanges <- function(HDA) {
  
  if(is.null(rowRanges(HDA))) {
    
    stop("'RowRages' have not been added to the 'HDA' object. See '?addMutantPositionalInformation'.")
    
  }
  
}