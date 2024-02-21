setMethod("[", "HDA", function(x, i, j, ..., drop = TRUE) {
  x <- callNextMethod()
  #Do not use 'rowRanges<-' or 'colRanges<-' as intermediate object may fail validity test
  if(missing(i)==F) {
    x@rowRanges <- rowRanges(x)[i]
  }
  if(missing(j)==F) {
    x@colRanges <- colRanges(x)[j]
    control(x) <- control(x)[which(control(x) %in% colnames(x))]
  }
  validObject(x)
  return(x)
})

#####Need to edit for new rowRanges
setMethod("rbind", "HDA", function(..., deparse.level = 1) {
  
  #Get args
  args <- unname(list(...))
  
  #Get the row ranges for each HDA object
  rowRangesList <- lapply(args, rowRanges)
  
  #If there are a mixture of classes (i.e. some have row ranges and some do not), put blank row ranges in place of missing row ranges
  if(length(unique(sapply(rowRangesList, class)))!=1) {
    stop("Some HDA objects contain rowRanges and some do not. Please ensure that either all HDA objects contain rowRanges or that no HDA objects contain rowRanges.")
  }
  
  #Remove row ranges from arguments
  args <- lapply(args, function(i) {
    rowRanges(i) <- NULL
    return(i)
  })
  
  #Combine HDA objects using SummarizedExperiment method
  HDA <- do.call("callNextMethod", c(args, deparse.level = deparse.level))
  
  #If all are GRangesLists, then combine and add
  if(all(sapply(rowRangesList, is, "GRangesList"))) {
    rowRanges(HDA) <- do.call("c", rowRangesList)
  }
  
  #Check validity
  validObject(HDA)
  
  #Return
  return(HDA)
  
})

setMethod("cbind", "HDA", function(..., deparse.level = 1) {
  
  #Get args
  args <- unname(list(...))
  
  
  
})
