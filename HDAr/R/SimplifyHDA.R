#'@export
SimplifyHDA <- function(HDA, what = "auto", exclude = "all", rowRanges = TRUE) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Extract the desired assay
  assay <- DataFrame(extractAssay(HDA, what, exclude))
  
  #Get the rowData
  rowData <- rowData(HDA)
  
  #Get the rowRanges if desired
  if(rowRanges) {
    if(is.null(rowRanges(HDA))==F) {
      rowRanges <- rowRanges(HDA)
      rowRanges <- DataFrame(seqnames=seqnames(rowRanges), start=start(rowRanges), end=end(rowRanges), strand=strand(rowRanges))
      rowData <- cbind(rowData, rowRanges)
    }
  }
  
  #Combine into one object
  assay <- cbind(rowData, assay)
  
  #Return
  return(assay)
  
}

#'@export
calculateInteractions <- function(HDA, what="auto", exclude="all", limits=c(-2,2), p.value=T, ...) {
  
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  interactions2calculate <- which(!colData(HDA)$control)
  controls <- sapply(interactions2calculate, function(i) {
    which(colData(HDA)$batch == colData(HDA)$batch[i] & colData(HDA)$control)
  })
  
  assay <- as.data.frame(extractAssay(HDA, what, exclude))
  assay <- split(assay, rowData(HDA)$replicate)
  
  results <- lapply(assay, function(i) {
    tests <- do.call("cbind", (lapply(1:length(interactions2calculate), function(j) {
      x <- i[,interactions2calculate[j]]
      y <- i[,controls[j]]
      interaction <- log(median(x, na.rm=T)/median(y, na.rm=T), 10)
      if(p.value) {
        if(all(is.na(x)) | all(is.na(y)) | length(x)==1 | length(y)==1) {
          p.val <- NA
        } else {
          p.val <- t.test(x, y, ...)$p.value
        }
        interaction <- c(interaction, p.val)
      }
      return(interaction)
    })))
  })
  
  interactions <- do.call("rbind", lapply(results, `[`, 1, ))
  interactions <- scales::squish(interactions, limits, only.finite = F)
  colnames(interactions) <- paste0(colnames(HDA)[interactions2calculate], "_Interaction")
  if(p.value) {
    p.values <- do.call("rbind", lapply(results, `[`, 2, ))
    colnames(p.values) <- paste0(colnames(HDA)[interactions2calculate], "_p.value")
    interactions <- cbind(interactions, p.values)
    interactions <- interactions[,rep(1:length(interactions2calculate), each=2) + rep(c(0, length(interactions2calculate)))]
  }
  
  return(DataFrame(interactions))
  
}
