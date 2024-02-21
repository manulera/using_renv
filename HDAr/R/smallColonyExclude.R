otsuThreshold <- function(colonySizes, log=TRUE) {
  
  #Remove NAs
  if(any(is.na(colonySizes))) {
    colonySizes <- colonySizes[-which(is.na(colonySizes))]
  }

  #Log transform data
  if(log) {
    colonySizes <- log(colonySizes+1)
  }
  
  #Create a vector of breaks (thresholds) to test
  uniqueSizes <- unique(colonySizes)
  uniqueSizes <- uniqueSizes[order(uniqueSizes)]
  breaks <- zoo::rollmean(uniqueSizes, 2)
  
  #Get the total number of data points
  counts <- length(colonySizes)
  
  #Calculate the between class variance for each possible threshold
  betweenClassVarianceVector <- sapply(breaks, function(i) {
    absent <- colonySizes<i
    nAbsent <- sum(absent)
    weightAbsent <- nAbsent/counts
    weightPresent <- 1-weightAbsent
    meanAbsent <- mean(colonySizes[absent])
    meanPresent <- mean(colonySizes[!absent])
    betweenClassVariance <- weightAbsent*weightPresent*((meanAbsent-meanPresent)^2)
  })
  
  #Set the threshold as the break value which maximises the between class variance
  threshold <- breaks[which.max(betweenClassVarianceVector)]
  threshold <- exp(threshold) - 1
  
  return(threshold)
  
}

#'@export
#'@title Exclusion of Small Colonies
#'@description Extremely small colonies in the control dataset should be removed from all datasets, as these often represent strains which did not wake up.
#'No useful interaction information can be obtained from these strains, whilst their inclusion will disrupt subsequent normalisation processes.
#'@details For simplicity, technical replicates are treated together.
#'Thus, if the median size of the technical replicate is below the threshold, all techincal replicates are excluded.
#'For details on how to specify technical replicates, see \code{\link{HDALibrary}}.
#'In addition, there may be cases where there are multiple batches, and each batch has a different control dataset.
#'In these cases, each control will be analysed separately, and then colonies identified as absent will only be excluded from that particular batch.
#'@param HDA An object of class \code{\link{HDA}}.
#'@param method The thresholding method. Currently supported options are "\code{otsu}" or "\code{manual}".
#'For otsu thresholding, the threshold is automatically decided based on minimising the within-class variance of the log-transformed colony sizes.
#'Defaults to "\code{otsu}".
#'@param manualThreshold In the case that \code{method} is set to "\code{manual}", then the manually chosen threshold should be supplied here.
#'@param what Which dataset of the \code{\link{HDA}} object should be thresholded? Defaults to "\code{RS}" (raw size).
#'@return \code{smallColonyExclude} will return an object of class \code{\link{HDA}}, with an additional dataset called "\code{SC}" (small colony).
#'This is a logical matrix, with code{FALSE} representing strains which have been excluded.
excludeSmallColonies <- function(HDA, what="auto", exclude="all", method="otsu", manualThreshold=NULL) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude, "SC")
  
  if(class(method)!="character") {
    stop("Class of 'method' must be character.")
  }
  if(method %in% c("otsu", "manual")==F) {
    stop("Supported threshold methods are 'otsu' or 'manual.")
  }
  if(is.null(control(HDA))) {
    stop("Small colonies can only be excluded if a control has been specified. See ?`control<-`.")
  }
  
  #Extract the dataset
  sizes <- extractAssay(HDA, what, exclude)
  
  #Create a vector of median colony sizes for the control(s)
  controlMedians <- stats::aggregate(sizes[,which(colnames(sizes) %in% control(HDA)), drop=F], list(replicate=rowData(HDA)$replicate), median)
  
  #Threshold based on method
  if(method=="otsu") {
    threshold <- apply(controlMedians[ ,-1, drop=F], 2, otsuThreshold)
  }
  if(method == "manual") {
    if(is.numeric(manualThreshold)==F) {
      stop("For manual thresholding, please enter a numeric value for 'manualThreshold'.")
    }
    threshold <- rep_len(manualThreshold, ncol(controlMedians)-1)
    names(threshold) <- colnames(controlMedians)[2:ncol(controlMedians)]
  }
  
  #Save the threshold
  thresholdSaved <- threshold
  
  #Extract the including NAs
  sizes <- extractAssay(HDA, what, NULL)
  
  #Create a vector of median colony sizes for the control(s)
  controlMedians <- stats::aggregate(sizes[,which(colnames(sizes) %in% control(HDA)), drop=F], list(replicate=rowData(HDA)$replicate), median)
  
  #Subset for just colony sizes
  controlMediansSub <- controlMedians[, -1, drop=F]
  
  #Threshold
  threshold <- matrix(rep(threshold, nrow(controlMediansSub)), ncol=ncol(controlMediansSub), nrow=nrow(controlMediansSub), byrow = T)
  threshold <- ifelse(controlMediansSub < threshold, F, T)
  
  #Match back to the library file
  threshold <- cbind(controlMedians[, 1, drop=F], as.data.frame(threshold))
  exclusions <- threshold[match(rowData(HDA)$replicate, controlMedians$replicate),]
  
  #Create a matrix showing which to keep
  exclusions <- exclusions[, -1, drop=F]
  controls <- colData(HDA)$batch[match(colnames(exclusions), rownames(colData(HDA)))]
  exclusions <- exclusions[, match(colData(HDA)$batch, controls), drop=F]
  colnames(exclusions) <- rownames(colData(HDA))
  
  #Add to the HDA object
  exclusions <- as.matrix(exclusions)
  rownames(exclusions) <- NULL
  assays(HDA)$SC <- exclusions
  
  #Mark SC as a quality control dataset
  qualityControl(HDA) <- c(qualityControl(HDA), "SC")
  
  #Add to the pipeline
  pipeline(HDA)$SC <- SimpleList(threshold=thresholdSaved)
  
  #Return
  return(HDA)
  
}