getTotalGRanges <- function(HDA) {
  
  #Get all chromosomes
  totalGRanges <- range(rowRanges(HDA), ignore.strand=T)
  
  #Remove any chromosomes called 'NA'
  if(any(as.vector(seqnames(totalGRanges))=="NA")) {
    NAs <- as.vector(seqnames(totalGRanges)=="NA")
    totalGRanges <- totalGRanges[-which(NAs)]
  }
  
  return(totalGRanges)
  
}

prepareLinkageData <- function(HDA, what, exclude) {
  
  #Get the data and aggregate medians
  sizes <- extractAssay(HDA, what, exclude)
  sizes <- aggregate(sizes, by=list(replicate=rowData(HDA)$replicate), median, na.rm=T)
  
  #Create a data frame of positional information for the library and linkage data and aggregate to 1 row per technical replicate
  positions <- data.frame(start=start(rowRanges(HDA)), chr=seqnames(rowRanges(HDA)))
  positions <- aggregate(positions, list(replicate=rowData(HDA)$replicate), function(i) {i[1]})
  
  #Return
  return(list(sizes=sizes, positions=positions))
  
}

#A functin to fit a median spline to size data for a chromosome for a sample (i and j)
fitMedianSplineToChr <- function(HDA, positions, sizes, totalGRanges, nPoints, i, j, ...) {
  
  chr <- as.character(seqnames(totalGRanges)[j])
  
  df <- data.frame(start=positions$start[positions$chr==chr], size=sizes[positions$chr==chr, i])
  df <- df[order(df$start),]
  
  invisible(capture.output({fit <- tryCatch({
    suppressWarnings({cobs::cobs(df$start, df$size, ...)})
  }, error=function(e) {
    msg <- paste0(
      "Unable to fit a spline to chromosome '",
      chr,
      "' for sample '",
      i,
      "'."
    )
    warning(msg, call. = F)
    return(NULL)
  })}))
  if(is.null(fit)) {
    return(NULL)
  }
  
  pos <- seq(from=start(totalGRanges[j]), to=end(totalGRanges[j]), length.out=nPoints[j])
  
  fit <- predict(fit, pos)
  
  return(fit)
  
}

fitMedianSplineToAllChr <- function(HDA, positions, sizes, i, totalGRanges, ...) {
  
  #Set the number of points to be used to fit a spline for each chromosome
  chrSizes <- abs(end(totalGRanges) - start(totalGRanges))
  genomeSize <- sum(chrSizes)
  nLoci <- nrow(sizes)
  nPoints <- as.integer((nLoci/genomeSize)*chrSizes)
  
  #Fit a spline to each chromosome
  fits <- lapply(1:length(totalGRanges), function(j) {
    fitMedianSplineToChr(HDA, positions, sizes, totalGRanges, nPoints, i, j, ...)
  })
  
  #Return
  return(fits)
  
}

#'@export
excludeLinkedLoci <- function(HDA, distance=250000) {
  
  #Check arguments
  checkHDA(HDA)
  checkRowRanges(HDA)
  
  #Get all chromosomes
  totalGRanges <- getTotalGRanges(HDA)
  
  #For each SGA
  exclusionData <- lapply(1:ncol(HDA), function(i) {
    
    #Get the positions of the query loci
    colRange <- colRanges(HDA)[[i]]
    
    #If there are no colRanges for this gene, then skip
    if(length(colRange)==0) {
      
      #Send a warning that no linkage exclusion was performed
      warning(paste0(
        "There are no 'colRanges' indicated for '",
        colnames(HDA)[i],
        "'. Unable to exclude linked loci for this sample."
      ), call. = F)
      
      #Set everything as unlinked
      mat <- rep(T, nrow(HDA))
      linkage <- totalGRanges
      linkage$linked <- F
      
      #Return
      return(list(mat, linkage))
      
    }
    
    #Use only start
    end(colRange) <- start(colRange)
    
    #Expand the query loci 
    start(colRange) <- start(colRange) - distance
    end(colRange) <- end(colRange) + distance
    
    #Reduce to a single GRange
    colRange <- GenomicRanges::reduce(colRange)
    
    #Find the overlaps with the library mutants
    hits <- GenomicRanges::findOverlaps(colRange, rowRanges(HDA), ignore.strand=T)
    
    #Mark the overlaps
    mat <- rep(T, nrow(HDA))
    mat[is.na(rowData(HDA)$library_mutant)] <- NA
    mat[subjectHits(hits)] <- F
    
    #Create a GRanges object showing the regions which are unlinked
    hits <- GenomicRanges::findOverlaps(totalGRanges, colRange, ignore.strand=T)
    unlinked <- unlist(GenomicRanges::psetdiff(totalGRanges, IRanges::extractList(colRange, as(hits, "List"))))
    unlinked$linked <- F
    
    #Create a GRanges object showing all chromosomal regions are linked or unlinked
    colRange$gene <- NULL
    colRange$linked <- T
    linkage <- c(unlinked, colRange)
    strand(linkage) <- "*"
    
    return(list(mat, linkage))
    
  })
  
  exclusions <- do.call("cbind", lapply(exclusionData, `[[`, 1))
  linkage <- lapply(exclusionData, `[[`, 2)
  names(linkage) <- colnames(HDA)
  
  assays(HDA)$L <- exclusions
  
  qualityControl(HDA) <- c(qualityControl(HDA), "L")
  
  pipeline(HDA)$L <- SimpleList(linkage=SimpleList(linkage))
  
  return(HDA)
  
}
