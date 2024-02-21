#'@export
plotSmallColonyExclusion <- function(HDA, what="auto", exclude="all", method="otsu", manualThreshold=NULL, facet=T) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude, "SC")
  
  #Perform exclusion
  HDA <- excludeSmallColonies(HDA, what, exclude, method, manualThreshold)
  
  #Extract data according to the pipeline
  sizes <- extractAssay(HDA, what, exclude)
  
  #Subset for controls
  sizes <- sizes[,colnames(sizes) %in% names(pipeline(HDA)$SC$threshold), drop=F]
  
  #Calculate median colony sizes of each technical replicate
  sizes <- stats::aggregate(sizes, list(replicate=rowData(HDA)$replicate), median)
  sizes <- sizes[, -1, drop=F]
  
  #Convert to long format
  sizes <- reshape2::melt(sizes, id.vars=character(0), variable.name="Sample")
  sizes$Sample <- factor(sizes$Sample, levels=unique(sizes$Sample))
  
  #Turn threshold into a data.frame
  threshold <- data.frame(Sample=names(pipeline(HDA)$SC$threshold), threshold=unname(pipeline(HDA)$SC$threshold))
  threshold$Sample <- factor(threshold$Sample, levels=threshold$Sample)
  
  #Plot
  g <- ggplot2::ggplot(sizes, ggplot2::aes(x=value, fill=Sample))
  if(nrow(threshold)==1) {
    g <- g + ggplot2::geom_histogram(bins = nrow(sizes)/nrow(threshold)/50, colour="white")
  } else {
    g <- g + ggplot2::geom_histogram(position="dodge", bins = nrow(sizes)/nrow(threshold)/50)
  }
  g <- g + ggplot2::geom_vline(data=threshold, ggplot2::aes(xintercept=threshold, colour=Sample), linetype="dashed")
  if(facet) {
    g <- g + ggplot2::facet_grid(Sample~.)
  }
  
  #Return
  return(g)
  
}

#Function to create a plot of a plate
plateHeatmap <- function(sizeData, median, excludedColour) {
  
  colnames(sizeData)[1:2] <- c("ymax", "xmax")
  sizeData$ymin <- sizeData$ymax - 1
  sizeData$xmin <- sizeData$xmax - 1
  
  #Map row and col data
  g <- ggplot2::ggplot(sizeData, ggplot2::aes(ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax))
  
  #Add z - colony size
  g <- g + ggplot2::geom_rect(ggplot2::aes_string(fill=colnames(sizeData)[3]))
  
  #Put minimal theme
  g <- g + ggplot2::theme_minimal()
  
  #Fix the aspect ratio of the axes to 1
  g <- g + ggplot2::coord_fixed()
  
  #Create a colour gradient from blue to black to yellow
  myPalette <- colorRampPalette(c("dodgerblue", "dodgerblue4", "black", "gold4", "gold1"))
  
  #Apply the gradient such that blue is 0, black is the median colony size, and yellow is twice the median or greater
  #Also set NAs as a different colour if desired
  if(is.null(excludedColour)) {
    g <- g + ggplot2::scale_fill_gradientn(colours=myPalette(100), limits=c(0,median*2), oob=scales::squish)
  } else {
    g <- g + ggplot2::scale_fill_gradientn(colours=myPalette(100), limits=c(0,median*2), oob=scales::squish, na.value = excludedColour)
  }
  
  #Return the plot
  return(g)
  
}

#Function to check that HDA has been subsetted to a single sample and single plate
HDAplatecheck <- function(HDA) {
  
  if(ncol(HDA)!=1) {
    
    stop("Please select a single sample to plot.")
    
  }
  
  if(length(unique(rowData(HDA)$plate))!=1) {
    
    stop("Please select a single plate to plot.")
    
  }
  
}

#'@export
plotPlate <- function(HDA, what="auto", exclude="all", sample, plate, footprints=T, excludedColour=NULL) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Max dims
  maxDims <- c(max(rowData(HDA)$row), max(rowData(HDA)$col))
  
  #Subset for plate and sample
  control(HDA) <- NULL
  HDA <- HDA[rowData(HDA)$plate==plate,]
  HDA <- HDA[,sample]
  
  #Check that a single plate has been selected
  HDAplatecheck(HDA)
  
  #Get the colony sizes
  if(is.null(excludedColour)==F) {
    sizes <- extractAssay(HDA, what, exclude)
  } else {
    sizes <- extractAssay(HDA, what, NULL)
  }
  colnames(sizes) <- what
  
  #Get row/col data
  sizes <- as.data.frame(cbind(DataFrame(rowData(HDA))[, c("row", "col")], sizes))
  
  #Plot the plate
  g <- plateHeatmap(sizes, median(extractAssay(HDA, what, exclude)[,sample], na.rm=T), excludedColour)
  
  #Add labs and limits
  g <- g + ggplot2::xlab("col") + ggplot2::ylab("row")
  g <- g + ggplot2::xlim(c(0, maxDims[2])) + ggplot2::scale_y_reverse(limits=c(maxDims[1], 0))
  
  #Add footprints if desired
  if(footprints) {
    footprints <- as.data.frame(DataFrame(rowData(HDA))[is.na(rowData(HDA)$replicate),c("row", "col")])
    if(nrow(footprints)>1) {
      g <- g + ggplot2::geom_rect(data=footprints, ggplot2::aes(xmin=col-1, xmax=col, ymin=row-1, ymax=row), colour="white", fill="transparent")
    }
  }
  
  #Return the plot
  return(g)
  
}

#'@export
plotNormalisationSurface <- function(HDA, what="auto", exclude="all", sample, plate, method=c("median", "rowcol", "spatial"), lowess=T, span=5) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Max dims
  maxDims <- c(max(rowData(HDA)$row), max(rowData(HDA)$col))
  
  #Subset for plate and sample
  control(HDA) <- NULL
  HDA <- HDA[rowData(HDA)$plate==plate,]
  HDA <- HDA[,sample]
  
  #Check that a single plate has been selected
  HDAplatecheck(HDA)
  
  #Construct the normalisation surface
  normalisationSurface <- constructNormalisationSurface(HDA, what, exclude, method, lowess, span)
  colnames(normalisationSurface) <- what
  
  #Get row/col data
  normalisationSurface <- as.data.frame(cbind(DataFrame(rowData(HDA))[, c("row", "col")], normalisationSurface))
  
  #Plot the plate
  g <- plateHeatmap(normalisationSurface, median(extractAssay(HDA, what, exclude)[,sample], na.rm=T), NULL)
  
  #Add labs and limits
  g <- g + ggplot2::xlab("col") + ggplot2::ylab("row")
  g <- g + ggplot2::xlim(c(0, maxDims[2])) + ggplot2::scale_y_reverse(limits=c(maxDims[1], 0))
  
  return(g)
  
}

#'@export
plotPlateMultiple <- function(HDA, what="auto", exclude="all", sample, plate, footprints=T, excludedColour=NULL) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Get arguments for sample and plate
  if(missing(sample)) {
    sample <- colnames(HDA)
  }
  if(missing(plate)) {
    plate <- unique(rowData(HDA)$plate)
  }
  
  #Create a list to save plots
  plotList <- list()
  
  #Create each plot
  counter <- 0
  for(j in 1:length(plate)) {
    for(i in 1:length(sample)) {
      
      counter <- counter + 1
      
      plotList[[counter]] <- plotPlate(HDA, what, exclude, sample[i], plate[j], footprints, excludedColour)
      
      #Add row/column titles
      if(i==1 & j==1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top=sample[i], left=paste0("Plate ", plate[j]))
      }
      if(i==1 & j>1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top="", left=paste0("Plate ", plate[j]))
      }
      if(i>1 & j==1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top=sample[i], left="")
      }
      if(i>1 & j>1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top="", left="")
      }
      
    }
  }
  
  #Plot
  plotList$ncol <- length(sample)
  grid.arrange <- gridExtra::grid.arrange
  do.call("grid.arrange", plotList)
  
}

#'@export
plotNormalisationSurfaceMultiple <- function(HDA, what="auto", exclude="all", sample, plate, method=c("median", "rowcol", "spatial"), lowess=T, span=5) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Get arguments for sample and plate
  if(missing(sample)) {
    sample <- colnames(HDA)
  }
  if(missing(plate)) {
    plate <- unique(rowData(HDA)$plate)
  }
  
  #Create a list to save plots
  plotList <- list()
  
  #Create each plot
  counter <- 0
  for(j in 1:length(plate)) {
    for(i in 1:length(sample)) {
      
      counter <- counter + 1
      
      plotList[[counter]] <- plotNormalisationSurface(HDA, what, exclude, sample[i], plate[j], method, lowess, span)
      
      #Add row/column titles
      if(i==1 & j==1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top=sample[i], left=paste0("Plate ", plate[j]))
      }
      if(i==1 & j>1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top="", left=paste0("Plate ", plate[j]))
      }
      if(i>1 & j==1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top=sample[i], left="")
      }
      if(i>1 & j>1) {
        plotList[[counter]] <- gridExtra::arrangeGrob(plotList[[counter]], top="", left="")
      }
      
    }
  }
  
  #Plot
  plotList$ncol <- length(sample)
  grid.arrange <- gridExtra::grid.arrange
  do.call("grid.arrange", plotList)
  
}

#'@export
plotLinkage <- function(HDA, what="auto", exclude="all", sample, distance=250000, ...) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude, "L")
  checkRowRanges(HDA)
  
  #Perform exclusion
  control(HDA) <- NULL
  HDA <- HDA[, sample]
  HDA <- excludeLinkedLoci(HDA, distance)
  
  #Create limits for y axis
  assay <- extractAssay(HDA, what, exclude)
  assay <- assay[, sample]
  limits <- c(0, max(c(2*(median(assay, na.rm=T)), quantile(assay, 0.995, na.rm=T))))
  
  #Get all chromosomes
  totalGRanges <- getTotalGRanges(HDA)
  
  #Aggregate linkage data
  linkageData <- prepareLinkageData(HDA, what, exclude)
  sizes <- linkageData$sizes
  positions <- linkageData$positions
  
  #Fit a median spline to each chromosome
  fits <- fitMedianSplineToAllChr(HDA, positions, sizes, sample, totalGRanges, ...)
  
  #Convert the fits in GRanges and mark everything as unlinked
  fits <- lapply(1:length(totalGRanges), function(j) {
    if(is.null(fits[[j]])) {
      return(NULL)
    }
    fit <- fits[[j]]
    fit <- GenomicRanges::GRanges(seqnames(totalGRanges[j]), IRanges::IRanges(fit[,"z"], fit[,"z"]), "*", fit=fit[,"fit"])
    fit$linkage <- "3"
    return(fit)
  })
  
  #Get the linkage information and mark parts of the fit within the linked region as linked
  linkage <- pipeline(HDA)$L$linkage[[sample]]
  linkageLinked <- linkage[linkage$linked]
  fits <- lapply(fits, function(fit) {
    if(is.null(fit)) {
      return(NULL)
    }
    fit$linkage[subjectHits(GenomicRanges::findOverlaps(linkageLinked, fit))] <- "2"
    return(fit)
  })
  
  #Convert positions into a GRanges object containing size and linkage information for the selected sample
  positions <- GenomicRanges::GRanges(Rle(positions$chr), IRanges::IRanges(positions$start, positions$start), "*", sample=sizes[,sample])
  
  #Mark the linked genes as linked and split each fit based on the linakge classification
  linkageClassification <- extractAssay(HDA, "L", exclude)[,sample]
  linkageClassification <- aggregate(linkageClassification, list(replicate=rowData(HDA)$replicate), function(i) {i[1]})$x
  positions$linkage <- linkageClassification
  fits <- do.call("c", fits)
  if(is.null(fits)==F) {
    fits <- lapply(1:length(linkage), function(i) {
      fit <- fits[subjectHits(GenomicRanges::findOverlaps(linkage[i], fits))]
    })
  }
  positions$linkage <- factor(as.numeric(positions$linkage), levels=as.character(0:3))
  
  #Plot gene start against size data, colour based on linkage classification
  g <- ggplot2::ggplot(as.data.frame(positions), ggplot2::aes(start, sample, colour=linkage))
  g <- g + ggplot2::geom_point()
  g <- g + ggplot2::guides(colour="none")
  g <- g + ggplot2::xlab("Start / bp")
  g <- g + ggplot2::ylab(what)
  
  #Split for each chromosome
  g <- g + ggplot2::facet_grid(seqnames~.)
  
  #Set manual colours based on linkage classification
  cols <- c("firebrick4", "dodgerblue4", "firebrick1", "dodgerblue1")
  names(cols) <- as.character(0:3)
  g <- g + ggplot2::scale_colour_manual(values=cols)
  
  #Get chromosome positions for each query and plot straight line to show location
  lociData <- as.data.frame(colRanges(HDA)[[sample]])
  if(nrow(lociData)>0) {
    lociData$linkage <- "2"
    g <- g + ggplot2::geom_vline(data=lociData, ggplot2::aes(xintercept=start, colour=linkage),  linetype="longdash")
  }
  
  #Construct a line for each fit and add to the plot
  lines <- lapply(1:length(fits), function(i) {
    if(is.null(fits[[i]])) {
      return(NULL)
    }
    ggplot2::geom_line(data=as.data.frame(fits[[i]]), ggplot2::aes(start, fit, colour=linkage))
  })
  for(i in 1:length(lines)) {
    if(is.null(lines[[i]])) {
      next
    }
    g <- g + lines[[i]]
  }
  
  #Add limits
  g <- g + ggplot2::ylim(limits)
  
  #Return the plot
  return(g)
  
}

#'@export
plotLinkageMultiple <- function(HDA, what="auto", exclude="all", sample, distance=250000, ...) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Get arguments for sample
  if(missing(sample)) {
    sample <- colnames(HDA)
  }
  
  #Create limits for y axis
  assay <- extractAssay(HDA, what, exclude)
  assay <- assay[, sample]
  limits <- c(0, max(c(2*(median(assay, na.rm=T)), quantile(assay, 0.995, na.rm=T))))
  
  #Create a list to save plots
  plotList <- list()
  
  #Create each plot
  for(i in 1:length(sample)) {
    
    plotList[[i]] <- plotLinkage(HDA, what, exclude, sample[i], distance, ...)
    suppressMessages({
      plotList[[i]] <- plotList[[i]] + ggplot2::ylim(limits)
    })
    
    #Add column titles
    plotList[[i]] <- gridExtra::arrangeGrob(plotList[[i]], top=sample[i])
    
  }
  
  #Plot
  plotList$ncol <- length(sample)
  grid.arrange <- gridExtra::grid.arrange
  return(do.call("grid.arrange", plotList))
  
}
