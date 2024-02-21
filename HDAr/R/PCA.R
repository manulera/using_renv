#'@export
plotPCA <- function(HDA, what="auto", exclude="all", impute = T, maxMissing=0.5, replicate="separate", x=1, y=2, ...) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  if(class(replicate)!="character") {
    stop("Class of 'replicate' must be character.")
  }
  if(replicate %in% c("separate", "median", "mean", "none")==F) {
    stop("Supported methods for dealing with replicates are 'separate', 'median', 'mean' and 'none'.")
  }
  
  #Remove footprints and order by replicate
  if(any(is.na(rowData(HDA)$replicate))) {
    HDA <- HDA[-which(is.na(rowData(HDA)$replicate)),]
  }
  HDA <- HDA[order(rowData(HDA)$replicate, rowData(HDA)$plate, rowData(HDA)$row, rowData(HDA)$col),]
  
  #Get data
  sizes <- extractAssay(HDA, what, exclude)

  #Get coldata
  coldata <- colData(HDA)
  
  #Decide what to do with replicates
  if(replicate=="separate") {
    nReplicates <- unique(table(rownames(sizes)))
    if(length(nReplicates)!=1) {
      stop("The number of technical replicates for each mutant is not equal. Unable to treat each replicate separately.")
    }
    indicies <- 1:nReplicates%%nReplicates
    sizes <- do.call("cbind", lapply(1:nReplicates, function(i) {
      sizes[1:nrow(sizes)%%nReplicates==indicies[i], , drop=F]
    }))
    coldata <- do.call("rbind", rep(list(coldata), nReplicates))
    coldata$replicate <- as.factor(rep(1:nReplicates, each=ncol(HDA)))
  }
  if(replicate=="median") {
    sizes <- stats::aggregate(sizes, list(replicate=rowData(HDA)$replicate), median)
    rownames(sizes) <- sizes$replicate
    sizes$replicate <- NULL
    sizes <- as.matrix(sizes)
  }
  if(replicate=="mean") {
    sizes <- stats::aggregate(sizes, list(replicate=rowData(HDA)$replicate), mean)
    rownames(sizes) <- sizes$replicate
    sizes$replicate <- NULL
    sizes <- as.matrix(sizes)
  }
  
  #Remove data for which too much is missing
  nMissing <- apply(sizes, 1, function(i) {sum(is.na(i))})
  remove <- ifelse(nMissing/ncol(sizes)>maxMissing, T, F)
  if(any(remove)) {
    sizes <- sizes[!remove,]
  }
  
  #Deal with missing values
  if(impute) {
    
    invisible(capture.output({sizes <- impute::impute.knn(sizes)$data}))
    
  } else {
    
    ind2remove <- which(apply(sizes, 1, function(i) {any(is.na(i))}))
    
    if(length(ind2remove)>0) {
      
      sizes <- sizes[-ind2remove,]
      
    }
    
  }
  
  #Transpose
  sizes <- t(sizes)
  
  #Perform PCA
  pca <- prcomp(sizes)
  
  #Plot
  g <- ggfortify:::autoplot.prcomp(pca, x=x, y=y, data=as.data.frame(coldata))
  
  #Add additional aesthetics
  args <- as.list(sapply(match.call()[-1], deparse))
  args <- args[-which(names(args) %in% c("HDA", "what", "exclude", "impute", "maxMissing", "replicate", "x", "y"))]
  if(length(args)>0) {
    g <- gginnards::delete_layers(g, "GeomPoint")
    aes_string <- ggplot2::aes_string
    g <- g + ggplot2::geom_point(do.call("aes_string", args))
  }
  
  return(g)
  
}

#'@export
batchNormalise <- function(HDA, what="auto", exclude="all", replicate=T, model=model.matrix(~1, data=colData(HDA))) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Check that there are batches
  if(is.null(batch(HDA))) {
    stop("No batches have been indicated for the HDA object. Unable to correct for batch effetcs.")
  }
  
  #Order by replicate
  HDA <- HDA[order(rowData(HDA)$replicate, rowData(HDA)$plate, rowData(HDA)$row, rowData(HDA)$col),]
  
  #Get data
  sizes <- extractAssay(HDA, what, exclude)

  #Get batch
  batch1 <- colData(HDA)$batch
  
  #Decide what to do with replicates
  if(replicate) {
    nReplicates <- unique(table(rownames(sizes)))
    if(length(nReplicates)!=1) {
      stop("The number of technical replicates for each mutant is not equal. Unable to treat each replicate separately.")
    }
    indicies <- 1:nReplicates%%nReplicates
    sizes <- do.call("cbind", lapply(1:nReplicates, function(i) {
      sizes[1:nrow(sizes)%%nReplicates==indicies[i], , drop=F]
    }))
    model <- do.call("rbind", rep(list(model), nReplicates))
    batch1 <- rep(batch1, nReplicates)
    batch2 <- as.factor(rep(1:nReplicates, each=ncol(HDA)))
  } else {
    batch2 <- NULL
  }
  
  #Deal with missing values
  invisible(capture.output({sizes <- impute::impute.knn(sizes)$data}))
  
  #Perform batch correction
  sizes <- limma::removeBatchEffect(sizes, batch1, batch2, design=model)
  
  #Put things replicates back in the same columns
  if(replicate) {
    sizes2 <- matrix(0, nrow=nrow(HDA), ncol=ncol(HDA))
    for(i in 1:nReplicates) {
      rows1 <- 1:nrow(sizes)
      cols1 <- seq(from=(ncol(HDA)*(i-1)+1), to=ncol(HDA)*i)
      rows2 <- (1:nrow(sizes2))[(1:nrow(sizes2))%%nReplicates==indicies[i]]
      cols2 <- 1:ncol(sizes2)
      sizes2[rows2, cols2] <- sizes[rows1, cols1]
    }
  } else {
    sizes2 <- sizes
    rownames(sizes2) <- NULL
  }
  
  #Put back in HDA
  assays(HDA)$BNS <- sizes2
  
  return(HDA)
  
}
