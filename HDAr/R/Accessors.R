#'@include HDA.R
#'@export
#'@title HDA Metadata Accessors
#'@description Utility functions to access or set the value of slots of a \code{\link{HDAMetadata}} object.
setGeneric("control", function(object) {standardGeneric("control")})
setMethod("control", "list", function(object) {object$control})
setMethod("control", "HDA", function(object) {
  object <- metadata(object)
  callGeneric()
})

#'@export
#'@rdname control
setGeneric("control<-", function(object, value) {standardGeneric("control<-")})
setMethod("control<-", "list", function(object, value) {
  object <- AddHDAMetadata(object, control=value)
})
setMethod("control<-", signature(object="HDA", value="vector"), function(object, value) {
  #Do not use `metadata<-` as object will fail `validObject` unless `colData` is changed
  object@metadata <- AddHDAMetadata(metadata(object), control=value)
  colData(object)$control <- as.logical(ifelse(rownames(colData(object)) %in% value, 1, 0))
  return(object)
})
setMethod("control<-", signature(object="HDA", value="NULL"), function(object, value) {
  #Do not use `metadata<-` as object will fail `validObject` unless `colData` is changed
  object@metadata <- AddHDAMetadata(metadata(object), control=value)
  colData(object)$control <- NULL
  return(object)
})

#'@export
#'@rdname control
setGeneric("batch", function(object) {standardGeneric("batch")})
setMethod("batch", "HDA", function(object) {
  colData(object)$batch
})

#'@export
#'@rdname control
setGeneric("batch<-", function(object, value) {standardGeneric("batch<-")})
setMethod("batch<-", "HDA", function(object, value) {
  colData(object)$batch <- as.factor(value)
  return(object)
})

#'@export
#'@rdname control
setGeneric("activeDataset", function(object) {standardGeneric("activeDataset")})
setMethod("activeDataset", "list", function(object) {object$activeDataset})
setMethod("activeDataset", "HDA", function(object) {
  object <- metadata(object)
  callGeneric()
})

#'@export
#'@rdname control
setGeneric("activeDataset<-", function(object, value) {standardGeneric("activeDataset<-")})
setMethod("activeDataset<-", "list", function(object, value) {
  AddHDAMetadata(object, activeDataset=value)
})
setMethod("activeDataset<-", "HDA", function(object, value) {
  metadata(object) <- AddHDAMetadata(metadata(object), activeDataset=value)
  return(object)
})

#'@export
#'@rdname control
setGeneric("qualityControl", function(object) {standardGeneric("qualityControl")})
setMethod("qualityControl", "list", function(object) {object$qualityControl})
setMethod("qualityControl", "HDA", function(object) {
  object <- metadata(object)
  callGeneric()
})

#'@export
#'@rdname control
setGeneric("qualityControl<-", function(object, value) {standardGeneric("qualityControl<-")})
setMethod("qualityControl<-", "list", function(object, value) {
  AddHDAMetadata(object, qualityControl=value)
})
setMethod("qualityControl<-", "HDA", function(object, value) {
  metadata(object) <- AddHDAMetadata(metadata(object), qualityControl=value)
  return(object)
})

setMethod("organism", "list", function(object) {object$organism})
setMethod("organism", "HDA", function(object) {
  object <- metadata(object)
  callGeneric()
})

setMethod("organism<-", "list", function(object, value) {
  AddHDAMetadata(object, organism=value)
})
setMethod("organism<-", "HDA", function(object, value) {
  metadata(object) <- AddHDAMetadata(metadata(object), organism=value)
  return(object)
})

setMethod("metadata<-", signature(x="HDA", value="list"), function(x, ..., value) {
  x <- callNextMethod()
  validObject(x)
  return(x)
})

setMethod("rowData<-", signature(x="HDA", value="DataFrame"), function(x, ..., value) {
  x <- callNextMethod()
  validObject(x)
  return(x)
})
setMethod("colData<-", signature(x="HDA", value="DataFrame"), function(x, ..., value) {
  x <- callNextMethod()
  validObject(x)
  return(x)
})

setMethod("rowRanges", "HDA", function(x, ...) {
  return(x@rowRanges)
})
setMethod("rowRanges<-", signature(x="HDA", value="GRanges"), function(x, ..., value) {
  x@rowRanges <- value
  validObject(x)
  return(x)
})
setMethod("rowRanges<-", signature(x="HDA", value="NULL"), function(x, ..., value) {
  x@rowRanges <- NULL
  validObject(x)
  return(x)
})

#'@export
setGeneric("colRanges", function(x) {standardGeneric("colRanges")})
setMethod("colRanges", "HDA", function(x) {
  x@colRanges
})

#'@export
setGeneric("colRanges<-", function(x, value) {standardGeneric("colRanges<-")})
setMethod("colRanges<-", signature(x="HDA", value="GRangesList"), function(x, value) {
  x@colRanges <- value
  if(is.null(names(x@colRanges))) {
    names(x@colRanges) <- colnames(x)
  }
  validObject(x)
  return(x)
})
setMethod("colRanges<-", signature(x="HDA", value="ANY"), function(x, value) {
  if(is.null(organism(x))) {
    stop("Unable to automatically convert gene names to GRanges if the organism is not specified. See ?'organism<-'.")
  }
  value <- convertHDAColList2GRangesList(HDAColList(value), organism(x))
  callGeneric()
})

setMethod("colnames<-", signature(x="HDA"), function(x, value) {
  #Do not use accessors as object will fail `validObject` unless `colData` and `colRanges` are changed
  coldata <- colData(x)
  coldata$sample <- value
  x@colData <- coldata
  colranges <- colRanges(x)
  names(colranges) <- value
  x@colRanges <- colranges
  x <- callNextMethod()
  validObject(x)
  return(x)
})

setMethod("rownames<-", signature(x="HDA"), function(x, value) {
  x <- callNextMethod()
  names(rowRanges(x)) <- value
  validObject(x)
  return(x)
})

#'@export
#'@rdname control
setGeneric("pipeline", function(object) {standardGeneric("pipeline")})
setMethod("pipeline", "list", function(object) {object$pipeline})
setMethod("pipeline", "HDA", function(object) {
  object <- metadata(object)
  callGeneric()
})

#'@export
#'@rdname control
setGeneric("pipeline<-", function(object, value) {standardGeneric("pipeline<-")})
setMethod("pipeline<-", "list", function(object, value) {
  AddHDAMetadata(object, pipeline=value)
})
setMethod("pipeline<-", "HDA", function(object, value) {
  metadata(object) <- AddHDAMetadata(metadata(object), pipeline=value)
  return(object)
})
