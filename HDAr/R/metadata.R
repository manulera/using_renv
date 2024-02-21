AddHDAMetadata <- function(metadata, ...) {
  newMetadata <- list(...)
  if(is.null(names(newMetadata)) | any(names(newMetadata)=="")) {
    stop("Please supply named arguments to when adding metadata.")
  }
  duplicatedElements <- names(newMetadata)[names(newMetadata) %in% names(metadata)]
  if(any(names(metadata) %in% duplicatedElements)) {
    metadata <- metadata[-which(names(metadata) %in% duplicatedElements)]
  }
  allMetadata <- c(metadata, newMetadata)
  allMetadata[which(sapply(allMetadata, function(x) {is.null(x)==F}))]
}
