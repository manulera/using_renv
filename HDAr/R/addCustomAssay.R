#'@export
addCustomAssay <- function(HDA, assay, assayName, qualityControl) {
  
  assays(HDA)[[assayName]] <- assay
  
  if(qualityControl) {
    metadata(HDA)$qualityControl <- unique(c(metadata(HDA)$qualityControl, assayName))
  }
  
  return(HDA)
  
}