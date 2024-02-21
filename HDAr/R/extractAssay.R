#'@export
extractAssay <- function(HDA, assay, exclude) {
  
  assay <- assays(HDA)[[assay]]
  
  qcMat <- excludedData(HDA, exclude)
  
  assay[qcMat==F] <- NA
  
  rownames(assay) <- rowData(HDA)$replicate
  
  return(assay)
  
}
