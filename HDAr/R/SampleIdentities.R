setClassUnion("character_Or_NULL", c("character", "NULL"))

validHDAColList <- function(object) {
  
  #Create a vector to store error messages
  error <- character(0)
  
  #Check that the class of every element of the list is 'character'.
  if(all(sapply(object, is, "character_Or_NULL"))==F) {
    msg <- "The class of every element in a 'HDAColList' must be 'character' or 'NULL'."
    error <- c(error, msg)
  }
  
  #Return any errors if necessary
  if(length(error)==0) {return(T)} else {return(error)}
  
}

setClass("HDAColList", contains = "SimpleList", validity = validHDAColList)

setGeneric("HDAColList", function(conditions) {standardGeneric("HDAColList")})
setMethod("HDAColList", "character", function(conditions) {
  conditions <- as.list(conditions)
  callGeneric()
})
setMethod("HDAColList", "list", function(conditions) {
  conditions <- SimpleList(conditions)
  callGeneric()
})
setMethod("HDAColList", "SimpleList", function(conditions) {
  new("HDAColList", conditions)
})

convertHDAColList2GRangesList <- function(HDAColList, organism) {
  
  genes <- unique(unlist(HDAColList))
  
  genePositions <- downloadGenes(genes, organism)

  missedGenes <- genes[genes %in% genePositions$ensembl_gene_id == F]
  
  if(length(missedGenes)>0) {
    stop(paste0(
      "Unable to find positional information for the following gene(s): ",
      paste0(missedGenes, collapse=", "),
      ". Please check that systematic gene names have been specified correctly, or manually construct the 'GRangesList' object."
    ))
  }
  
  GRangesList(lapply(HDAColList, function(x) {
    if(is.null(x)) {
      return(GenomicRanges::GRanges())
    }
    indicies <- match(x, genePositions$ensembl_gene_id)
    GRanges(Rle(genePositions$chromosome_name[indicies]), IRanges(genePositions$start_position[indicies], genePositions$end_position[indicies]), Rle(genePositions$strand[indicies]), gene=genePositions$ensembl_gene_id[indicies])
  }))

}

