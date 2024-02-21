downloadGenes <- function(genes, organism) {
  
  attributesToDownload <- c("chromosome_name", "start_position", "end_position", "strand")
  
  if(organism=="Saccharomyces cerevisiae") {
    mart <- "ENSEMBL_MART_ENSEMBL"
    dataset <- "scerevisiae_gene_ensembl"
  }
  if(organism=="Schizosaccharomyces pombe") {
    mart <- "fungi_mart"
    dataset <- "spombe_eg_gene"
  }
  
  genePositions <- biomartr::biomart(genes, mart, dataset, attributes=attributesToDownload, "ensembl_gene_id")
  
  genePositions$chromosome_name <- as.factor(genePositions$chromosome_name)
  genePositions$strand <- as.character(genePositions$strand)
  genePositions$strand[genePositions$strand=="1"] <- "+"
  genePositions$strand[genePositions$strand=="-1"] <- "-"
  genePositions$strand[genePositions$strand %in% c("+", "-") ==F] <- "*"
  genePositions$strand <- as.factor(genePositions$strand)
  levels(genePositions$strand) <- c("+", "-", "*")
  
  return(genePositions)
  
}

#'@export
addMutantPositionalInformation <- function(HDA) {
  
  if(is.null(organism(HDA))) {
    stop("Unable to automatically convert gene names to GRanges if the organism is not specified. See ?'organism<-'.")
  }
  
  genes <- rowData(HDA)$library_mutant
  genes <- unique(genes)
  genes <- genes[-which(is.na(genes))]
  
  genePositions <- downloadGenes(genes, organism(HDA))
  
  genePositions <- genePositions[match(rowData(HDA)$library_mutant, genePositions$ensembl_gene_id),]
  
  genePositions$strand[is.na(genePositions$strand)] <- "*"
  levels(genePositions$chromosome_name) <- c(levels(genePositions$chromosome_name), "NA")
  genePositions$chromosome_name[is.na(genePositions$chromosome_name)] <- "NA"
  genePositions$start_position[is.na(genePositions$start_position)] <- 0
  genePositions$end_position[is.na(genePositions$end_position)] <- 0
  
  genePositions <- GRanges(Rle(genePositions$chromosome_name), IRanges(genePositions$start_position, genePositions$end_position), Rle(genePositions$strand))
  
  rowRanges(HDA) <- genePositions
  
  return(HDA)
  
}

