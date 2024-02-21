#'@export
plotColonySizes <- function(HDA, what="auto", exclude="all", gene, ...) {
  
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  sizes <- extractAssay(HDA, what, exclude)[which(rowData(HDA)$replicate==gene), , drop=F]
  sizes <- melt(sizes)
  colnames(sizes) <- c("replicate", "sample", "value")
  
  for(i in 1:ncol(colData(HDA))) {
    sizes[[colnames(colData(HDA)[i])]] <- colData(HDA)[[i]][match(sizes$sample, rownames(colData(HDA)))]
  }
  
  args <- as.list(sapply(match.call()[-1], deparse))
  args <- args[-which(names(args) %in% c("HDA", "what", "exclude", "gene"))]
  
  if(length(args)>0) {
    sizes <- sizes[do.call("order", list(sizes[,unlist(args)])),]
  }
  
  sizes$sample <- as.character(sizes$sample)
  sizes$sample <- factor(sizes$sample, levels=unique(sizes$sample))
  
  args <- c(list(y="value", x="sample", group="replicate"), args)
  
  ggplot(sizes, do.call("aes_string", args)) + geom_bar(stat="identity", position="dodge")
  
}
