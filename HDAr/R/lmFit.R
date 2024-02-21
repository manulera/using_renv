#'@export
#'@title lmFit
setGeneric("lmFit", function(object, design, ndups=1, block=NULL, correlation, weights=NULL, method="ls", what="auto", exclude="all", ...) {standardGeneric("lmFit")})

setMethod("lmFit", c(object="HDA", design="missing"), function(object, design, ndups=1, block=NULL, correlation, weights=NULL, method="ls", what="auto", exclude="all", ...) {
  
  #Get the design matrix from the HDA object
  design <- designMatrix(object)
  
  #Need to explicitly specify arguments so that `design` is no longer treated as `missing` in next method
  if(missing(correlation)) {
    callNextMethod(object=object, design=design, ndups=ndups, block=block, weights=weights, method=method, what=what, exclude=exclude, ...)
  } else {
    callNextMethod(object=object, design=design, ndups=ndups, block=block, correlation=correlation, weights=weights, method=method, what=what, exclude=exclude, ...)
  }

})

setMethod("lmFit", c(object="HDA", design="ANY"), function(object, design, ndups=1, block=NULL, correlation, weights=NULL, method="ls", what="auto", exclude="all", ...) {
  
  #Check arguments
  what <- checkHDAwhat(object, what)
  exclude <- checkExclude(object, exclude)
  
  #Get names
  rownames <- rowData(object)$replicate
  
  #Get the matrix and replicates
  object <- extractAssay(object, what, exclude)
  rownames(object) <- rownames
  
  #Remomve footprints
  if(any(is.na(rownames))) {
    object <- object[-which(is.na(rownames)),]
  }
  
  #Remove genes which are all missing
  missingGenes <- apply(object, 1, function(i) {all(is.na(i))})
  if(any(missingGenes)) {
    object <- object[-which(missingGenes),]
  }
  
  #Order by replicate
  object <- object[order(rownames(object)),]
  
  #Average replicates
  #object <- limma::avereps(object)
  
  #Call next method (object="ANY")
  fit <- callNextMethod()
  
  #Set class of fit as HDALM
  new("HDALM", unclass(fit))
  
})

setMethod("lmFit", c(object="ANY", design="missing"), function(object, design, ndups=1, block=NULL, correlation, weights=NULL, method="ls", ...) {
  
  #Set design as NULL - default for limma version of lmFit
  design <- NULL
  
  #Need to explicitly specify arguments so that `design` is no longer treated as `missing` in next method
  if(missing(correlation)) {
    callNextMethod(object=object, design=design, ndups=ndups, block=block, weights=weights, method=method, what=what, exclude=exclude, ...)
  } else {
    callNextMethod(object=object, design=design, ndups=ndups, block=block, correlation=correlation, weights=weights, method=method, what=what, exclude=exclude, ...)
  }
  
})

setMethod("lmFit", c(object="ANY", design="ANY"), function(object, design, ndups=1, block=NULL, correlation, weights=NULL, method="ls", ...) {
  
  #Call lmFit from limma package
  if(missing(correlation)) {
    suppressWarnings(limma::lmFit(object, design, ndups, block, weights=weights, method=method, ...))
  } else {
    suppressWarnings(limma::lmFit(object, design, ndups, block, correlation, weights, method, ...))
  }
  
})

#'@export
topTable <- function(fit, coef) {
  table <- limma::topTable(fit, coef, nrow(coef(fit)), sort.by = "none")
  colnames(table)[1] <- "Colony.Size.Difference"
  table <- table[,-2]
  table <- table[,-5]
  return(table)
}

#'@export
volcanoplot <- function(fit, coef) {
  
  table <- topTable(fit, coef)
  table$logP <- -log10(table$adj.P.Val)
  
  ggplot(table, aes(Colony.Size.Difference, logP)) + geom_point() + xlab("Colony Size Difference") + ylab("-log10(P-value)")
  
}
