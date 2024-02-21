excludedData <- function(HDA, exclude) {
  
  if(is.null(exclude)==F) {
    
    if("footprints" %in% exclude) {
      
      footprints <- T
      
      exclude <- exclude[-which(exclude == "footprints")]
      
    } else {
      
      footprints <- F
      
    }
    
  } else {
    
    footprints <- F
    
  }
  
  exclusions <- list()
  
  if(footprints) {
    exclusions[[1]] <- matrix(rep(is.na(rowData(HDA)$replicate), ncol(HDA)), ncol=ncol(HDA))
  }
  
  if(is.null(exclude)==F) {
    
    exclusions <- c(exclusions, lapply(exclude, function(x) {assays(HDA)[[x]]==F}))
    
  }
  
  if(length(exclusions)>0) {
    
    exclusionsAll <- exclusions[[1]]
    
    if(length(exclusions)>1) {
      
      for(i in 2:length(exclusions)) {
        
        exclusionsAll <- exclusionsAll | exclusions[[i]]
        
      }
      
    }
    
    exclusionsAll <- !exclusionsAll
    
  } else {
    
    exclusionsAll <- matrix(T, nrow(HDA), ncol(HDA))
    
  }
  
  exclusionsAll[is.na(exclusionsAll)] <- T
  
  return(exclusionsAll)
  
}