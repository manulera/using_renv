#A function to construct a normalisation surface for plate median normalisation
medianNormaliseColonySizes <- function(HDA, what, qcMat) {
  
  #Get the colony sizes to be normalised
  inputColonySizes <- assays(HDA)[[what]]
  
  #Mark excluded data as NA
  inputColonySizesQC <- inputColonySizes
  inputColonySizesQC[qcMat==F] <- NA
  
  #Calculate medians for each plate
  plateMedians <- stats::aggregate(inputColonySizesQC, by=list(plate=rowData(HDA)$plate), FUN=median, na.rm=T)
  
  #Match each colony to its respective plate median
  normalisationSurface <- sapply(colnames(inputColonySizes), function(i) {
    plateMedians[[i]][match(rowData(HDA)$plate, plateMedians$plate)]
  })
  
  #Return the normalisation surface
  return(normalisationSurface)
  
}

#A function to take xy coordinates and return smoothed y coordinates
rowcolLowess <- function(x, y, ...) {
  lowessOutput <- lowess(x, y, ...)
  return(lowessOutput$y)
}

#A function to construct a normalisation surface for rowcol normalisation
rowcolNormaliseColonySizes <- function(HDA, what, qcMat, lowess, span) {
  
  #Get the colony sizes to be normalised
  inputColonySizes <- assays(HDA)[[what]]
  
  #Mark excluded data as NA
  inputColonySizesQC <- inputColonySizes
  inputColonySizesQC[qcMat==F] <- NA
  
  #Calculate row and col medians for each plate
  rowMedians <- stats::aggregate(inputColonySizesQC, by=list(plate=rowData(HDA)$plate, row=rowData(HDA)$row), FUN=median, na.rm=T)
  colMedians <- stats::aggregate(inputColonySizesQC, by=list(plate=rowData(HDA)$plate, col=rowData(HDA)$col), FUN=median, na.rm=T)
  
  #If performing smoothing, then smooth the row and col medians
  if(lowess) {
    
    rowMedians <- plyr::ddply(
      rowMedians,
      plyr::.(plate),
      function(x) {
        smoothed <- apply(x[,colnames(HDA), drop=F], 2, rowcolLowess, x=x$row, f=span/nrow(x))
        rows <- matrix(x$row)
        colnames(rows) <- "row"
        smoothed <- cbind(rows, smoothed)
        return(smoothed)
      }
    )
    
    colMedians <- plyr::ddply(
      colMedians,
      plyr::.(plate),
      function(x) {
        smoothed <- apply(x[,colnames(HDA), drop=F], 2, rowcolLowess, x=x$col, f=span/nrow(x))
        cols <- matrix(x$col)
        colnames(cols) <- "col"
        smoothed <- cbind(cols, smoothed)
        return(smoothed)
      }
    )
    
  }
  
  #Calculate mean row and col medians for each plate
  rowMediansMeans <- stats::aggregate(rowMedians[,colnames(HDA), drop=F], by=list(plate=rowMedians$plate), FUN=mean)
  colMediansMeans <- stats::aggregate(colMedians[,colnames(HDA), drop=F], by=list(plate=colMedians$plate), FUN=mean)
  
  #Construct the normalisation surface by matching plate, row and col for each colony to the appropriate row/col median and row/col mean on medians for each plate
  #Use product of row/col median divided by product of mean of row/col medians for plate
  normalisationSurface <- sapply(colnames(inputColonySizes), function(i) {
    
    rows <- rowMedians[[i]][match(interaction(rowData(HDA)$plate, rowData(HDA)$row), interaction(rowMedians$plate, rowMedians$row))]
    rowMeans <- rowMediansMeans[[i]][match(rowData(HDA)$plate, rowMediansMeans$plate)]
    cols <- colMedians[[i]][match(interaction(rowData(HDA)$plate, rowData(HDA)$col), interaction(colMedians$plate, colMedians$col))]
    colMeans <- colMediansMeans[[i]][match(rowData(HDA)$plate, colMediansMeans$plate)]
    (rows*cols)/(rowMeans*colMeans)
    
  })
  
  #Return the normalisation surface
  return(normalisationSurface)
  
}

###PUBLIC DOMAIN CODE TAKEN FROM SGATOOLS BY CHARLIE BOONE'S GROUP###
###SEE http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3692131/ ###
###FUNCTIONS USED FOR SPATIAL EFFECT NORMALISATION###

# Spatial normalization: normalizes any gradient effect on the plate via median smoothing
# @param plate.data: SGA formatted data frame
# @param field.to.normalize: name of the column in the data to normalize
# @param ignore.ind: logical for any rows to be ignored 
# @return: vector of normalized values
spatialNormalization <- function(plate.data, field.to.normalize, ignore.ind, returnfilter=F) {  
  
  #num.rows = attr(plate.data, 'num.rows')
  #num.cols = attr(plate.data, 'num.cols')
  num.rows <- max(plate.data$Row)
  num.cols <- max(plate.data$Column)
  
  # Get gaussian/average filters
  gaussian.filt = fgaussian(7,2)
  average.filt = faverage(9)
  
  # Data to be normalized before ignored
  before.ignore = plate.data[[field.to.normalize]]
  
  # Data to be normalized after ignored (used in the analysis)
  after.ignore = before.ignore
  after.ignore[ignore.ind] = NA
  
  # Construct plate matrix
  plate.mat = matrix(NA, num.rows, num.cols)
  #rc.mat = as.matrix(plate.data[,1:2])
  rc.mat = as.matrix(cbind(plate.data$Row, plate.data$Column))
  plate.mat[rc.mat] = after.ignore
  
  # Fill NA with a placeholder (mean of all colonies) 
  t = plate.mat
  ind.na = which(is.na(t))
  t[ind.na] = mean(plate.mat, na.rm=TRUE)
  
  # Fill in NA with smoothed version of neighbors using gaussian blur
  filt.g = applyfilter(t, gaussian.filt)
  t[ind.na] = filt.g[ind.na]
  
  # Apply median/average filters
  filtered = medianfilter2d(t, 7, padding_type='replicate')
  filtered = applyfilter(filtered, average.filt, 'replicate')
  
  # Subtract the mean of the filtered data from the filtered data
  f = filtered / mean(filtered)
  
  #If the filter is to be returned, do so
  if(returnfilter==T) {
    return(f[rc.mat])
  }
  
  #Otherwise, return the normalised plate
  #If there are any correction factors which are 0, then apply no correction - to avoid dividing by 0 and generating NaN or Inf
  f <- ifelse(f==0, 1, f)
  
  # Subtract filtered - mean from  
  before.ignore = before.ignore / f[rc.mat]
  
  return(before.ignore)
  
}

# Filter functions used in spaital normalization: rewritten from matlab for R 

#Returns a gaussian filter matrix with dimensions x by x: equal to fspecial function in matlab
#Inputs:
#  x = dimensions (number of rows/cols) of the returned gaussian filter
#  sigma = standard deviation 
fgaussian <- function(x, sigma){
  x = gradientseq(x)
  mat = matrix(NA, length(x),length(x));
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      n1 = x[i];
      n2 = x[j];
      mat[i,j] = exp(-(n1^2+n2^2)/(2*sigma^2));
    }
  }
  mat = mat/sum(mat)
  return(mat)
}

#Helper function for fgaussian - given some value x, returns a gradient array begining with 0 on the inside and increasing outwards. Example: x = 7 returns [3,2,1,0,1,2,3] 
#Inputs:
#  x = number of elements in the returned array
gradientseq <- function(x){
  n = x;
  x = c(1:x)
  if(n%%2){
    rhs = x[1:floor(length(x)/2)];
    lhs = rev(rhs);
    return(c(lhs,0,rhs))
  }else{
    rhs = x[1:floor(length(x)/2)] - 0.5;
    lhs = rev(rhs);
    return(c(lhs,rhs))
  }
}

#Average filter
faverage <- function(size){
  x = 1/(size*size)
  ret = matrix(rep(x, size*size), size,size)
  return(ret)
}

#Applies a filter to a matrix: see imfilter (matlab) with replicate option
#Inputs:
#  mat = matrix to which the filter is applied
#   filter = a square matrix filter to be applied to the matrix 
applyfilter <- function(mat, filter, padding_type = 'zeros'){
  mat2 = mat
  fs = dim(filter);
  if(fs[1] != fs[2])
    stop('Filter must be a square matrix')
  if(fs[1] %% 2 == 0)
    stop('Filter dimensions must be odd')
  if(fs[1] == 1)
    stop('Filter dimensions must be greater than one')
  
  x = fs[1];
  a = (x-1)/2;
  
  s = dim(mat2)
  r = matrix(0, s[1], s[2])
  
  start = 1+a;
  end_1 = s[1]+a;
  end_2 = s[2]+a;
  
  mat2 = padmatrix(mat, a, padding_type)
  
  for(i in start:end_1){
    for(j in start:end_2){
      temp = mat2[(i-a):(i+a), (j-a):(j+a)] * filter;
      r[(i-a),(j-a)] = sum(temp)
    }
  }
  return(r)
}

#Applies a filter to a matrix: see imfilter (matlab) with replicate option
#Inputs:
#  mat = matrix to which the filter is applied
#   dim = number of rows/cols of window
medianfilter2d <- function(mat, dim, padding_type = 'zeros'){
  mat2 = mat
  fs = c()
  fs[1] = dim
  fs[2] = dim
  
  if(fs[1] != fs[2])
    stop('Filter must be a square matrix')
  if(fs[1] %% 2 == 0)
    stop('Filter dimensions must be odd')
  if(fs[1] == 1)
    stop('Filter dimensions must be greater than one')
  
  x = fs[1];
  a = (x-1)/2;
  
  s = dim(mat2)
  r = matrix(0, s[1], s[2])
  
  start = 1+a;
  end_1 = s[1]+a;
  end_2 = s[2]+a;
  
  mat2 = padmatrix(mat, a, padding_type)
  
  for(i in start:end_1){
    for(j in start:end_2){
      temp = mat2[(i-a):(i+a), (j-a):(j+a)];
      r[(i-a),(j-a)] = median(temp)
    }
  }
  return(r)
}

#Adds a padding to some matrix mat such that the padding is equal to the value of the nearest cell
#Inputs:
#	mat = matrix to which the padding is added
#	lvl = number of levels (rows/columns) of padding to be added
#	padding = type of padding on the matrix, zero will put zeros as borders, replicate will put the value of the nearest cell
padmatrix <- function(mat, lvl, padding){
  s = dim(mat);
  row_up = mat[1,]
  row_down = mat[s[1],]
  
  if(padding == 'zeros'){
    row_up = rep(0, length(row_up))
    row_down = rep(0, length(row_down))
  }
  #Add upper replicates
  ret = t(matrix(rep(row_up, lvl), length(as.vector(row_up))))
  #Add matrix itself
  ret = rbind(ret, mat)
  #Add lower replicates
  ret = rbind(ret, t(matrix(rep(row_down, lvl), length(as.vector(row_down)))))
  
  #Add columns
  s = dim(ret);
  col_left = ret[,1]
  col_right = ret[,s[2]]
  
  if(padding == 'zeros'){
    col_left = rep(0, length(col_left))
    col_right = rep(0, length(col_right))
  }
  
  #Add left columns
  ret2 = matrix(rep(col_left, lvl), length(as.vector(col_left)))
  #Add matrix itself
  ret2 = cbind(ret2, ret)
  #Add right columns
  ret2 = cbind(ret2, matrix(rep(col_right, lvl), length(as.vector(col_right))))
  
  #return 
  return(ret2)
}

###END OF FUNCTIONS TAKEN FROM SGATOOLS###

spatialNormaliseColonySizes <- function(HDA, what, qcMat) {
  
  sapply(colnames(HDA), function(i) {
    
    do.call("c", lapply(unique(rowData(HDA)$plate), function(j) {
      
      assaySub <- assays(HDA)[[what]][rowData(HDA)$plate==j, i, drop=F]
      assaySub <- cbind(data.frame(DataFrame(rowData(HDA))[rowData(HDA)$plate==j, c("row", "col")]), assaySub)
      colnames(assaySub)[1:2] <- c("Row", "Column")
      
      ignore.ind <- qcMat[rowData(HDA)$plate==j, i]==F
      
      return(spatialNormalization(assaySub, i, ignore.ind, T))
      
    }))
    
  })
  
}

#Create a list containing all of the normalisation methods currently available
normalisationMethods <- list(
  median = medianNormaliseColonySizes,
  rowcol = rowcolNormaliseColonySizes,
  spatial = spatialNormaliseColonySizes
)

#A function to construct a normalisation surface based on multiple normalisations
constructNormalisationSurface <- function(HDA, what, exclude, method=c("median", "rowcol", "spatial"), lowess=T, span=5) {
  
  #Get the excluded data
  qcMat <- excludedData(HDA, exclude)
  
  #Create a list to store the sequential normalisation surfaces which are constructed
  normalisationSurface <- list()
  
  #For each normalisation
  for(i in 1:length(method)) {
    
    #If it is the first normalisation, then normalise the dataset indicated by the user (normally "RS")
    #If subsequent, the normalise the output from the previous normalisation
    if(i > 1) {
      what <- "NS"
    }
    
    #Get all the arguments
    args <- list(HDA, what, qcMat)
    
    #If performing a rowcol normalisation, add the lowess and span arguments
    if(method[i] == "rowcol") {
      args <- c(args, list(lowess=lowess, span=span))
    }
    
    #Construct the appropriate normalisation surface and add it to the list of normalisation surfaces
    if(method[i] %in% names(normalisationMethods)==F) {
      stop(paste0("'", method[i], "' is not one of the currently available normalisation methods."))
    }
    normalisationSurface <- c(normalisationSurface, list(do.call(normalisationMethods[[method[i]]], args)))
    
    #Apply the current normalisation to the data
    #These normalised values will not be returned, by are required for subsequent iterations of the loop
    #This ensures that subsequent normalisations will be performed on top of the previously normalised data, not on the raw input indicated by the user (usually RS)
    assays(HDA)$NS <- assays(HDA)[[what]]/normalisationSurface[[i]]
    
  }
  
  #Combine all normalisation surfaces together into a single surface
  normalisationSurfaceAll <- normalisationSurface[[1]]
  if(length(normalisationSurface) > 1) {
    for(i in 2:length(normalisationSurface)) {
      normalisationSurfaceAll <- matrixcalc::hadamard.prod(normalisationSurfaceAll, normalisationSurface[[i]])
    }
  }
  
  #Return the normalisation surfaces
  return(normalisationSurfaceAll)
  
}

#'@export
normaliseColonySizes <- function(HDA, what="auto", exclude="all", method=c("median", "rowcol", "spatial"), lowess=T, span=5) {
  
  #Check arguments
  what <- checkHDAwhat(HDA, what)
  exclude <- checkExclude(HDA, exclude)
  
  #Construct the normalisation surface and normalise the data
  assays(HDA)$NS <- assays(HDA)[[what]]/constructNormalisationSurface(HDA, what, exclude, method, lowess, span)
  
  #Mark the normalised data as the active dataset
  activeDataset(HDA) <- "NS"
  
  #Add to the pipeline
  pipeline(HDA) <- c(pipeline(HDA), SimpleList(NS=SimpleList(method=method, what=what, exclude=exclude, lowess=lowess, span=span)))
  
  #Return the HDA object with the normalised colony sizes
  return(HDA)
  
}
