library(xlsx)
library(devtools)

setwd("~/Google_Drive/SGA Analysis Pipeline/SGA-Analysis-Pipeline/HDAr")

libraryFileList <- list()
libraryFileTransformations <- list()

#Create the Bioneer V5 library file
libraryFile <- read.xlsx("data-raw/Bioneer_collection_version_5.xlsx", 1, rowIndex = 1:3421, colIndex = 1:2)
libraryFile <- libraryFile[,c(2,1)]
colnames(libraryFile) <- c("library_mutant", "position")
libraryFile$library_mutant <- as.character(libraryFile$library_mutant)
libraryFile$replicate <- libraryFile$library_mutant
libraryFile <- libraryFile[,c(1,3,2)]
libraryFile$position <- as.character(libraryFile$position)
libraryFile$position <- gsub("V5-", "", libraryFile$position)
libraryFile <- cbind(libraryFile, matrix(unlist(strsplit(libraryFile$position, "-")), nrow=nrow(libraryFile), byrow=T))
libraryFile$position <- NULL
colnames(libraryFile)[3:4] <- c("plate", "position")
libraryFile$plate <- as.numeric(sub("P", "", libraryFile$plate))
libraryFile$position <- as.numeric(libraryFile$position)
plateFormat <- matrix(1:96, ncol=12, byrow=T)
libraryFile <- cbind(libraryFile, matrix(unlist(lapply(libraryFile$position, function(x) {which(plateFormat==x, arr.ind=T)})), nrow=nrow(libraryFile), byrow=T))
libraryFile$position <- NULL
colnames(libraryFile)[4:5] <- c("row", "col")
newLibraryFile <- expand.grid(plate=1:max(libraryFile$plate), row=1:max(libraryFile$row), col=1:max(libraryFile$col))
newLibraryFile <- merge(newLibraryFile, libraryFile, all=T)
newLibraryFile <- newLibraryFile[order(newLibraryFile$plate, newLibraryFile$row, newLibraryFile$col),]
newLibraryFile <- newLibraryFile[, c("library_mutant", "replicate", "plate", "row", "col")]

name <- "Fission_Yeast-Bioneer_V5-96_Well_Format"

libraryFileList[[name]] <- newLibraryFile
libraryFileTransformations[[name]] <- list(list("condense", c(2,2)), list("replicate", c(2,2)))

#Create the Bahler Lab ncRNA Library File
posInfo <- expand.grid(plate=1:8, row=1:8, col=1:12)
posInfo <- posInfo[order(posInfo$plate, posInfo$row, posInfo$col),]
libraryFile <- posInfo
libraryFile$library_mutant <- read.table("data-raw/Bahler_ncRNA_Library.txt", stringsAsFactors = F)[,1]
libraryFile <- libraryFile[,c(4, 1:3)]
libraryFile$replicate <- libraryFile$library_mutant
libraryFile <- libraryFile[,c(1, 5, 2:4)]

name <- "Fission_Yeast-Bahler_Lab_ncRNA_Deletion_Collection-96_Well_Format"

libraryFileList[[name]] <- libraryFile
libraryFileTransformations[[name]] <- list(list("condense", c(2,2)))

#Object to map the index of a library file to the source file
librarySourceFile <- do.call("c", lapply(1:length(libraryFileTransformations), function(i) {rep(i, 1+length(libraryFileTransformations[[i]]))}))

#Simple conversion of numeric to character describing how many replicates
nReplicatesVector <- c("singlicate", "duplicate", "triplicate", "quadruplicate")

use_data(libraryFileList, libraryFileTransformations, librarySourceFile, nReplicatesVector, overwrite = T, internal = T)
