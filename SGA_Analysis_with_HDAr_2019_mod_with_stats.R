#############################################
###Analysis of S. pombe SGA data with gene deletions 
##library with St John's refurbished custard.
#####################################################
###Date: 05.04.2019
##By: Shajahan Anver
###################################################
###Sample pipelien @
## /Users/shajahananver/Bioinfo_Soft/Custard_Refurbsihed/Sample_Pipeline_Simple.R
###Look at the detailed .r file @
# /Users/shajahananver/Bioinfo_Soft/Custard_Refurbsihed/SGA_Analysis_with_HDAr.R
#######################################################


######################################################
####The scanned plate images need to be processed with
##gitter first to get the colony sizes
##.jpg format images are good
###put all the SGA image files in to one DIR and
##seperate the control and the query images in to
##their own subdir
##e.g. Main Dir ="SGA", subdir "Ade6", "aal1D"
##each has 9 images if you use Bioneer deletion library
##in 384 well-format
####################################################

###Install gitter
#install.packages ("gitter", dep=T) 

###call the gitter library now
library (gitter)

##Install the gitter dependency which doesn't get 
## installed with dep=T
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")

###call the gitter library again
library (gitter)

###Set the working directory to the .jpg image files
###gitter will save all the output files (2 per each imaage). 
## 1. .dat (tab delimited .txt file) and
## 2. a gridded image file.
###check the gridded image file for propper gridding
######################################################
###change the wd to the folder with images
setwd(".")

###Batch analysis for 384 well plate.format
###define the file path
f=file.path (getwd())
#f=file.path ("~/Research_Bahler_Lab/ncRNA/aal1_aal2_Expts/aal1_aal2_SGA/SGA_aal1oxNAT/24.05.2018/20180524_aal1oxNAT_Ade6_YES_24h_for_Analysis/Ade6")
gitter.batch (f, plate.format = c(16, 24), verbose = "l", inverse=T)
####need "inverse=T" as the scanned images are inversed, 
## 96_well_plates=8rows*12Cols

###to process individual images
gitter(image.file = file.choose(), plate.format = c(16, 24), remove.noise = F, autorotate = F, inverse = T, verbose = "l", contrast = NULL, fast = NULL, plot = F, grid.save = getwd(), dat.save = getwd(), .is.ref = F, .params = NULL)

###Files are renamed from .jpg.dat to .jpg.txt to 
## match the pattern (see below) but could also change 
## the pattern to match your file names.
##The .jpg.dat files must be seperated by the mutant
##in to their own subdirectory (if you haven't 
## done it before) with in the same main dir
### e.g. Analyzed SGA of aal1∆, with ade6∆ control then
###I have two subdirectories "aal1D" and "Ade6" under
### the directory "aal1D_SGA_03052019"

############################################################
#############################################
####set the working dir to the DIR with sample subdirectories
### "aal1D" and "Ade6" etc
setwd("~/Research_Bahler_Lab/ncRNA/aal1_aal2_Expts/aal1_aal2_SGA/SGA_aal1oxNAT/24.05.2018/20180524_aal1oxNAT_Ade6_EMM_30h_for_Analysis")



####################################################
###SGA analysis with St John's refurbished custard
##############################################
###Install refurbished custard fron the local repository
###Dependencies won't be installed if you are using a local
###repository. Deps needed to be installed first.

# install.packages("gitter", dep=T)
install.packages("matrixcalc", dep=T)
install.packages("biomartr", dep=T)
install.packages("cobs", dep=T)

###Now install the package
## Manu -> Add package here
install.packages ("", repos = NULL, type="source", dep=T)

##load the library
library (HDAr)
setwd("~/Research_Bahler_Lab/ncRNA/aal1_aal2_Expts/aal1_aal2_SGA/aal1_aal2_Mutant_SGA_YES_3%Glucose/SGA1.3.5_for_combined_Analysis")
#myDirectory <- "~/Research_Bahler_Lab/ncRNA/aal1_aal2_Expts/aal1_aal2_SGA/SGA_aal1oxNAT/24.05.2018/20180524_aal1oxNAT_Ade6_YES_24h_for_Analysis"
myDirectory <- getwd()

#Show all files
list.files(myDirectory, recursive = T)

#Create a FileExtractor to identify files with 
## the extension '.jpg.txt', where the files for 
#each sample are contained within their own folder.
myFileExtractor <- FileExtractor(regularExpression="(.*)/(.*)\\.jpg\\.txt")

#Create a SimpleList containing all files which are matched by myFileExtractor in the specified directory.
myFileList <- createFileList(myDirectory, myFileExtractor)

######################################################
###specify the library...check the available libraries
##############################
#listLibraryFiles()
myHDA <- HDA(myFileList, myDirectory, "gitter", listLibraryFiles()[2])

#Set control plate
batch(myHDA) <- rep(1:3, 2)
control(myHDA) <- c("Ade6_1", "Ade6_2", "Ade6_3")

#Plot the exclusions
plotSmallColonyExclusion(myHDA)

#Exclude small colonies
myHDA <- excludeSmallColonies(myHDA)

#Perform median and rowcol normalisations
plotPlate(myHDA, sample="Ade6_1", plate=9)
plotNormalisationSurface(myHDA, sample="Ade6_1", plate=9)
myHDA <- normaliseColonySizes(myHDA)
###normalizaes all-samples/all-plates

###If you want to check any plate, e.g. plate=9
#plotPlate(myHDA, sample="Ade6", plate=9)
#plotPlate(myHDA, sample="aal1OE", plate=9)

#Mark organism
organism(myHDA) <- "Schizosaccharomyces pombe"

#Add mutant positional information
myHDA <- addMutantPositionalInformation(myHDA)

#Add colRanges ... define your genes, Ade6 and gene of interest
##for linkage mapping
colRanges(myHDA) <- rep(c("SPBC2G5.02c", "SPCC1322.13"), each=3)

#View the linakge
plotLinkage(myHDA, sample="Ade6_1")###looks right
plotLinkedLociExclusion(myHDA, sample="Ade6_1", method="manual", distance=500000)
plotLinkage(myHDA, sample="aal1D_1") ##looks right
# + ggplot2::ylim(c(0,2))
plotLinkedLociExclusion(myHDA, sample="aal1D_1", method="manual", distance=500000) 

#Perform linkage exclusion
myHDA <- excludeLinkedLoci(myHDA, method="manual", distance = 500000)

#Set up models
library(limma)
treatment <- rep(c("Aal1", "Ade6"), each=3)
batch <- rep(1:3, 2)
design <- model.matrix(~0+treatment+batch)

#Run fits
sizes <- log2(extractAssay(myHDA, "NS", c("footprints", "SC", "L")))
fit <- lmFit(sizes, design)
fit <- contrasts.fit(fit, makeContrasts(treatmentAal1 - treatmentAde6, levels = design))
fit <- eBayes(fit)
hits <- topTable(fit, number = Inf, sort.by = "logFC")
plot(hits$logFC, -log10(hits$adj.P.Val), cex=0.1)

NS <- SimplifyHDA(myHDA)

Int <- log2(apply(NS[, c("aal1D_1", "aal1D_2", "aal1D_3")], 1, mean, na.rm=T)/apply(NS[, c("Ade6_1", "Ade6_2", "Ade6_3")], 1, mean, na.rm=T))
PVal <- sapply(1:nrow(NS), function(i) {
  tryCatch({
    t.test(as.vector(as.data.frame(NS[i, c("aal1D_1", "aal1D_2", "aal1D_3")])), as.vector(as.data.frame(NS[i, c("Ade6_1", "Ade6_2", "Ade6_3")])))$p.value
  }, error=function(e) {
    return(NA)
  })
})

plot(Int, -log10(PVal), cex=0.1)

#Create data output
MyDF <- SimplifyHDA(myHDA) ###get the object as a dataframe
MyDF <- as.data.frame(MyDF)
head (MyDF)
colnames (MyDF)

###Calculate Interactions
MyInt <- calculateInteractions(myHDA) 
head (MyInt)
MyInt <- as.data.frame(MyInt)

##merge to write everything to a file 
###################################
myData1 <- merge (MyDF, MyInt, by.x="library_mutant", 
                      by.y="row.names")
head (myData1)
summary(myData1)
write.table (myData1, 
             file="aal1D_SGA6_interactions_new_custard.txt", 
             sep="\t", quote=F, row.names=F)


